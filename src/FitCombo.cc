// -*- lsst-c++ -*-
/* 
 * LSST Data Management System
 * Copyright 2012 LSST Corporation.
 * 
 * This product includes software developed by the
 * LSST Project (http://www.lsst.org/).
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the LSST License Statement and 
 * the GNU General Public License along with this program.  If not, 
 * see <http://www.lsstcorp.org/LegalNotices/>.
 */

#include "lsst/utils/ieee.h"
#include "lsst/meas/extensions/multiShapelet/FitCombo.h"
#include "lsst/shapelet/ModelBuilder.h"
#include "lsst/afw/math/LeastSquares.h"
#include "lsst/afw/detection/FootprintArray.h"
#include "lsst/afw/detection/FootprintArray.cc"


namespace lsst { namespace meas { namespace extensions { namespace multiShapelet {

//------------ FitComboControl ----------------------------------------------------------------------------

PTR(algorithms::AlgorithmControl) FitComboControl::_clone() const {
    return boost::make_shared<FitComboControl>(*this);
}

PTR(algorithms::Algorithm) FitComboControl::_makeAlgorithm(
    afw::table::Schema & schema,
    PTR(daf::base::PropertyList) const & metadata,
    algorithms::AlgorithmMap const & others
) const {
    return boost::make_shared<FitComboAlgorithm>(*this, boost::ref(schema), others);
}

//------------ FitComboModel ------------------------------------------------------------------------------

FitComboModel::FitComboModel(FitComboControl const & ctrl) :
    components(ndarray::allocate(ctrl.componentNames.size())),
    flux(std::numeric_limits<double>::quiet_NaN()), fluxErr(std::numeric_limits<double>::quiet_NaN()),
    chisq(std::numeric_limits<double>::quiet_NaN())
{}

FitComboModel::FitComboModel(FitComboModel const & other) :
    components(ndarray::copy(other.components)), flux(other.flux), fluxErr(other.fluxErr),
    chisq(other.chisq)
{}

FitComboModel & FitComboModel::operator=(FitComboModel const & other) {
    if (&other != this) {
        components = ndarray::copy(other.components);
        flux = other.flux;
        fluxErr = other.fluxErr;
        chisq = other.chisq;
    }
    return *this;
}

//------------ FitComboAlgorithm --------------------------------------------------------------------------

FitComboAlgorithm::FitComboAlgorithm(
    FitComboControl const & ctrl,
    afw::table::Schema & schema,
    algorithms::AlgorithmMap const & others
) :
    algorithms::Algorithm(ctrl),
    _fluxKeys(
        afw::table::addFluxFields(
            schema, ctrl.name + ".flux",
            "combined flux of the linear combination of fixed-profile models"
        )
    ),
    _fluxCorrectionKeys(ctrl.name, schema),
    _componentsKey(
        schema.addField< afw::table::Array<float> >(
            ctrl.name + ".components",
            "relative fluxes of fixed-profile components, normalized to sum to one",
            ctrl.componentNames.size()
        )),
    _chisqKey(
        schema.addField<float>(
            ctrl.name + ".chisq",
            "reduced chi^2"
        ))
{
    algorithms::AlgorithmMap::const_iterator i = others.find(ctrl.psfName);
    if (i == others.end()) {
        throw LSST_EXCEPT(
            pex::exceptions::LogicErrorException,
            (boost::format("FitPsf with name '%s' not found; needed by FitCombo.") % ctrl.psfName).str()
        );
    }
    _psfCtrl = boost::dynamic_pointer_cast<FitPsfControl const>(i->second->getControl().clone());
    if (!_psfCtrl) {
        throw LSST_EXCEPT(
            pex::exceptions::LogicErrorException,
            (boost::format("Algorithm with name '%s' is not FitPsf.") % ctrl.psfName).str()
        );
    }
    for (
        std::vector<std::string>::const_iterator nameIter = ctrl.componentNames.begin();
        nameIter != ctrl.componentNames.end();
        ++nameIter
    ) {
        algorithms::AlgorithmMap::const_iterator i = others.find(*nameIter);
        if (i == others.end()) {
            throw LSST_EXCEPT(
                pex::exceptions::LogicErrorException,
                (boost::format("FitProfile with name '%s' not found; needed by FitCombo.")
                 % (*nameIter)).str()
            );
        }
        _componentCtrl.push_back(
            boost::dynamic_pointer_cast<FitProfileControl const>(i->second->getControl().clone())
        );
        if (!_componentCtrl.back()) {
            throw LSST_EXCEPT(
                pex::exceptions::LogicErrorException,
                (boost::format("Algorithm with name '%s' is not FitProfile.") % (*nameIter)).str()
            );
        }
    }
}

template <typename PixelT>
ModelInputHandler FitComboAlgorithm::adjustInputs(
    FitComboControl const & ctrl,
    FitPsfModel const & psfModel,
    std::vector<FitProfileModel> const & components,
    afw::detection::Footprint const & footprint,
    afw::image::MaskedImage<PixelT> const & image,
    afw::geom::Point2D const & center
) {
    afw::image::MaskPixel badPixelMask(0);
        badPixelMask = afw::image::Mask<>::getPlaneBitMask(ctrl.badMaskPlanes);
    if (ctrl.radiusInputFactor > 0.0) {
        std::vector<afw::geom::ellipses::Ellipse> boundsEllipses;
        for (std::size_t n = 0; n < components.size(); ++n) {
            boundsEllipses.push_back(afw::geom::ellipses::Ellipse(components[n].ellipse, center));
            boundsEllipses.back().getCore().scale(ctrl.radiusInputFactor);
        }
        return ModelInputHandler(image, center,
                                 boundsEllipses, footprint, ctrl.growFootprint, 
                                 badPixelMask, ctrl.usePixelWeights);
    } else {
        return ModelInputHandler(image, center, footprint, ctrl.growFootprint, 
                                 badPixelMask, ctrl.usePixelWeights);
    }
}

FitComboModel FitComboAlgorithm::apply(
    FitComboControl const & ctrl,
    FitPsfModel const & psfModel,
    std::vector<FitProfileModel> const & components,
    ModelInputHandler const & inputs
) {
    if (components.size() != 2u) {
        throw LSST_EXCEPT(
            pex::exceptions::LogicErrorException, "Only 2-component combo model is current implemented"
        );
    }
    FitComboModel model(ctrl);
    typedef shapelet::MultiShapeletFunction MSF;
    shapelet::ModelBuilder<double> builder(inputs.getX(), inputs.getY());
    ndarray::Array<double,2,2> matrixT = ndarray::allocate(components.size(), inputs.getSize());
    ndarray::Array<double,2,-2> matrix(matrixT.transpose());
    matrixT.deep() = 0.0;
    for (int n = 0; n < matrixT.getSize<0>(); ++n) {
        MSF msf = components[n].asMultiShapelet(afw::geom::Point2D())
            .convolve(psfModel.asMultiShapelet());
        msf.normalize();
        for (
            MSF::ElementList::const_iterator i = msf.getElements().begin();
            i != msf.getElements().end();
            ++i
        ) {
            builder.update(i->getEllipse().getCore());
            builder.addModelVector(i->getOrder(), i->getCoefficients(), matrixT[n]);
        }
        if (!inputs.getWeights().isEmpty()) {
            matrixT[n].asEigen<Eigen::ArrayXpr>() *= inputs.getWeights().asEigen<Eigen::ArrayXpr>();
        }
    }
    // We should really do constrained linear least squares to get the errors right, but this
    // produces the same result for the fluxes, and we don't have a constrained solver handy.
    afw::math::LeastSquares lstsq = afw::math::LeastSquares::fromDesignMatrix(matrix, inputs.getData());
    if (lstsq.getSolution()[0] < 0.0) {
        if (lstsq.getSolution()[1] < 0.0) {
            throw LSST_EXCEPT(pex::exceptions::RuntimeErrorException, "measured negative flux");
        }
        model.components[0] = 0.0;
        model.components[1] = 1.0;
        model.flux = components[1].flux;
        model.fluxErr = components[1].fluxErr;
    } else if (lstsq.getSolution()[1] < 0.0) {
        model.components[0] = 1.0;
        model.components[1] = 0.0;
        model.flux = components[0].flux;
        model.fluxErr = components[0].fluxErr;
    } else {
        model.flux = lstsq.getSolution().asEigen().sum();
        model.components.deep() = lstsq.getSolution();
        model.components.asEigen() /= model.flux;
#if 0 // don't know why this doesn't work; numbers are way too big
        model.fluxErr = std::sqrt(
            lstsq.getSolution().asEigen().dot(
                lstsq.getCovariance().asEigen() * lstsq.getSolution().asEigen()
            )
        );
#else // this is incorrect, but a good-enough workaround for now: weighted average in quadrature
        model.fluxErr = std::sqrt(components[0].fluxErr * components[0].fluxErr * model.components[0]
                                  + components[1].fluxErr * components[1].fluxErr * model.components[1]);
#endif
    }
    model.chisq =
        (matrix.asEigen() * lstsq.getSolution().asEigen() - inputs.getData().asEigen()).squaredNorm()
        / (matrix.getSize<0>() - matrix.getSize<1>());
    return model;
}

template <typename PixelT>
void FitComboAlgorithm::_apply(
    afw::table::SourceRecord & source,
    afw::image::Exposure<PixelT> const & exposure,
    afw::geom::Point2D const & center
) const {
    source.set(_fluxKeys.flag, true);
    if (!exposure.hasPsf()) {
        throw LSST_EXCEPT(
            pex::exceptions::LogicErrorException,
            "Cannot run FitComboAlgorithm without a PSF."
        );
    }
    FitPsfModel psfModel(*_psfCtrl, source);
    std::vector<FitProfileModel> components;
    for (std::size_t n = 0; n < _componentCtrl.size(); ++n) {
        components.push_back(FitProfileModel(*_componentCtrl[n], source));
        if (components.back().fluxFlag) {
            return; // Don't bother trying linear components if one of the inputs failed.
        }
        assert(lsst::utils::isfinite(components.back().ellipse.getArea()));
    }
    ModelInputHandler inputs = adjustInputs(
        getControl(), psfModel, components, *source.getFootprint(), exposure.getMaskedImage(), center);
    FitComboModel model = apply(getControl(), psfModel, components, inputs);

    source.set(_componentsKey, model.components);
    source.set(_fluxKeys.meas, model.flux);
    source.set(_fluxKeys.err, model.fluxErr);
    source.set(_fluxKeys.flag, false);
    source.set(_chisqKey, model.chisq);

    source.set(_fluxCorrectionKeys.psfFactorFlag, true);
    PTR(afw::image::Image<afw::math::Kernel::Pixel>) psfImage = exposure.getPsf()->computeImage(center);
    ModelInputHandler psfInputs(*psfImage, center, psfImage->getBBox(afw::image::PARENT));
    std::vector<FitProfileModel> psfComponents;
    for (std::size_t n = 0; n < _componentCtrl.size(); ++n) {
        psfComponents.push_back(FitProfileModel(*_componentCtrl[n], source, true));
        if (psfComponents.back().fluxFlag) {
            return; // Don't bother trying linear components if one of the inputs failed.
        }
        assert(lsst::utils::isfinite(psfComponents.back().ellipse.getArea()));
    }
    FitComboModel psfProfileModel = apply(getControl(), psfModel, psfComponents, psfInputs);
    source.set(_fluxCorrectionKeys.psfFactor, psfProfileModel.flux);
    source.set(_fluxCorrectionKeys.psfFactorFlag, false);
    
}


LSST_MEAS_ALGORITHM_PRIVATE_IMPLEMENTATION(FitComboAlgorithm);

}}}} // namespace lsst::meas::extensions::multiShapelet
