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
    algorithms::AlgorithmControlMap const & others
) const {
    return boost::make_shared<FitComboAlgorithm>(*this, boost::ref(schema), others);
}

//------------ FitComboModel ------------------------------------------------------------------------------

FitComboModel::FitComboModel(FitComboControl const & ctrl) :
    components(ndarray::allocate(ctrl.componentNames.size())),
    flux(std::numeric_limits<double>::quiet_NaN()), fluxErr(std::numeric_limits<double>::quiet_NaN()),
    chisq(std::numeric_limits<double>::quiet_NaN()), psfFactor(1.0)
{}

FitComboModel::FitComboModel(FitComboModel const & other) :
    components(ndarray::copy(other.components)), flux(other.flux), fluxErr(other.fluxErr),
    chisq(other.chisq), psfFactor(other.psfFactor)
{}

FitComboModel & FitComboModel::operator=(FitComboModel const & other) {
    if (&other != this) {
        components = ndarray::copy(other.components);
        flux = other.flux;
        fluxErr = other.fluxErr;
        chisq = other.chisq;
        psfFactor = other.psfFactor;
    }
    return *this;
}

//------------ FitComboAlgorithm --------------------------------------------------------------------------

FitComboAlgorithm::FitComboAlgorithm(
    FitComboControl const & ctrl,
    afw::table::Schema & schema,
    algorithms::AlgorithmControlMap const & others
) :
    algorithms::Algorithm(ctrl),
    _fluxKeys(
        afw::table::addFluxFields(
            schema, ctrl.name + ".flux",
            "combined flux of the linear combination of fixed-profile models"
        )
    ),
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
    if (ctrl.scaleByPsfFit) {
        _psfFactorKey = schema.addField<float>(
            ctrl.name + ".psf.factor",
            "PSF flux correction factor; multiply flux by this get uncorrected value"
        );  
    }
    algorithms::AlgorithmControlMap::const_iterator i = others.find(ctrl.psfName);
    if (i == others.end()) {
        throw LSST_EXCEPT(
            pex::exceptions::LogicErrorException,
            (boost::format("FitPsf with name '%s' not found; needed by FitCombo.") % ctrl.psfName).str()
        );
    }
    _psfCtrl = boost::dynamic_pointer_cast<FitPsfControl const>(i->second);
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
        algorithms::AlgorithmControlMap::const_iterator i = others.find(*nameIter);
        if (i == others.end()) {
            throw LSST_EXCEPT(
                pex::exceptions::LogicErrorException,
                (boost::format("FitProfile with name '%s' not found; needed by FitCombo.")
                 % (*nameIter)).str()
            );
        }
        _componentCtrl.push_back(boost::dynamic_pointer_cast<FitProfileControl const>(i->second));
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
    afw::image::Exposure<PixelT> const & image,
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
        return ModelInputHandler(image.getMaskedImage(), center,
                                 boundsEllipses, footprint, ctrl.growFootprint, 
                                 badPixelMask, ctrl.usePixelWeights);
    } else {
        return ModelInputHandler(image.getMaskedImage(), center, footprint, ctrl.growFootprint, 
                                 badPixelMask, ctrl.usePixelWeights);
    }
}

FitComboModel FitComboAlgorithm::apply(
    FitComboControl const & ctrl,
    FitPsfModel const & psfModel,
    std::vector<FitProfileModel> const & components,
    ModelInputHandler const & inputs,
    bool fitPsf
) {
    if (components.size() != 2u) {
        throw LSST_EXCEPT(
            pex::exceptions::LogicErrorException, "Only 2-component combo model is current implemented"
        );
    }
    FitComboModel model(ctrl);
    typedef shapelet::MultiShapeletFunction MSF;
    shapelet::ModelBuilder builder(inputs.getX(), inputs.getY());
    ndarray::Array<double,2,2> matrixT = ndarray::allocate(components.size(), inputs.getSize());
    ndarray::Array<double,2,-2> matrix(matrixT.transpose());
    matrixT.deep() = 0.0;
    for (int n = 0; n < matrixT.getSize<0>(); ++n) {
        MSF msf = components[n].asMultiShapelet(afw::geom::Point2D(), fitPsf)
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
    afw::math::LeastSquares lstsq = afw::math::LeastSquares::fromDesignMatrix(matrix, inputs.getData());
    if (lstsq.getSolution()[0] < 0.0) {
        if (lstsq.getSolution()[1] < 0.0) {
            throw LSST_EXCEPT(pex::exceptions::RuntimeErrorException, "measured negative flux");
        }
        model.components[0] = 0.0;
        model.components[1] = 1.0;
        if (fitPsf) {
            model.flux = components[1].psfFactor;
        } else {
            // remove PSF fit correction from single-profile result
            model.flux = components[1].flux * components[1].psfFactor;
            model.fluxErr = components[1].fluxErr * components[1].psfFactor;
        }
    } else if (lstsq.getSolution()[1] < 0.0) {
        model.components[0] = 1.0;
        model.components[1] = 0.0;
        if (fitPsf) {
            model.flux = components[0].psfFactor;
        } else {
            // remove PSF fit correction from single-profile result
            model.flux = components[0].flux * components[0].psfFactor;
            model.fluxErr = components[0].fluxErr * components[0].psfFactor;
        }
    } else {
        model.flux = lstsq.getSolution().asEigen().sum();
        model.fluxErr = std::sqrt(
            lstsq.getSolution().asEigen().dot(
                lstsq.getCovariance().asEigen() * lstsq.getSolution().asEigen()
            )
        );
        model.components.deep() = lstsq.getSolution();
        model.components.asEigen() /= model.flux;
    }
    model.chisq =
        (matrix.asEigen() * lstsq.getSolution().asEigen() - inputs.getData().asEigen()).squaredNorm()
        / (matrix.getSize<0>() - matrix.getSize<1>());
    return model;
}

template <typename PixelT>
FitComboModel FitComboAlgorithm::apply(
    FitComboControl const & ctrl,
    FitPsfModel const & psfModel,
    std::vector<FitProfileModel> const & components,
    afw::detection::Footprint const & footprint,
    afw::image::Exposure<PixelT> const & image,
    afw::geom::Point2D const & center
) {
    ModelInputHandler inputs = adjustInputs(
        ctrl, psfModel, components, footprint, image, center
    );
    FitComboModel model = apply(ctrl, psfModel, components, inputs);
    if (ctrl.scaleByPsfFit) {
        CONST_PTR(afw::detection::Psf) psf = image.getPsf();
        PTR(afw::image::Image<afw::math::Kernel::Pixel>) psfImage = psf->computeImage(center);
        double s = psfImage->getArray().asEigen().sum();
        psfImage->getArray().asEigen() /= s;
        ModelInputHandler psfInputs(*psfImage, center, psfImage->getBBox(afw::image::PARENT));
        FitComboModel psfProfileModel = apply(ctrl, psfModel, components, psfInputs, true);
        model.psfFactor = psfProfileModel.flux;
        model.flux /= model.psfFactor;
        model.fluxErr /= model.psfFactor;
    }
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
        if (components.back().flagFailed) {
            return; // Don't bother trying linear components if one of the inputs failed.
                    // This saves us from slowdown caused by bogus huge initial ellipses/
        }
        assert(lsst::utils::isfinite(components.back().ellipse.getArea()));
    }
    FitComboModel model = apply(
        getControl(), psfModel, components,
        *source.getFootprint(),
        exposure, center
    );
    source.set(_componentsKey, model.components);
    source.set(_fluxKeys.meas, model.flux);
    source.set(_fluxKeys.err, model.fluxErr);
    source.set(_fluxKeys.flag, false);
    if (getControl().scaleByPsfFit) {
        source.set(_psfFactorKey, model.psfFactor);
    }
    source.set(_chisqKey, model.chisq);
}


LSST_MEAS_ALGORITHM_PRIVATE_IMPLEMENTATION(FitComboAlgorithm);

}}}} // namespace lsst::meas::extensions::multiShapelet
