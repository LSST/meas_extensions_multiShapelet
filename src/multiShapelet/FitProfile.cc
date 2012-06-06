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

#include "lsst/meas/extensions/multiShapelet/FitProfile.h"
#include "lsst/meas/extensions/multiShapelet/MultiGaussianObjective.h"
#include "lsst/meas/extensions/multiShapelet/MultiGaussianRegistry.h"
#include "lsst/shapelet/ModelBuilder.h"
#include "lsst/afw/math/LeastSquares.h"
#include "lsst/afw/detection/FootprintArray.h"
#include "lsst/afw/detection/FootprintArray.cc"

namespace lsst { namespace meas { namespace extensions { namespace multiShapelet {

//------------ FitProfileControl ----------------------------------------------------------------------------

PTR(algorithms::AlgorithmControl) FitProfileControl::_clone() const {
    return boost::make_shared<FitProfileControl>(*this);
}

PTR(algorithms::Algorithm) FitProfileControl::_makeAlgorithm(
    afw::table::Schema & schema,
    PTR(daf::base::PropertyList) const & metadata,
    algorithms::AlgorithmControlMap const & others
) const {
    return boost::make_shared<FitProfileAlgorithm>(*this, boost::ref(schema), others);
}

//------------ FitProfileModel ------------------------------------------------------------------------------

FitProfileModel::FitProfileModel(
    FitProfileControl const & ctrl,
    double amplitude,
    ndarray::Array<double const,1,1> const & parameters
) :
    profile(ctrl.profile), flux(amplitude), fluxErr(0.0),
    ellipse(MultiGaussianObjective::EllipseCore(parameters[0], parameters[1], parameters[2])),
    chisq(std::numeric_limits<double>::quiet_NaN()), psfFactor(1.0),
    flagMaxIter(false), flagTinyStep(false), flagMinRadius(false), flagMinAxisRatio(false)
{}

FitProfileModel::FitProfileModel(
    FitProfileControl const & ctrl, afw::table::SourceRecord const & source
) :
    profile(ctrl.profile), flux(1.0), fluxErr(0.0), ellipse(),
    chisq(std::numeric_limits<double>::quiet_NaN()), psfFactor(1.0),
    flagMaxIter(false), flagTinyStep(false), flagMinRadius(false), flagMinAxisRatio(false)
{
    afw::table::SubSchema s = source.getSchema()[ctrl.name];
    flux = source.get(s.find< double >("flux").key);
    fluxErr = source.get(s.find< double >("flux.err").key);
    ellipse = source.get(s.find< afw::table::Moments<float> >("ellipse").key);
    chisq = source.get(s.find<float>("chisq").key);
    flagMaxIter = source.get(s.find<afw::table::Flag>("flags.maxiter").key);
    flagTinyStep = source.get(s.find<afw::table::Flag>("flags.tinystep").key);
    if (ctrl.scaleByPsfFit) {
        try {
            psfFactor = source.get(s.find<float>("psf.factor").key);
            psfEllipse = boost::make_shared<afw::geom::ellipses::Quadrupole>(
                source.get(s.find< afw::table::Moments<float> >("psf.ellipse").key)
            );
        } catch (pex::exceptions::NotFoundException &) {}
    }
}

FitProfileModel::FitProfileModel(FitProfileModel const & other) :
    profile(other.profile), flux(other.flux), fluxErr(other.fluxErr), ellipse(other.ellipse),
    chisq(other.chisq), psfFactor(other.psfFactor),
    flagMaxIter(other.flagMaxIter),
    flagTinyStep(other.flagTinyStep),
    flagMinRadius(other.flagMinRadius),
    flagMinAxisRatio(other.flagMinAxisRatio)
{
    if (other.psfEllipse) {
        psfEllipse = boost::make_shared<afw::geom::ellipses::Quadrupole>(*other.psfEllipse);
    }
}

FitProfileModel & FitProfileModel::operator=(FitProfileModel const & other) {
    if (&other != this) {
        profile = other.profile;
        flux = other.flux;
        fluxErr = other.fluxErr;
        ellipse = other.ellipse;
        chisq = other.chisq;
        psfFactor = other.psfFactor;
        flagMaxIter = other.flagMaxIter;
        flagTinyStep = other.flagTinyStep;
        flagMinRadius = other.flagMinRadius;
        flagMinAxisRatio = other.flagMinAxisRatio;
        if (other.psfEllipse) {
            psfEllipse = boost::make_shared<afw::geom::ellipses::Quadrupole>(*other.psfEllipse);
        } else {
            psfEllipse.reset();
        }
    }
    return *this;
}

shapelet::MultiShapeletFunction FitProfileModel::asMultiShapelet(
    afw::geom::Point2D const & center
) const {
    shapelet::MultiShapeletFunction::ElementList elements;
    MultiGaussian const & multiGaussian = MultiGaussianRegistry::lookup(profile);
    for (MultiGaussian::const_iterator i = multiGaussian.begin(); i != multiGaussian.end(); ++i) {
        afw::geom::ellipses::Ellipse fullEllipse(ellipse, center);
        elements.push_back(i->makeShapelet(fullEllipse));
        elements.back().getCoefficients().asEigen() *= flux;
    }
    return shapelet::MultiShapeletFunction(elements);
}

//------------ FitProfileAlgorithm --------------------------------------------------------------------------

FitProfileAlgorithm::FitProfileAlgorithm(
    FitProfileControl const & ctrl,
    afw::table::Schema & schema,
    algorithms::AlgorithmControlMap const & others
) :
    algorithms::Algorithm(ctrl),
    _fluxKeys(
        afw::table::addFluxFields(
            schema, ctrl.name + ".flux",
            "flux of multi-Gaussian approximation to the profile, integrated to infinite radius"
        )
    ),
    _ellipseKey(
        schema.addField< afw::table::Moments<float> >(
            ctrl.name + ".ellipse",
            "half-light radius ellipse"
        )),
    _chisqKey(
        schema.addField<float>(
            ctrl.name + ".chisq",
            "reduced chi^2"
        )),
    _flagMaxIterKey(
        schema.addField<afw::table::Flag>(
            ctrl.name + ".flags.maxiter",
            "set if the optimizer ran into the maximum number of iterations limit"
        )),
    _flagTinyStepKey(
        schema.addField<afw::table::Flag>(
            ctrl.name + ".flags.tinystep",
            "set if the optimizer step or trust region got so small no progress could be made"
        )),
    _flagMinRadiusKey(
        schema.addField<afw::table::Flag>(
            ctrl.name + ".flags.constraint.r",
            "set if the best-fit radius was the minimum allowed by the constraint (not a failure)"
        )),
    _flagMinAxisRatioKey(
        schema.addField<afw::table::Flag>(
            ctrl.name + ".flags.constraint.q",
            "set if the best-fit ellipticity was the maximum allowed by the constraint"
        )),
    _psfCtrl()
{
    if (ctrl.scaleByPsfFit) {
        _psfEllipseKey = schema.addField< afw::table::Moments<float> >(
            ctrl.name + ".psf.ellipse",
            "ellipse from fitting the profile model to the PSF model"
        );
        _psfFactorKey = schema.addField<float>(
            ctrl.name + ".psf.factor",
            "PSF flux correction factor; multiply flux by this get uncorrected value"
        );  
    }
    algorithms::AlgorithmControlMap::const_iterator i = others.find(ctrl.psfName);
    if (i == others.end()) {
        throw LSST_EXCEPT(
            pex::exceptions::LogicErrorException,
            (boost::format("FitPsf with name '%s' not found; needed by FitProfile.") % ctrl.psfName).str()
        );
    }
    _psfCtrl = boost::dynamic_pointer_cast<FitPsfControl const>(i->second);
    if (!_psfCtrl) {
        throw LSST_EXCEPT(
            pex::exceptions::LogicErrorException,
            (boost::format("Algorithm with name '%s' is not FitPsf.") % ctrl.psfName).str()
        );
    }
}

PTR(MultiGaussianObjective) FitProfileAlgorithm::makeObjective(
    FitProfileControl const & ctrl,
    FitPsfModel const & psfModel,
    ModelInputHandler const & inputs
) {
    return boost::make_shared<MultiGaussianObjective>(
        inputs, ctrl.getMultiGaussian(), psfModel.getMultiGaussian(), psfModel.ellipse,
        ctrl.minRadius, ctrl.minAxisRatio
    );
}

HybridOptimizer FitProfileAlgorithm::makeOptimizer(
    FitProfileControl const & ctrl,
    FitPsfModel const & psfModel,
    MultiGaussianObjective::EllipseCore const & ellipse,
    ModelInputHandler const & inputs
) {
    PTR(Objective) obj = makeObjective(ctrl, psfModel, inputs);
    ndarray::Array<double,1,1> initial = ndarray::allocate(obj->getParameterSize());
    ellipse.writeParameters(initial.getData());
    HybridOptimizerControl optCtrl; // TODO: nest this in FitProfileControl
    optCtrl.tau = 1E-2;
    optCtrl.useCholesky = true;
    optCtrl.gTol = 1E-4;
    return HybridOptimizer(obj, initial, optCtrl);
}

template <typename PixelT>
ModelInputHandler FitProfileAlgorithm::adjustInputs(
    FitProfileControl const & ctrl,
    FitPsfModel const & psfModel,
    afw::geom::ellipses::Quadrupole & shape,
    afw::detection::Footprint const & footprint,
    afw::image::Exposure<PixelT> const & image,
    afw::geom::Point2D const & center
) {
   
    MultiGaussianObjective::EllipseCore ellipse(shape);
    if (ctrl.deconvolveShape) {
        try {
            ellipse = ctrl.getMultiGaussian().deconvolve(
                shape, psfModel.ellipse, psfModel.getMultiGaussian()
            );
        } catch (pex::exceptions::InvalidParameterException &) {
            ellipse = psfModel.ellipse;
            ellipse.scale(ctrl.minInitialRadius);
        }
    }
    // We never want to start with an ellipse smaller than the PSF or an ellipticity
    // on the constraint, because we might never find our way out.
    std::pair<bool,bool> constrained = MultiGaussianObjective::constrainEllipse(
        ellipse, psfModel.ellipse.getTraceRadius() * ctrl.minInitialRadius, ctrl.minAxisRatio
    );
    if (constrained.first || constrained.second) {
        ellipse = psfModel.ellipse;
        ellipse.scale(ctrl.minInitialRadius);
    }
    shape = ellipse;
    afw::image::MaskPixel badPixelMask = afw::image::Mask<>::getPlaneBitMask(ctrl.badMaskPlanes);
    if (ctrl.radiusInputFactor > 0.0) {
        afw::geom::ellipses::Ellipse boundsEllipse(shape, center);
        boundsEllipse.getCore().scale(ctrl.radiusInputFactor);
        return ModelInputHandler(image.getMaskedImage(), boundsEllipse, footprint, ctrl.growFootprint, 
                                 badPixelMask, ctrl.usePixelWeights);
    } else {
        return ModelInputHandler(image.getMaskedImage(), center, footprint, ctrl.growFootprint, 
                                 badPixelMask, ctrl.usePixelWeights);
    }
}

void FitProfileAlgorithm::fitShapeletTerms(
    FitProfileControl const & ctrl,
    FitPsfModel const & psfModel,
    ModelInputHandler const & inputs,
    FitProfileModel & model
) {
    typedef shapelet::MultiShapeletFunction MSF; 
    MSF msf = model.asMultiShapelet().convolve(psfModel.asMultiShapelet());
    ndarray::Array<double,1,1> vector = ndarray::allocate(inputs.getSize());
    vector.deep() = 0.0;
    shapelet::ModelBuilder builder(inputs.getX(), inputs.getY());
    for (MSF::ElementList::const_iterator i = msf.getElements().begin(); i != msf.getElements().end(); ++i) {
        builder.addModelVector(i->getOrder(), i->getCoefficients(), vector);
    }
    if (!inputs.getWeights().isEmpty()) {
        vector.asEigen<Eigen::ArrayXpr>() *= inputs.getWeights().asEigen<Eigen::ArrayXpr>();
    }
    // the following is just linear least squares with one free parameter
    double variance = 1.0 / vector.asEigen().squaredNorm();
    model.flux = vector.asEigen().dot(inputs.getData().asEigen()) * variance;
    model.fluxErr = std::sqrt(variance);
    model.chisq = (model.flux * vector.asEigen() - inputs.getData().asEigen()).squaredNorm()
        / (inputs.getSize() - 4);
}

FitProfileModel FitProfileAlgorithm::apply(
    FitProfileControl const & ctrl,
    FitPsfModel const & psfModel,
    MultiGaussianObjective::EllipseCore const & inEllipse,
    ModelInputHandler const & inputs
) {
    HybridOptimizer opt = makeOptimizer(ctrl, psfModel, inEllipse, inputs);
    opt.run();
    Model model(
        ctrl, 
        boost::static_pointer_cast<MultiGaussianObjective const>(opt.getObjective())->getAmplitude(),
        opt.getParameters()
    );
    MultiGaussianObjective::EllipseCore ellipse = MultiGaussianObjective::readParameters(opt.getParameters());
    std::pair<bool,bool> constrained 
        = MultiGaussianObjective::constrainEllipse(ellipse, ctrl.minRadius, ctrl.minAxisRatio);
    model.flagMaxIter = opt.getState() & HybridOptimizer::FAILURE_MAXITER;
    model.flagTinyStep = (opt.getState() & HybridOptimizer::FAILURE_MINSTEP)
        || (opt.getState() & HybridOptimizer::FAILURE_MINTRUST);
    model.flagMinRadius = constrained.first;
    model.flagMinAxisRatio = constrained.second;
    if (ctrl.usePsfShapeletTerms) {
        fitShapeletTerms(ctrl, psfModel, inputs, model);
    } else {
        model.chisq = opt.getChiSq() / (inputs.getSize() - 4);
    }
    return model;
}

template <typename PixelT>
FitProfileModel FitProfileAlgorithm::apply(
    FitProfileControl const & ctrl,
    FitPsfModel const & psfModel,
    afw::geom::ellipses::Quadrupole & shape,
    afw::detection::Footprint const & footprint,
    afw::image::Exposure<PixelT> const & image,
    afw::geom::Point2D const & center
) {
    ModelInputHandler inputs = adjustInputs(
        ctrl, psfModel, shape, footprint, image, center
    );
    MultiGaussianObjective::EllipseCore ellipse(shape);
    FitProfileModel model = apply(ctrl, psfModel, ellipse, inputs);
    if (ctrl.scaleByPsfFit) {
        CONST_PTR(afw::detection::Psf) psf = image.getPsf();
        PTR(afw::image::Image<afw::math::Kernel::Pixel>) psfImage = psf->computeImage(center);
        double s = psfImage->getArray().asEigen().sum();
        psfImage->getArray().asEigen() /= s;
        ModelInputHandler psfInputs(*psfImage, center, psfImage->getBBox(afw::image::PARENT));
        MultiGaussianObjective::EllipseCore psfEllipse(psfModel.ellipse);
        psfEllipse.scale(ctrl.minInitialRadius);
        FitProfileModel psfProfileModel = apply(ctrl, psfModel, psfEllipse, psfInputs);
        model.psfFactor = psfProfileModel.flux;
        model.flux /= model.psfFactor;
        model.fluxErr /= model.psfFactor;
        model.psfEllipse = boost::make_shared<afw::geom::ellipses::Quadrupole>(psfProfileModel.ellipse);
    }
    return model;
}

template <typename PixelT>
void FitProfileAlgorithm::_apply(
    afw::table::SourceRecord & source,
    afw::image::Exposure<PixelT> const & exposure,
    afw::geom::Point2D const & center
) const {
    source.set(_fluxKeys.flag, true);
    if (!exposure.hasPsf()) {
        throw LSST_EXCEPT(
            pex::exceptions::LogicErrorException,
            "Cannot run FitProfileAlgorithm without a PSF."
        );
    }
    FitPsfModel psfModel(*_psfCtrl, source);
    afw::geom::ellipses::Quadrupole shape = source.getShape();
    if (source.getShapeFlag()) {
        shape = psfModel.ellipse;
    }
    FitProfileModel model = apply(
        getControl(), psfModel,
        shape, *source.getFootprint(),
        exposure, center
    );
    source.set(_fluxKeys.meas, model.flux);
    source.set(_fluxKeys.err, model.fluxErr);
    source.set(_fluxKeys.flag, false);
    source.set(_ellipseKey, model.ellipse);
    if (getControl().scaleByPsfFit) {
        source.set(_psfEllipseKey, *model.psfEllipse);
        source.set(_psfFactorKey, model.psfFactor);
    }
    source.set(_chisqKey, model.chisq);
    source.set(_flagMaxIterKey, model.flagMaxIter);
    source.set(_flagTinyStepKey, model.flagTinyStep);
    source.set(_flagMinRadiusKey, model.flagMinRadius);
    source.set(_flagMinAxisRatioKey, model.flagMinAxisRatio);
}


LSST_MEAS_ALGORITHM_PRIVATE_IMPLEMENTATION(FitProfileAlgorithm);

}}}} // namespace lsst::meas::extensions::multiShapelet
