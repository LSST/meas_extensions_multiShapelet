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
    chisq(std::numeric_limits<double>::quiet_NaN()),
    failedMaxIter(false), failedTinyStep(false), atMinRadius(false), failedMinAxisRatio(false)
{}

FitProfileModel::FitProfileModel(
    FitProfileControl const & ctrl, afw::table::SourceRecord const & source
) :
    profile(ctrl.profile), flux(1.0), fluxErr(0.0), ellipse(),
    chisq(std::numeric_limits<double>::quiet_NaN()),
    failedMaxIter(false), failedTinyStep(false), atMinRadius(false), failedMinAxisRatio(false)
{
    afw::table::SubSchema s = source.getSchema()[ctrl.name];
    flux = source.get(s.find< double >("flux").key);
    fluxErr = source.get(s.find< double >("flux.err").key);
    ellipse = source.get(s.find< afw::table::Moments<float> >("ellipse").key);
    chisq = source.get(s.find<float>("chisq").key);
    failedMaxIter = source.get(s.find<afw::table::Flag>("flags.maxiter").key);
    failedTinyStep = source.get(s.find<afw::table::Flag>("flags.tinystep").key);
    atMinRadius = source.get(s.find<afw::table::Flag>("flags.constraint.r").key);
    failedMinAxisRatio = source.get(s.find<afw::table::Flag>("flags.constraint.q").key);
}

FitProfileModel::FitProfileModel(FitProfileModel const & other) :
    profile(other.profile), flux(other.flux), fluxErr(other.fluxErr), ellipse(other.ellipse),
    chisq(other.chisq),
    failedMaxIter(other.failedMaxIter),
    failedTinyStep(other.failedTinyStep),
    atMinRadius(other.atMinRadius),
    failedMinAxisRatio(other.failedMinAxisRatio)    
{}

FitProfileModel & FitProfileModel::operator=(FitProfileModel const & other) {
    if (&other != this) {
        profile = other.profile;
        flux = other.flux;
        fluxErr = other.fluxErr;
        ellipse = other.ellipse;
        chisq = other.chisq;
        failedMaxIter = other.failedMaxIter;
        failedTinyStep = other.failedTinyStep;
        atMinRadius = other.atMinRadius;
        failedMinAxisRatio = other.failedMinAxisRatio;
    }
    return *this;
}

shapelet::MultiShapeletFunction FitProfileModel::asMultiShapelet(
    afw::geom::Point2D const & center
) const {
    shapelet::MultiShapeletFunction::ElementList elements;
    MultiGaussian const & components = MultiGaussianRegistry::lookup(profile);
    for (MultiGaussian::const_iterator i = components.begin(); i != components.end(); ++i) {
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
    _fluxKey(
        schema.addField<double>(
            ctrl.name + ".flux", "surface brightness at half-light radius", "dn/pix^2"
        )),
    _fluxErrKey(
        schema.addField<double>(
            ctrl.name + ".flux.err", "uncertainty on flux", "dn/pix^2"
        )),
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
    _flagKey(
        schema.addField< afw::table::Flag >(
            ctrl.name + ".flags",
            "error flag; set if model fit failed in any way"
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
        inputs, ctrl.getComponents(), psfModel.getComponents(), psfModel.ellipse,
        ctrl.minRadius, ctrl.minAxisRatio
    );
}

HybridOptimizer FitProfileAlgorithm::makeOptimizer(
    FitProfileControl const & ctrl,
    FitPsfModel const & psfModel,
    afw::geom::ellipses::Quadrupole const & shape,
    ModelInputHandler const & inputs
) {
    MultiGaussianObjective::EllipseCore ellipse(shape);
    if (ctrl.deconvolveShape) {
        ellipse = ctrl.getComponents().deconvolve(
            shape, psfModel.ellipse, psfModel.getComponents()
        );
    }
    // We never want to start with an ellipse smaller than the PSF or an ellipticity
    // on the constraint, because we might never find our way out.
    std::pair<bool,bool> constrained = MultiGaussianObjective::constrainEllipse(
        ellipse, psfModel.ellipse.getTraceRadius(), ctrl.minAxisRatio
    );
    if (constrained.first || constrained.second) {
        ellipse = psfModel.ellipse;
    }
    PTR(Objective) obj = makeObjective(ctrl, psfModel, inputs);
    ndarray::Array<double,1,1> initial = ndarray::allocate(obj->getParameterSize());
    ellipse.writeParameters(initial.getData());
    HybridOptimizerControl optCtrl; // TODO: nest this in FitProfileControl
    optCtrl.tau = 1E-2;
    optCtrl.useCholesky = true;
    optCtrl.gTol = 1E-4;
    return HybridOptimizer(obj, initial, optCtrl);
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
    afw::geom::ellipses::Quadrupole const & shape,
    ModelInputHandler const & inputs
) {
    HybridOptimizer opt = makeOptimizer(ctrl, psfModel, shape, inputs);
    opt.run();
    Model model(
        ctrl, 
        boost::static_pointer_cast<MultiGaussianObjective const>(opt.getObjective())->getAmplitude(),
        opt.getParameters()
    );
    MultiGaussianObjective::EllipseCore ellipse = MultiGaussianObjective::readParameters(opt.getParameters());
    std::pair<bool,bool> constrained 
        = MultiGaussianObjective::constrainEllipse(ellipse, ctrl.minRadius, ctrl.minAxisRatio);
    model.failedMaxIter = opt.getState() & HybridOptimizer::FAILURE_MAXITER;
    model.failedTinyStep = (opt.getState() & HybridOptimizer::FAILURE_MINSTEP)
        || (opt.getState() & HybridOptimizer::FAILURE_MINTRUST);
    model.atMinRadius = constrained.first;
    model.failedMinAxisRatio = constrained.second;
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
    afw::geom::ellipses::Quadrupole const & shape,
    afw::detection::Footprint const & footprint,
    afw::image::MaskedImage<PixelT> const & image,
    afw::geom::Point2D const & center
) {
    afw::image::MaskPixel badPixelMask = afw::image::Mask<>::getPlaneBitMask(ctrl.badMaskPlanes);
    ModelInputHandler inputs(image, center, footprint, ctrl.growFootprint, 
                             badPixelMask, ctrl.usePixelWeights);
    return apply(ctrl, psfModel, shape, inputs);
}

template <typename PixelT>
void FitProfileAlgorithm::_apply(
    afw::table::SourceRecord & source,
    afw::image::Exposure<PixelT> const & exposure,
    afw::geom::Point2D const & center
) const {
    source.set(_flagKey, true);
    if (!exposure.hasPsf()) {
        throw LSST_EXCEPT(
            pex::exceptions::LogicErrorException,
            "Cannot run FitProfileAlgorithm without a PSF."
        );
    }
    FitPsfModel psfModel(*_psfCtrl, source);
    FitProfileModel model = apply(
        getControl(), psfModel,
        source.getShape(), *source.getFootprint(),
        exposure.getMaskedImage(), center
    );
    source.set(_fluxKey, model.flux);
    source.set(_fluxErrKey, model.fluxErr);
    source.set(_ellipseKey, model.ellipse);
    source.set(_chisqKey, model.chisq);
    source.set(_flagMaxIterKey, model.failedMaxIter);
    source.set(_flagTinyStepKey, model.failedTinyStep);
    source.set(_flagMinRadiusKey, model.atMinRadius);
    source.set(_flagMinAxisRatioKey, model.failedMinAxisRatio);
    // Note: don't consider tinystep to be a failure mode right now; in practice
    // it doesn't seem to indicate a problem.
    bool failed = model.failedMinAxisRatio || model.failedMaxIter;
    source.set(_flagKey, failed);
}


LSST_MEAS_ALGORITHM_PRIVATE_IMPLEMENTATION(FitProfileAlgorithm);

}}}} // namespace lsst::meas::extensions::multiShapelet
