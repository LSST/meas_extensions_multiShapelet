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

#include "lsst/meas/extensions/multiShapelet/FitMultiGaussian.h"
#include "lsst/meas/extensions/multiShapelet/MultiGaussianObjective.h"
#include "lsst/meas/extensions/multiShapelet/MultiGaussianRegistry.h"
#include "lsst/shapelet/ModelBuilder.h"
#include "lsst/afw/math/LeastSquares.h"
#include "lsst/afw/detection/FootprintArray.h"
#include "lsst/afw/detection/FootprintArray.cc"

namespace lsst { namespace meas { namespace extensions { namespace multiShapelet {

FitMultiGaussianModel::FitMultiGaussianModel(
    FitMultiGaussianControl const & ctrl,
    double amplitude_,
    ndarray::Array<double const,1,1> const & parameters
) :
    profile(ctrl.profile), amplitude(amplitude_),
    ellipse(MultiGaussianObjective::EllipseCore(parameters[0], parameters[1], parameters[2])),
    failed(false)
{}

FitMultiGaussianModel::FitMultiGaussianModel(
    FitControl const & ctrl, afw::table::SourceRecord const & source
) :
    profile(ctrl.profile), amplitude(1.0), ellipse(), failed(false)
{
    afw::table::SubSchema s = source.getSchema()[ctrl.name];
    afw::table::SubSchema sPsf = source.getSchema()[ctrl.psfName];
    ellipse = source.get(s.find< afw::table::Moments<float> >("ellipse").key);
    failed = source.get(s.find<afw::table::Flag>("flux.flags").key);
    amplitude = source.get(s.find<afw::table::Flux::MeasTag>("flux"));
    amplitude /= (ellipse.getArea() / afw::geom::PI);
    amplitude /= MultiGaussianComponent::integrate(MultiGaussianRegistry::lookup(profile));
    // TODO: PSF/aperture flux correction
}

FitMultiGaussianModel::FitMultiGaussianModel(FitMultiGaussianModel const & other) :
    profile(other.profile), amplitude(other.amplitude), ellipse(other.ellipse), failed(other.failed)
{}

FitMultiGaussianModel & FitMultiGaussianModel::operator=(FitMultiGaussianModel const & other) {
    if (&other != this) {
        profile = other.profile;
        amplitude = other.amplitude;
        ellipse = other.ellipse;
        failed = other.failed;
    }
    return *this;
}


shapelet::MultiShapeletFunction FitMultiGaussianModel::asMultiShapelet(
    afw::geom::Point2D const & center
) const {
    static double const NORM2 = shapelet::NORMALIZATION * shapelet::NORMALIZATION;
    shapelet::MultiShapeletFunction::ElementList elements;
    MultiGaussianList const & components = MultiGaussianRegistry::lookup(profile);
    for (MultiGaussianList::const_iterator i = components.begin(); i != components.end(); ++i) {
        afw::geom::ellipses::Ellipse fullEllipse(ellipse, center);
        fullEllipse.scale(i->radius);
        elements.push_back(
            shapelet::ShapeletFunction(
                0,
                shapelet::HERMITE,
                fullEllipse,
                i->amplitude * NORM2
            )
        );
    }
    return shapelet::MultiShapeletFunction(elements);
}

FitMultiGaussianAlgorithm::FitMultiGaussianAlgorithm(
    FitMultiGaussianControl const & ctrl,
    afw::table::Schema & schema
) :
    algorithms::Algorithm(ctrl),
    _fluxKeys(afw::table::addFluxFields(schema, ctrl + ".flux", "multi-Gaussian model flux")),
    _ellipseKey(schema.addField< afw::table::Moments<float> >(
                    ctrl.name + ".ellipse",
                    "half-light radius ellipse"
                ))
{}

template <typename PixelT>
ndarray::Array<double,1,1> FitMultiGaussianAlgorithm::processInputs(
    afw::detection::Footprint & footprint,
    afw::image::MaskedImage<PixelT> const & image,
    bool usePixelWeights, afw::image::MaskPixel badPixelMask, int growFootprint
) {
    footprint = *afw::detection::growFootprint(footprint, growFootprint);
    footprint.intersectMask(*image.getMask(), badPixelMask);
    ndarray::Array<double,2,2> result = ndarray::allocate(2,footprint.getArea());
    flattenArray(footprint, image.getImage()->getArray(), result[0].shallow(), image.getXY0());
    flattenArray(footprint, image.getVariance()->getArray(), result[1].shallow(), image.getXY0());
    ndarray::EigenView<double,1,1,Eigen::ArrayXpr> weights(result[1]);
    weights = weights.sqrt().inverse();
    if (!usePixelWeights) {
        weights = weights.mean(); // TODO: mean(weights) or 1/RMS (or ...) here?
    }
    return result;
}

PTR(MultiGaussianObjective) FitMultiGaussianAlgorithm::makeObjective(
    FitMultiGaussianControl const & ctrl,
    FitPsfModel const & psfModel,
    afw::geom::ellipses::Quadrupole const & shape,
    afw::detection::Footprint const & footprint,
    afw::image::MaskedImage<double> const & image,
    afw::geom::Point2D const & center
) {
    PTR(afw::detection::Footprint) grown = afw::detection::growFootprint(footprint, ctrl.growFootprint);
    ndarray::Array<double,1,1> data = afw::detection::flattenArray(
        footprint, image.getImage()->getArray(), image.getXY0()
    );
    return boost::make_shared<MultiGaussianObjective>(
        MultiGaussianRegistry::lookup(ctrl.profile), center,
        footprint,
        ndarray::flatten<1>(ndarray::copy(image.getArray()))
    );
}

HybridOptimizer FitPsfAlgorithm::makeOptimizer(
    FitPsfControl const & ctrl,
    afw::image::Image<double> const & image,
    afw::geom::Point2D const & center
) {
    PTR(Objective) obj = makeObjective(ctrl, image, center);
    HybridOptimizerControl optCtrl; // TODO: nest this in FitPsfControl
    optCtrl.tau = 1E-6;
    optCtrl.useCholesky = true;
    optCtrl.gTol = 1E-6;
    ndarray::Array<double,1,1> initial = ndarray::allocate(obj->getParameterSize());
    MultiGaussianObjective::EllipseCore ellipse(0.0, 0.0, ctrl.initialRadius);
    ellipse.writeParameters(initial.getData());
    return HybridOptimizer(obj, initial, optCtrl);
}

template <typename PixelT>
void FitPsfAlgorithm::_apply(
    afw::table::SourceRecord & source,
    afw::image::Exposure<PixelT> const & exposure,
    afw::geom::Point2D const & center
) const {
    if (!exposure.hasPsf()) {
        throw LSST_EXCEPT(
            pex::exceptions::LogicErrorException,
            "Cannot run FitPsfAlgorithm without a PSF."
        );
    }
    source.set(_flagKey, true);
    FitPsfModel model = apply(getControl(), *exposure.getPsf(), center);
    source[_innerKey] = model.inner;
    source[_outerKey] = model.outer;
    source.set(_ellipseKey, model.ellipse);
    source.set(_flagKey, model.failed);
}

void FitPsfAlgorithm::fitShapeletTerms(
    FitPsfControl const & ctrl,
    afw::image::Image<double> const & image,
    afw::geom::Point2D const & center,
    FitPsfModel & model
) {
    afw::geom::Box2I bbox = image.getBBox(afw::image::PARENT);
    int innerCoeffs = shapelet::computeSize(ctrl.innerOrder);
    int outerCoeffs = shapelet::computeSize(ctrl.outerOrder);
    ndarray::Array<double,2,-2> matrix = ndarray::allocate(bbox.getArea(), innerCoeffs + outerCoeffs);
    afw::geom::ellipses::Ellipse innerEllipse(model.ellipse, center);
    afw::geom::ellipses::Ellipse outerEllipse(model.ellipse, center);
    outerEllipse.getCore().scale(ctrl.radiusRatio);
    shapelet::ModelBuilder innerShapelets(ctrl.innerOrder, innerEllipse, bbox);
    shapelet::ModelBuilder outerShapelets(ctrl.outerOrder, outerEllipse, bbox);
    matrix[ndarray::view()(0, innerCoeffs)] = innerShapelets.getModel();
    matrix[ndarray::view()(innerCoeffs, innerCoeffs + outerCoeffs)] = outerShapelets.getModel();
    ndarray::Array<double,1,1> data = ndarray::flatten<1>(ndarray::copy(image.getArray()));
    afw::math::LeastSquares lstsq = afw::math::LeastSquares::fromDesignMatrix(matrix, data);
    model.inner.deep() = lstsq.getSolution()[ndarray::view(0, innerCoeffs)];
    model.outer.deep() = lstsq.getSolution()[ndarray::view(innerCoeffs, innerCoeffs + outerCoeffs)];
}

FitPsfModel FitPsfAlgorithm::apply(
    FitPsfControl const & ctrl,
    afw::image::Image<double> const & image,
    afw::geom::Point2D const & center
) {
    // First, we fit an elliptical double-Gaussian with fixed radius and amplitude ratios;
    // the free parameters are the amplitude and ellipse core of the inner component.
    // We intentionally separate these steps into public functions so we can reproduce
    // the same results with a pure-Python implementation that lets us visualize what's
    // going on.
    HybridOptimizer opt = makeOptimizer(ctrl, image, center);
    opt.run();
    Model model(
        ctrl, 
        boost::static_pointer_cast<MultiGaussianObjective const>(opt.getObjective())->getAmplitude(),
        opt.getParameters()
    );
    model.failed = !(opt.getState() & HybridOptimizer::SUCCESS);
    fitShapeletTerms(ctrl, image, center, model);
    return model;
}

PTR(algorithms::AlgorithmControl) FitPsfControl::_clone() const {
    return boost::make_shared<FitPsfControl>(*this);
}

PTR(algorithms::Algorithm) FitPsfControl::_makeAlgorithm(
    afw::table::Schema & schema,
    PTR(daf::base::PropertyList) const & metadata
) const {
    return boost::make_shared<FitPsfAlgorithm>(*this, boost::ref(schema));
}

LSST_MEAS_ALGORITHM_PRIVATE_IMPLEMENTATION(FitPsfAlgorithm);


}}}} // namespace lsst::meas::extensions::multiShapelet
