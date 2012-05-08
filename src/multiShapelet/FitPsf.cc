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

#include "lsst/meas/extensions/multiShapelet/FitPsf.h"
#include "lsst/meas/extensions/multiShapelet/MultiGaussianObjective.h"
#include "lsst/shapelet/ModelBuilder.h"
#include "lsst/afw/math/LeastSquares.h"

namespace lsst { namespace meas { namespace extensions { namespace multiShapelet {

namespace {

int computeOrder(int size2d) {
    // Could just use f.p. square roots here, but we only care about small numbers
    // and I don't want to think about how to be robust against round-off error issues.
    for (int o=0, s=1; s <= size2d; s += (++o + 1)) {
        if (s == size2d) return o;
    }
    throw LSST_EXCEPT(
        pex::exceptions::LogicErrorException,
        (boost::format("%d is not a valid shapelet vector size") % size2d).str()
    );
}

} // anonymous

FitPsfModel::FitPsfModel(
    FitPsfControl const & ctrl,
    double amplitude,
    ndarray::Array<double const,1,1> const & parameters
) :
    inner(ndarray::allocate(shapelet::computeSize(ctrl.innerOrder))),
    outer(ndarray::allocate(shapelet::computeSize(ctrl.outerOrder))),
    ellipse(),
    radiusRatio(ctrl.radiusRatio),
    failed(false)
{
    static double const NORM2 = shapelet::NORMALIZATION * shapelet::NORMALIZATION;
    inner.deep() = 0.0;
    outer.deep() = 0.0;
    ellipse = MultiGaussianObjective::EllipseCore(parameters[0], parameters[1], parameters[2]);
    inner[0] = amplitude / NORM2;
    outer[0] = amplitude * ctrl.amplitudeRatio / NORM2;
}

FitPsfModel::FitPsfModel(FitPsfControl const & ctrl, afw::table::SourceRecord const & source) :
    inner(ndarray::allocate(shapelet::computeSize(ctrl.innerOrder))),
    outer(ndarray::allocate(shapelet::computeSize(ctrl.outerOrder))),
    ellipse(),
    radiusRatio(ctrl.radiusRatio),
    failed(false)
{
    afw::table::SubSchema s = source.getSchema()[ctrl.name];
    afw::table::Key< afw::table::Array<float> > innerKey = s.find< afw::table::Array<float> >("inner").key;
    if (innerKey.getSize() <= inner.getSize<0>()) {
        inner.deep() = 0.0;
        inner[ndarray::view(0, innerKey.getSize())] = source.get(innerKey);
    } else {
        inner.deep() = source.get(innerKey)[ndarray::view(0, inner.getSize<0>())];
    }
    afw::table::Key< afw::table::Array<float> > outerKey = s.find< afw::table::Array<float> >("outer").key; 
    if (outerKey.getSize() <= outer.getSize<0>()) {
        outer.deep() = 0.0;
        outer[ndarray::view(0, outerKey.getSize())] = source.get(outerKey);
    } else {
        outer.deep() = source.get(outerKey)[ndarray::view(0, outer.getSize<0>())];
    }
    ellipse = source.get(s.find< afw::table::Moments<float> >("ellipse").key);
    failed = source.get(s.find<afw::table::Flag>("flags").key);
}

FitPsfModel::FitPsfModel(FitPsfModel const & other) :
    inner(ndarray::copy(other.inner)),
    outer(ndarray::copy(other.outer)),
    ellipse(other.ellipse),
    radiusRatio(other.radiusRatio),
    failed(other.failed)
{}


FitPsfModel & FitPsfModel::operator=(FitPsfModel const & other) {
    if (&other != this) {
        inner = ndarray::copy(other.inner);
        outer = ndarray::copy(other.outer);
        ellipse = other.ellipse;
        radiusRatio = other.radiusRatio;
        failed = other.failed;
    }
    return *this;
}


shapelet::MultiShapeletFunction FitPsfModel::asMultiShapelet(
    afw::geom::Point2D const & center
) const {
    shapelet::MultiShapeletFunction::ElementList elements;
    afw::geom::ellipses::Ellipse fullEllipse(ellipse, center);
    // FIXME: should have data-type cast functionality in ndarray
    ndarray::Array<shapelet::Pixel,1,1> sInner(ndarray::allocate(inner.getSize<0>()));
    ndarray::Array<shapelet::Pixel,1,1> sOuter(ndarray::allocate(outer.getSize<0>()));
    sInner.deep() = inner;
    sOuter.deep() = outer;
    elements.push_back(
        shapelet::ShapeletFunction(
            computeOrder(inner.getSize<0>()),
            shapelet::HERMITE,
            fullEllipse,
            sInner
        )
    );
    fullEllipse.scale(radiusRatio);
    elements.push_back(
        shapelet::ShapeletFunction(
            computeOrder(outer.getSize<0>()),
            shapelet::HERMITE,
            fullEllipse,
            sOuter
        )
    );
    return shapelet::MultiShapeletFunction(elements);
}

FitPsfAlgorithm::FitPsfAlgorithm(FitPsfControl const & ctrl, afw::table::Schema & schema) :
    algorithms::Algorithm(ctrl),
    _innerKey(schema.addField< afw::table::Array<float> >(
                  ctrl.name + ".inner",
                  "Gauss-Hermite coefficients of the inner expansion (see shapelet)",
                  shapelet::computeSize(ctrl.innerOrder)
              )),
    _outerKey(schema.addField< afw::table::Array<float> >(
                  ctrl.name + ".outer",
                  "Gauss-Hermite coefficients of the outer expansion (see shapelet)",
                  shapelet::computeSize(ctrl.outerOrder)
              )),
    _ellipseKey(schema.addField< afw::table::Moments<float> >(
                    ctrl.name + ".ellipse",
                    "Ellipse corresponding to the inner expansion"
                )),
    _flagKey(schema.addField<afw::table::Flag>(
                 ctrl.name + ".flags",
                 "set if the multi-shapelet PSF fit was unsuccessful"
             ))
{}

PTR(MultiGaussianObjective) FitPsfAlgorithm::makeObjective(
    FitPsfControl const & ctrl,
    ModelInputHandler const & inputs
) {
    MultiGaussianList components;
    components.push_back(MultiGaussianComponent(1.0, 1.0));
    components.push_back(MultiGaussianComponent(ctrl.amplitudeRatio, ctrl.radiusRatio));
    return boost::make_shared<MultiGaussianObjective>(inputs, components);
}

HybridOptimizer FitPsfAlgorithm::makeOptimizer(
    FitPsfControl const & ctrl,
    ModelInputHandler const & inputs
) {
    PTR(Objective) obj = makeObjective(ctrl, inputs);
    HybridOptimizerControl optCtrl; // TODO: nest this in FitPsfControl
    optCtrl.tau = 1E-6;
    optCtrl.useCholesky = true;
    optCtrl.gTol = 1E-6;
    ndarray::Array<double,1,1> initial = ndarray::allocate(obj->getParameterSize());
    MultiGaussianObjective::EllipseCore ellipse(0.0, 0.0, ctrl.initialRadius);
    ellipse.writeParameters(initial.getData());
    return HybridOptimizer(obj, initial, optCtrl);
}

void FitPsfAlgorithm::fitShapeletTerms(
    FitPsfControl const & ctrl,
    ModelInputHandler const & inputs,
    FitPsfModel & model
) {
    int innerCoeffs = shapelet::computeSize(ctrl.innerOrder);
    int outerCoeffs = shapelet::computeSize(ctrl.outerOrder);
    ndarray::Array<double,2,-2> matrix = ndarray::allocate(inputs.getSize(), innerCoeffs + outerCoeffs);
    shapelet::ModelBuilder innerShapelets(ctrl.innerOrder, inputs.getX(), inputs.getY());
    shapelet::ModelBuilder outerShapelets(ctrl.outerOrder, inputs.getX(), inputs.getY());
    innerShapelets.update(model.ellipse);
    model.ellipse.scale(ctrl.radiusRatio);
    outerShapelets.update(model.ellipse);
    model.ellipse.scale(1.0 / ctrl.radiusRatio);
    matrix[ndarray::view()(0, innerCoeffs)] = innerShapelets.getModel();
    matrix[ndarray::view()(innerCoeffs, innerCoeffs + outerCoeffs)] = outerShapelets.getModel();
    afw::math::LeastSquares lstsq = afw::math::LeastSquares::fromDesignMatrix(matrix, inputs.getData());
    model.inner.deep() = lstsq.getSolution()[ndarray::view(0, innerCoeffs)];
    model.outer.deep() = lstsq.getSolution()[ndarray::view(innerCoeffs, innerCoeffs + outerCoeffs)];
}

FitPsfModel FitPsfAlgorithm::apply(
    FitPsfControl const & ctrl,
    ModelInputHandler const & inputs
) {
    // First, we fit an elliptical double-Gaussian with fixed radius and amplitude ratios;
    // the free parameters are the amplitude and ellipse core of the inner component.
    // We intentionally separate these steps into public functions so we can reproduce
    // the same results with a pure-Python implementation that lets us visualize what's
    // going on.
    HybridOptimizer opt = makeOptimizer(ctrl, inputs);
    opt.run();
    Model model(
        ctrl, 
        boost::static_pointer_cast<MultiGaussianObjective const>(opt.getObjective())->getAmplitude(),
        opt.getParameters()
    );
    model.failed = !(opt.getState() & HybridOptimizer::SUCCESS);
    fitShapeletTerms(ctrl, inputs, model);
    return model;
}

FitPsfModel FitPsfAlgorithm::apply(
    FitPsfControl const & ctrl,
    afw::detection::Psf const & psf,
    afw::geom::Point2D const & center
) {
    PTR(afw::image::Image<afw::math::Kernel::Pixel>) image = psf.computeImage(center);
    ModelInputHandler inputs(*image, center, image->getBBox(afw::image::PARENT));
    return apply(ctrl, inputs);
}

template <typename PixelT>
void FitPsfAlgorithm::_apply(
    afw::table::SourceRecord & source,
    afw::image::Exposure<PixelT> const & exposure,
    afw::geom::Point2D const & center
) const {
    source.set(_flagKey, true);
    if (!exposure.hasPsf()) {
        throw LSST_EXCEPT(
            pex::exceptions::LogicErrorException,
            "Cannot run FitPsfAlgorithm without a PSF."
        );
    }
    FitPsfModel model = apply(getControl(), *exposure.getPsf(), center);
    source[_innerKey] = model.inner;
    source[_outerKey] = model.outer;
    source.set(_ellipseKey, model.ellipse);
    source.set(_flagKey, model.failed);
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
