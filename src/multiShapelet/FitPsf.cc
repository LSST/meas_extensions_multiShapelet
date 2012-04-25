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
#include "lsst/meas/extensions/multiShapelet/GaussianObjective.h"
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

template <typename PixelT>
void fillMultiShapeletImage(
    ndarray::Array<PixelT,2,1> const & array, 
    shapelet::MultiShapeletFunction const & msf
) {
    // TODO: could use shapelet::ModelBuilder here; probably faster
    // (low priority because this function is mostly intended for testing purposes).
    shapelet::MultiShapeletFunctionEvaluator ev(msf);
    afw::geom::Point2D point;
    for (
        typename ndarray::Array<PixelT,2,1>::Iterator rowIter = array.begin();
        rowIter != array.end();
        ++rowIter, point.setY(point.getY() + 1.0) // FIXME: overload Point::getY to return reference
    ) {
        for (
            typename ndarray::Array<PixelT,2,1>::Reference::Iterator pixIter = rowIter->begin();
            pixIter != rowIter->end();
            ++pixIter, point.setX(point.getX() + 1.0) // FIXME: overload Point::getX to return reference
        ) {
            *pixIter = ev(point);
        }
        point.setX(0.0);
    }
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
    ellipse = GaussianObjective::EllipseCore(parameters[0], parameters[1], parameters[2]);
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
    inner.deep() = source.get(s.find< afw::table::Array<float> >("inner").key);
    outer.deep() = source.get(s.find< afw::table::Array<float> >("outer").key);
    ellipse = source.get(s.find< afw::table::Moments<float> >("ellipse").key);
    failed = source.get(s.find<afw::table::Flag>("flags").key);
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

template <typename PixelT>
void FitPsfModel::evaluate(
    ndarray::Array<PixelT,2,1> const & array, afw::geom::Point2D const & center
) const {
    shapelet::MultiShapeletFunction msf = asMultiShapelet(center);
    fillMultiShapeletImage(array, msf);
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

PTR(GaussianObjective) FitPsfAlgorithm::makeObjective(
    FitPsfControl const & ctrl,
    afw::image::Image<double> const & image,
    afw::geom::Point2D const & center
) {
    GaussianObjective::ComponentList components;
    components.push_back(GaussianObjective::Component(1.0, 1.0));
    components.push_back(GaussianObjective::Component(ctrl.amplitudeRatio, ctrl.radiusRatio));
    return boost::make_shared<GaussianObjective>(
        components, center, image.getBBox(afw::image::PARENT),
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
    optCtrl.useCholesky = false;
    ndarray::Array<double,1,1> initial = ndarray::allocate(obj->getParameterSize());
    GaussianObjective::EllipseCore ellipse(0.0, 0.0, ctrl.initialRadius);
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
        boost::static_pointer_cast<GaussianObjective const>(opt.getObjective())->getAmplitude(),
        opt.getParameters()
    );
    model.failed = !(opt.getState() & HybridOptimizer::SUCCESS);
    // Now, we free up all the amplitudes (including higher-order shapelet terms), and do a linear-only fit.
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

template
void FitPsfModel::evaluate(ndarray::Array<float,2,1> const & array, afw::geom::Point2D const & center) const;

template
void FitPsfModel::evaluate(ndarray::Array<double,2,1> const & array, afw::geom::Point2D const & center) const;


}}}} // namespace lsst::meas::extensions::multiShapelet
