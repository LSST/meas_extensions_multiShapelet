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
    FitMultiGaussianControl const & ctrl, afw::table::SourceRecord const & source
) :
    profile(ctrl.profile), amplitude(1.0), ellipse(), failed(false)
{
    afw::table::SubSchema s = source.getSchema()[ctrl.name];
    afw::table::SubSchema sPsf = source.getSchema()[ctrl.psfName];
    ellipse = source.get(s.find< afw::table::Moments<float> >("ellipse").key);
    failed = source.get(s.find<afw::table::Flag>("flux.flags").key);
    amplitude = source.get(s.find<afw::table::Flux::MeasTag>("flux").key);
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
                fullEllipse
            )
        );
        elements.back().getCoefficients()[0] = i->amplitude * NORM2;
    }
    return shapelet::MultiShapeletFunction(elements);
}

FitMultiGaussianAlgorithm::FitMultiGaussianAlgorithm(
    FitMultiGaussianControl const & ctrl,
    afw::table::Schema & schema
) :
    algorithms::Algorithm(ctrl),
    _fluxKeys(afw::table::addFluxFields(schema, ctrl.name + ".flux", "multi-Gaussian model flux")),
    _ellipseKey(schema.addField< afw::table::Moments<float> >(
                    ctrl.name + ".ellipse",
                    "half-light radius ellipse"
                ))
{}

template <typename PixelT>
void FitMultiGaussianAlgorithm::_apply(
    afw::table::SourceRecord & source,
    afw::image::Exposure<PixelT> const & exposure,
    afw::geom::Point2D const & center
) const {
    source.set(_fluxKeys.err, true);
    if (!exposure.hasPsf()) {
        throw LSST_EXCEPT(
            pex::exceptions::LogicErrorException,
            "Cannot run FitMultiGaussianAlgorithm without a PSF."
        );
    }
    // TODO
}

PTR(algorithms::AlgorithmControl) FitMultiGaussianControl::_clone() const {
    return boost::make_shared<FitMultiGaussianControl>(*this);
}

PTR(algorithms::Algorithm) FitMultiGaussianControl::_makeAlgorithm(
    afw::table::Schema & schema,
    PTR(daf::base::PropertyList) const & metadata
) const {
    return boost::make_shared<FitMultiGaussianAlgorithm>(*this, boost::ref(schema));
}

LSST_MEAS_ALGORITHM_PRIVATE_IMPLEMENTATION(FitMultiGaussianAlgorithm);

}}}} // namespace lsst::meas::extensions::multiShapelet
