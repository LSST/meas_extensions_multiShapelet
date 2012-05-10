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

FitProfileModel::FitProfileModel(
    FitProfileControl const & ctrl,
    double amplitude_,
    ndarray::Array<double const,1,1> const & parameters
) :
    profile(ctrl.profile), amplitude(amplitude_), amplitudeErr(0.0),
    ellipse(MultiGaussianObjective::EllipseCore(parameters[0], parameters[1], parameters[2])),
    failed(false)
{}

FitProfileModel::FitProfileModel(
    FitProfileControl const & ctrl, afw::table::SourceRecord const & source
) :
    profile(ctrl.profile), amplitude(1.0), amplitudeErr(0.0), ellipse(), failed(false)
{
    afw::table::SubSchema s = source.getSchema()[ctrl.name];
    amplitude = source.get(s.find< double >("amplitude").key);
    amplitudeErr = source.get(s.find< double >("amplitude.err").key);
    ellipse = source.get(s.find< afw::table::Moments<float> >("ellipse").key);
    failed = source.get(s.find< afw::table::Flag >("flags").key);
}

FitProfileModel::FitProfileModel(FitProfileModel const & other) :
    profile(other.profile), amplitude(other.amplitude), amplitudeErr(other.amplitudeErr),
    ellipse(other.ellipse), failed(other.failed)
{}

FitProfileModel & FitProfileModel::operator=(FitProfileModel const & other) {
    if (&other != this) {
        profile = other.profile;
        amplitude = other.amplitude;
        amplitudeErr = other.amplitudeErr;
        ellipse = other.ellipse;
        failed = other.failed;
    }
    return *this;
}

shapelet::MultiShapeletFunction FitProfileModel::asMultiShapelet(
    afw::geom::Point2D const & center
) const {
    shapelet::MultiShapeletFunction::ElementList elements;
    MultiGaussianList const & components = MultiGaussianRegistry::lookup(profile);
    for (MultiGaussianList::const_iterator i = components.begin(); i != components.end(); ++i) {
        afw::geom::ellipses::Ellipse fullEllipse(ellipse, center);
        elements.push_back(i->makeShapelet(fullEllipse));
        elements.back().getCoefficients().asEigen() *= amplitude;
    }
    return shapelet::MultiShapeletFunction(elements);
}

FitProfileAlgorithm::FitProfileAlgorithm(
    FitProfileControl const & ctrl,
    afw::table::Schema & schema,
    algorithms::AlgorithmControlMap const & others
) :
    algorithms::Algorithm(ctrl),
    _amplitudeKey(schema.addField<double>(
                      ctrl.name + ".amplitude", "surface brightness at half-light radius", "dn/pix^2"
                  )),
    _amplitudeErrKey(schema.addField<double>(
                         ctrl.name + ".amplitude.err", "uncertainty on amplitude", "dn/pix^2"
                     )),
    _ellipseKey(schema.addField< afw::table::Moments<float> >(
                    ctrl.name + ".ellipse",
                    "half-light radius ellipse"
                )),
    _flagKey(schema.addField< afw::table::Flag >(
                 ctrl.name + ".flags",
                 "error flags; set if model fit failed in any way"
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
    // TODO
}

PTR(MultiGaussianObjective) FitProfileAlgorithm::makeObjective(
    FitProfileControl const & ctrl,
    FitPsfModel const & psfModel,
    afw::geom::ellipses::Quadrupole const & shape,
    ModelInputHandler const & inputs
) {
    return boost::make_shared<MultiGaussianObjective>(
        inputs, ctrl.getComponents(), psfModel.getComponents(), psfModel.ellipse
    );
}


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

LSST_MEAS_ALGORITHM_PRIVATE_IMPLEMENTATION(FitProfileAlgorithm);

}}}} // namespace lsst::meas::extensions::multiShapelet
