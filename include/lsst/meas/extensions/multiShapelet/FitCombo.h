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
#ifndef MULTISHAPELET_FitCombo_h_INCLUDED
#define MULTISHAPELET_FitCombo_h_INCLUDED

#include "lsst/meas/extensions/multiShapelet/FitProfile.h"

namespace lsst { namespace meas { namespace extensions { namespace multiShapelet {

class FitComboAlgorithm;

class FitComboControl : public algorithms::AlgorithmControl {
public:

    LSST_CONTROL_FIELD(componentNames, std::vector<std::string>,
                       "Root names of the FitProfileAlgorithm fields to be fit together.");
    LSST_CONTROL_FIELD(psfName, std::string, "Root name of the FitPsfAlgorithm fields.");
    LSST_CONTROL_FIELD(usePixelWeights, bool,
                       "If true, individually weigh pixels using the variance image.");
    LSST_CONTROL_FIELD(badMaskPlanes, std::vector<std::string>,
                       "Mask planes that indicate pixels that should be ignored in the fit.");
    LSST_CONTROL_FIELD(growFootprint, int, "Number of pixels to grow the footprint by.");
    LSST_CONTROL_FIELD(radiusInputFactor, double,
                       "Number of half-light radii used to determine the pixels to fit");

    PTR(FitComboControl) clone() const {
        return boost::static_pointer_cast<FitComboControl>(_clone());
    }

    PTR(FitComboAlgorithm) makeAlgorithm(
        afw::table::Schema & schema,
        PTR(daf::base::PropertyList) const & metadata = PTR(daf::base::PropertyList)(),
        algorithms::AlgorithmMap const & others = algorithms::AlgorithmMap()
    ) const;

    FitComboControl() :
        algorithms::AlgorithmControl("multishapelet.combo", 2.6),
        componentNames(), psfName("multishapelet.psf"),
        usePixelWeights(false), badMaskPlanes(), growFootprint(5), radiusInputFactor(4.0)
    {
        componentNames.push_back("multishapelet.exp");
        componentNames.push_back("multishapelet.dev");
        badMaskPlanes.push_back("BAD");
        badMaskPlanes.push_back("SAT");
        badMaskPlanes.push_back("INTRP");
        badMaskPlanes.push_back("CR");
    }

private:
    virtual PTR(algorithms::AlgorithmControl) _clone() const;
    virtual PTR(algorithms::Algorithm) _makeAlgorithm(
        afw::table::Schema & schema,
        PTR(daf::base::PropertyList) const & metadata,
        algorithms::AlgorithmMap const & other
    ) const;
};

/**
 *  @brief A linear combination of several FitProfileModels
 */
struct FitComboModel {

    ndarray::Array<float,1,1> components; ///< relative flux of each component (always sums to one)
    double flux; ///< total flux of model, integrated to infinity (includes PSF factor, if enabled)
    double fluxErr; ///< uncertainty on flux
    double chisq; ///< reduced chi^2

    explicit FitComboModel(FitComboControl const & ctrl);

    /// @brief Deep copy constructor.
    FitComboModel(FitComboModel const & other);

    /// @brief Deep assignment operator.
    FitComboModel & operator=(FitComboModel const & other);

};

class FitComboAlgorithm :
        public algorithms::Algorithm 
#ifndef SWIG
        , public algorithms::ScaledFlux
#endif
{
public:

    typedef FitComboControl Control;
    typedef FitComboModel Model;

    /// @brief Construct an algorithm instance and register its fields with a Schema.
    FitComboAlgorithm(
        FitComboControl const & ctrl,
        afw::table::Schema & schema,
        algorithms::AlgorithmMap const & others
    );

    /// @brief Return the control object
    FitComboControl const & getControl() const {
        return static_cast<FitComboControl const &>(algorithms::Algorithm::getControl());
    }

#ifndef SWIG
    virtual afw::table::KeyTuple<afw::table::Flux> getFluxKeys(int n=0) const { return _fluxKeys; }
    virtual ScaledFlux::KeyTuple getFluxCorrectionKeys(int n=0) const { return _fluxCorrectionKeys; }
#endif

    template <typename PixelT>
    static ModelInputHandler adjustInputs(
        FitComboControl const & ctrl,
        FitPsfModel const & psfModel,
        std::vector<FitProfileModel> const & components,
        afw::detection::Footprint const & footprint,
        afw::image::MaskedImage<PixelT> const & image,
        afw::geom::Point2D const & center
    );

    static FitComboModel apply(
        FitComboControl const & ctrl,
        FitPsfModel const & psfModel,
        std::vector<FitProfileModel> const & components,
        ModelInputHandler const & inputs
    );

private:

    template <typename PixelT>
    void _apply(
        afw::table::SourceRecord & source,
        afw::image::Exposure<PixelT> const & exposure,
        afw::geom::Point2D const & center
    ) const;

    LSST_MEAS_ALGORITHM_PRIVATE_INTERFACE(FitComboAlgorithm);

    afw::table::KeyTuple< afw::table::Flux > _fluxKeys;
    algorithms::ScaledFlux::KeyTuple _fluxCorrectionKeys;
    afw::table::Key< afw::table::Array<float> > _componentsKey;
    afw::table::Key< float > _chisqKey;
    std::vector<CONST_PTR(FitProfileControl)> _componentCtrl;
    CONST_PTR(FitPsfControl) _psfCtrl;
};

inline PTR(FitComboAlgorithm) FitComboControl::makeAlgorithm(
    afw::table::Schema & schema,
    PTR(daf::base::PropertyList) const & metadata,
    algorithms::AlgorithmMap const & others
) const {
    return boost::static_pointer_cast<FitComboAlgorithm>(_makeAlgorithm(schema, metadata, others));
}

}}}} // namespace lsst::meas::extensions::multiShapelet

#endif // !MULTISHAPELET_FitCombo_h_INCLUDED
