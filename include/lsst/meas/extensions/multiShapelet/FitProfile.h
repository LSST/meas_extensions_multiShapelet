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
#ifndef MULTISHAPELET_FitProfile_h_INCLUDED
#define MULTISHAPELET_FitProfile_h_INCLUDED

#include "lsst/meas/extensions/multiShapelet/FitPsf.h"
#include "lsst/meas/extensions/multiShapelet/MultiGaussianRegistry.h"

namespace lsst { namespace meas { namespace extensions { namespace multiShapelet {

class FitProfileAlgorithm;

class FitProfileControl : public algorithms::AlgorithmControl {
public:

    LSST_CONTROL_FIELD(profile, std::string, "Name of a registered multi-Gaussian profile.");
    LSST_CONTROL_FIELD(psfName, std::string, "Root name of the FitPsfAlgorithm fields.");
    LSST_CONTROL_FIELD(usePsfShapeletTerms, bool,
                       "If true, use the full double-shapelet PSF model; if false, use double-Gaussian.");
    LSST_CONTROL_FIELD(minRadius, double, "Minimum half-light radius in units of PSF inner radius.");
    LSST_CONTROL_FIELD(minAxisRatio, double, "Minimum axis ratio for ellipse (b/a).");
    LSST_CONTROL_FIELD(deconvolveShape, bool,
                       "Attempt to approximately deconvolve the canonical shape before "
                       "using it to set the initial parameters.");
    LSST_CONTROL_FIELD(usePixelWeights, bool,
                       "If true, individually weigh pixels using the variance image.");
    LSST_CONTROL_FIELD(badMaskPlanes, std::vector<std::string>,
                       "Mask planes that indicate pixels that should be ignored in the fit.");
    LSST_CONTROL_FIELD(growFootprint, int,
                       "Number of pixels to grow the footprint by.");

    PTR(FitProfileControl) clone() const {
        return boost::static_pointer_cast<FitProfileControl>(_clone());
    }

    PTR(FitProfileAlgorithm) makeAlgorithm(
        afw::table::Schema & schema,
        PTR(daf::base::PropertyList) const & metadata = PTR(daf::base::PropertyList)(),
        algorithms::AlgorithmControlMap const & others = algorithms::AlgorithmControlMap()
    ) const;

    MultiGaussian const & getComponents() const { return MultiGaussianRegistry::lookup(profile); }

    FitProfileControl() :
        algorithms::AlgorithmControl("multishapelet.profile", 2.5),
        profile("tractor-exponential"), psfName("multishapelet.psf"),
        usePsfShapeletTerms(false), minRadius(0.01), minAxisRatio(0.05), deconvolveShape(true),
        usePixelWeights(false), badMaskPlanes(), growFootprint(3)
    {
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
        algorithms::AlgorithmControlMap const & other
    ) const;
};

/**
 *  @brief An elliptical model composed of several Gaussians.
 *
 *  All of the Gaussian components have the same ellipticity, and their fluxes and radius ratios
 *  are fixed relative to each other.  Different combinations of flux and radius ratios are
 *  registered with the MultiGaussianRegistry class, and are usually approximations to Sersic functions
 *  with the ellipse defined at the half-light radius of the exact Sersic function being approximated.
 */
struct FitProfileModel {

    std::string profile; ///< name of profile to look up in MultiGaussianRegistry
    double flux; ///< total flux of approximate model, integrated to infinity
    double fluxErr; ///< uncertainty on flux
    afw::geom::ellipses::Quadrupole ellipse; ///< half-light radius ellipse
    double chisq; ///< reduced chi^2
    bool failedMaxIter; ///< set to true if the optimizer hit the maximum number of iterations
    bool failedTinyStep; ///< set to true if the optimizer step size got too small to make progress
    bool atMinRadius; ///< set to true if the best-fit radius was at the minimum constraint (not a failure)
    bool failedMinAxisRatio; ///< set to true if the best-fit axis ratio was at the minimum constraint

    FitProfileModel(
        FitProfileControl const & ctrl,
        double amplitude,
        ndarray::Array<double const,1,1> const & parameters
    );

    /// @brief Construct by extracting saved values from a SourceRecord.
    FitProfileModel(FitProfileControl const & ctrl, afw::table::SourceRecord const & source);

    /// @brief Deep copy constructor.
    FitProfileModel(FitProfileModel const & other);

    /// @brief Deep assignment operator.
    FitProfileModel & operator=(FitProfileModel const & other);

    MultiGaussian const & getComponents() const { return MultiGaussianRegistry::lookup(profile); }

    /**
     *  @brief Return a MultiShapeletFunction representation of the model (unconvolved).
     */
    shapelet::MultiShapeletFunction asMultiShapelet(
        afw::geom::Point2D const & center = afw::geom::Point2D()
    ) const;

};

class FitProfileAlgorithm : public algorithms::Algorithm {
public:

    typedef FitProfileControl Control;
    typedef FitProfileModel Model;

    /// @brief Construct an algorithm instance and register its fields with a Schema.
    FitProfileAlgorithm(
        FitProfileControl const & ctrl,
        afw::table::Schema & schema,
        algorithms::AlgorithmControlMap const & others
    );

    /// @brief Return the control object
    FitProfileControl const & getControl() const {
        return static_cast<FitProfileControl const &>(algorithms::Algorithm::getControl());
    }

    /**
     *  @brief Return an Objective that can be used to fit the convolved model to an image.
     *
     *  The objective function will use only the Gaussian terms in the PSF for the convolution.
     *
     *  This is provided primarily for testing and debugging purposes.
     */
    static PTR(MultiGaussianObjective) makeObjective(
        FitProfileControl const & ctrl,
        FitPsfModel const & psfModel,
        ModelInputHandler const & inputs
    );

    /**
     *  @brief Return an optimizer that can be used to fit the convolved model image.
     *
     *  The optimizer's objective function will use only the Gaussian terms in the PSF for the convolution.
     *
     *  This is provided primarily for testing and debugging purposes; the user can create an optimizer,
     *  step through it, and use the FitProfileModel constructor that takes a parameter vector to
     *  visualize its progress.
     */
    static HybridOptimizer makeOptimizer(
        FitProfileControl const & ctrl,
        FitPsfModel const & psfModel,
        afw::geom::ellipses::Quadrupole const & shape,
        ModelInputHandler const & inputs
    );

    /**
     *  @brief Given a model computed using only the double-Gaussian PSF approximation,
     *         do a linear fit with additional shapelet terms in the PSF.
     *
     *  This is provided primarily for testing and debugging purposes; it is the second
     *  part of apply(), after the nonlinear double-Gaussian fit.
     */
    static void fitShapeletTerms(
        FitProfileControl const & ctrl,
        FitPsfModel const & psfModel,
        ModelInputHandler const & inputs,
        FitProfileModel & model
    );

    /**
     *  @brief Fit to a local PSF or kernel image.
     *
     *  @param[in]     ctrl           Details of the model to fit.
     *  @param[in]     psfModel       Localized double-shapelet PSF model.
     *  @param[in,out] shape          Shape measurement used to set initial ellipse parameters
     *                                (possibly modified as defined by ctrl data members).
     *  @param[in]     inputs         Inputs that determine the data to be fit.
     */
    static FitProfileModel apply(
        FitProfileControl const & ctrl,
        FitPsfModel const & psfModel,
        afw::geom::ellipses::Quadrupole const & shape,
        ModelInputHandler const & inputs
    );

    /**
     *  @brief Fit to a local PSF or kernel image.
     *
     *  @param[in]     ctrl           Details of the model to fit.
     *  @param[in]     psfModel       Localized double-shapelet PSF model.
     *  @param[in,out] shape          Shape measurement used to set initial ellipse parameters
     *                                (possibly modified as defined by ctrl data members).
     *  @param[in]     footprint      Region of the image to fit (after modification according to
     *                                ctrl parameters).
     *  @param[in]     image          Full masked image to fit.
     *  @param[in]     center         Center of the object in the image's PARENT coordinate system
     *                                (i.e. xy0 is used).
     */
    template <typename PixelT>
    static FitProfileModel apply(
        FitProfileControl const & ctrl,
        FitPsfModel const & psfModel,
        afw::geom::ellipses::Quadrupole const & shape,
        afw::detection::Footprint const & footprint,
        afw::image::MaskedImage<PixelT> const & image,
        afw::geom::Point2D const & center
    );

private:

    template <typename PixelT>
    void _apply(
        afw::table::SourceRecord & source,
        afw::image::Exposure<PixelT> const & exposure,
        afw::geom::Point2D const & center
    ) const;

    LSST_MEAS_ALGORITHM_PRIVATE_INTERFACE(FitProfileAlgorithm);

    afw::table::Key< double > _fluxKey;
    afw::table::Key< double > _fluxErrKey;
    afw::table::Key< afw::table::Moments<float> > _ellipseKey;
    afw::table::Key< float > _chisqKey;
    afw::table::Key< afw::table::Flag > _flagKey;
    afw::table::Key< afw::table::Flag > _flagMaxIterKey;
    afw::table::Key< afw::table::Flag > _flagTinyStepKey;
    afw::table::Key< afw::table::Flag > _flagMinRadiusKey;
    afw::table::Key< afw::table::Flag > _flagMinAxisRatioKey;
    CONST_PTR(FitPsfControl) _psfCtrl;
};

inline PTR(FitProfileAlgorithm) FitProfileControl::makeAlgorithm(
    afw::table::Schema & schema,
    PTR(daf::base::PropertyList) const & metadata,
    algorithms::AlgorithmControlMap const & others
) const {
    return boost::static_pointer_cast<FitProfileAlgorithm>(_makeAlgorithm(schema, metadata, others));
}

}}}} // namespace lsst::meas::extensions::multiShapelet

#endif // !MULTISHAPELET_FitProfile_h_INCLUDED
