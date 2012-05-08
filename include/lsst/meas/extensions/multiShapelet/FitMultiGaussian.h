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
#ifndef MULTISHAPELET_FitMultiGaussian_h_INCLUDED
#define MULTISHAPELET_FitMultiGaussian_h_INCLUDED

#include "lsst/meas/extensions/multiShapelet/FitPsf.h"

namespace lsst { namespace meas { namespace extensions { namespace multiShapelet {

class FitMultiGaussianAlgorithm;

class MultiGaussianObjective;

class FitMultiGaussianControl : public algorithms::AlgorithmControl {
public:

    LSST_CONTROL_FIELD(profile, std::string, "Name of a registered multi-Gaussian profile.");
    LSST_CONTROL_FIELD(psfName, std::string, "Root name of the FitPsfAlgorithm fields.");
    LSST_CONTROL_FIELD(useShapeletPsfTerms, bool,
                       "If true, use the full double-shapelet PSF model; if false, use double-Gaussian.");
    LSST_CONTROL_FIELD(deconvolveShape, bool,
                       "Attempt to approximately deconvolve the canonical shape before "
                       "using it to set the initial parameters.");
    LSST_CONTROL_FIELD(initialRadiusFactor, double,
                       "How to scale the initial ellipse from whatever the canonical shape "
                       "measures to the radius defined by the profile (usually half-light radius");
    LSST_CONTROL_FIELD(usePixelWeights, bool, 
                       "If true, individually weigh pixels using the variance image.");
    LSST_CONTROL_FIELD(badMaskPlanes, std::vector<std::string>,
                       "Mask planes that indicate pixels that should be ignored in the fit.");
    LSST_CONTROL_FIELD(growFootprint, int, 
                       "Number of pixels to grow the footprint by.");

    PTR(FitMultiGaussianControl) clone() const {
        return boost::static_pointer_cast<FitMultiGaussianControl>(_clone());
    }

    PTR(FitMultiGaussianAlgorithm) makeAlgorithm(
        afw::table::Schema & schema,
        PTR(daf::base::PropertyList) const & metadata = PTR(daf::base::PropertyList)()
    ) const;
    
    FitMultiGaussianControl() : 
        algorithms::AlgorithmControl("multishapelet.exp", 2.5), 
        profile("tractor-exponential"), psfName("multishapelet.psf"),
        useShapeletPsfTerms(true), deconvolveShape(true), initialRadiusFactor(1.0),
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
        afw::table::Schema & schema, PTR(daf::base::PropertyList) const & metadata
    ) const;
};

/**
 *  @brief An elliptical model composed of several Gaussians.
 *
 *  All of the Gaussian components have the same ellipticity, and their amplitudes and radius ratios
 *  are fixed relative to each other.  Different combinations of amplitude and radius ratios are
 *  registered with the MultiGaussianRegistry class, and are usually approximations to Sersic functions
 *  with the ellipse defined at the half-light radius of the exact Sersic function being approximated.
 */
struct FitMultiGaussianModel {

    std::string profile;
    double amplitude;
    afw::geom::ellipses::Quadrupole ellipse;
    bool failed;  ///< set to true if the measurement failed

    FitMultiGaussianModel(
        FitMultiGaussianControl const & ctrl,
        double amplitude,
        ndarray::Array<double const,1,1> const & parameters
    );

    /// @brief Construct by extracting saved values from a SourceRecord.
    FitMultiGaussianModel(FitMultiGaussianControl const & ctrl, afw::table::SourceRecord const & source);

    /// @brief Deep copy constructor.
    FitMultiGaussianModel(FitMultiGaussianModel const & other);

    /// @brief Deep assignment operator.
    FitMultiGaussianModel & operator=(FitMultiGaussianModel const & other);

    /**
     *  @brief Return a MultiShapeletFunction representation of the model (unconvolved).
     */
    shapelet::MultiShapeletFunction asMultiShapelet(
        afw::geom::Point2D const & center = afw::geom::Point2D()
    ) const;

};

class FitMultiGaussianAlgorithm : public algorithms::Algorithm {
public:

    typedef FitMultiGaussianControl Control;
    typedef FitMultiGaussianModel Model;

    /// @brief Construct an algorithm instance and register its fields with a Schema.
    FitMultiGaussianAlgorithm(FitMultiGaussianControl const & ctrl, afw::table::Schema & schema);

    /// @brief Return the control object
    FitMultiGaussianControl const & getControl() const {
        return static_cast<FitMultiGaussianControl const &>(algorithms::Algorithm::getControl());
    }

    /**
     *  @brief Return an Objective that can be used to fit the convolved model to an image.
     *
     *  The objective function will use only the Gaussian terms in the PSF for the convolution.
     *
     *  This is provided primarily for testing and debugging purposes.
     */
    template <typename PixelT>
    static PTR(MultiGaussianObjective) makeObjective(
        FitMultiGaussianControl const & ctrl,
        FitPsfModel const & psfModel,
        afw::geom::ellipses::Quadrupole const & shape,
        ModelInputHandler const & inputs
    );

    /**
     *  @brief Return an optimizer that can be used to fit the convolved model image.
     *
     *  The optimizer's objective function will use only the Gaussian terms in the PSF for the convolution.
     *
     *  This is provided primarily for testing and debugging purposes; the user can create an optimizer,
     *  step through it, and use the FitMultiGaussianModel constructor that takes a parameter vector to
     *  visualize its progress.
     */
    template <typename PixelT>
    static HybridOptimizer makeOptimizer(
        FitPsfControl const & ctrl,
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
    template <typename PixelT>
    static void fitShapeletTerms(
        FitMultiGaussianControl const & ctrl,
        FitPsfModel const & psfModel,
        ModelInputHandler const & inputs,
        FitMultiGaussianModel & model
    );

    /**
     *  @brief Fit to a local PSF or kernel image.
     *
     *  @param[in] ctrl           Details of the model to fit.
     *  @param[in] psfModel       Localized double-shapelet PSF model.
     *  @param[in] shape          Shape measurement used to set initial ellipse parameters
     *                            (interpreted as defined by ctrl data members).
     *  @param[in] inputs         Inputs that determine the data to be fit.
     */
    template <typename PixelT>
    static FitMultiGaussianModel apply(
        FitMultiGaussianControl const & ctrl,
        FitPsfModel const & psfModel,
        afw::geom::ellipses::Quadrupole const & shape,
        ModelInputHandler const & inputs
    );

    /**
     *  @brief Fit to a local PSF or kernel image.
     *
     *  @param[in] ctrl           Details of the model to fit.
     *  @param[in] psfModel       Localized double-shapelet PSF model.
     *  @param[in] shape          Shape measurement used to set initial ellipse parameters
     *                            (interpreted as defined by ctrl data members).
     *  @param[in] footprint      Region of the image to fit (after modification according to
     *                            ctrl parameters).
     *  @param[in] image          Full masked image to fit.
     *  @param[in] center         Center of the object in the image's PARENT coordinate system
     *                            (i.e. xy0 is used).
     */
    template <typename PixelT>
    static FitMultiGaussianModel apply(
        FitMultiGaussianControl const & ctrl,
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

    LSST_MEAS_ALGORITHM_PRIVATE_INTERFACE(FitMultiGaussianAlgorithm);

    afw::table::KeyTuple<afw::table::Flux> _fluxKeys;
    afw::table::Key< afw::table::Moments<float> > _ellipseKey;
};

inline PTR(FitMultiGaussianAlgorithm) FitMultiGaussianControl::makeAlgorithm(
    afw::table::Schema & schema,
    PTR(daf::base::PropertyList) const & metadata
) const {
    return boost::static_pointer_cast<FitMultiGaussianAlgorithm>(_makeAlgorithm(schema, metadata));
}

}}}} // namespace lsst::meas::extensions::multiShapelet

#endif // !MULTISHAPELET_FitMultiGaussian_h_INCLUDED
