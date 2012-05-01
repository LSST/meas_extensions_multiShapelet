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

    LSST_CONTROL_FIELD("profile", std::string, "Name of a registered multi-Gaussian profile.");


    PTR(FitMultiGaussianControl) clone() const {
        return boost::static_pointer_cast<FitMultiGaussianControl>(_clone());
    }

    PTR(FitMultiGaussianAlgorithm) makeAlgorithm(
        afw::table::Schema & schema,
        PTR(daf::base::PropertyList) const & metadata = PTR(daf::base::PropertyList)()
    ) const;
    
    FitMultiGaussianControl() : 
        algorithms::AlgorithmControl("multishapelet.exp", 2.5), 
        profile("tractor-exponential")
    {}

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

    double amplitude;
    afw::geom::ellipses::Quadrupole ellipse; ///< ellipse corresponding to inner expansion
    bool failed;  ///< set to true if the measurement failed

    /**
     *  @brief Construct a model from a double-Gaussian optimization parameter vector.
     *
     *  This is used internally by FitPsfAlgorithm to construct the model after fitting
     *  an elliptical double-Gaussian, and the 4-element parameter vector should correspond
     *  to the parameters used by FitPsfAlgorithm::makeObjective and
     *  FitPsfAlgorithm::makeOptimizer.  By using those functions and this constructor,
     *  we can inspect the optimizer at each step in Python.
     */
    FitMultiGaussianModel(
        FitMultiGaussianControl const & ctrl, double amplitude,
        ndarray::Array<double const,1,1> const & parameters
    );

    /// @brief Construct by extracting saved values from a SourceRecord.
    FitMultiGaussianModel(FitMultiGaussianControl const & ctrl, afw::table::SourceRecord const & source);

    /// @brief Deep copy constructor.
    FitMultiGaussianModel(FitMultiGaussianControl const & other);

    /// @brief Deep assignment operator.
    FitMultiGaussianModel & operator=(FitMultiGaussianModel const & other);

    /**
     *  @brief Return a MultiShapeletFunction representation of the model (unconvolved).
     *
     *  The elements will be [inner, outer].
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
    static PTR(MultiGaussianObjective) makeObjective(
        FitMultiGaussianControl const & ctrl,
        FitPsfModel const & psfModel,
        afw::image::Image<double> const & image,
        afw::geom::Point2D const & center
    );

    /**
     *  @brief Return an optimizer that can be used to fit an elliptical double-Gaussian to the image.
     *
     *  This is provided primarily for testing and debugging purposes; the user can create an optimizer,
     *  step through it, and use the FitPsfModel constructor that takes a parameter vector to
     *  visualize its progress.
     */
    static HybridOptimizer makeOptimizer(
        FitPsfControl const & ctrl,
        afw::image::Image<double> const & image,
        afw::geom::Point2D const & center
    );

    /**
     *  @brief Given a double-Gaussian-only model, fit additional shapelet terms.
     *
     *  This is provided primarily for testing and debugging purposes; it is the second
     *  part of apply(), after the nonlinear double-Gaussian fit.
     */
    static void fitShapeletTerms(
        FitPsfControl const & ctrl,
        afw::image::Image<double> const & image,
        afw::geom::Point2D const & center,
        FitPsfModel & model
    );

    /**
     *  @brief Fit to a local PSF or kernel image.
     *
     *  This overload accepts the configuration options (inner and outer order) as separate
     *  values, and hence does not require a control object or an Algorithm instance.
     *
     *  @param[in] ctrl           Details of the model to fit.
     *  @param[in] image          Postage-stamp image of the PSF.
     *  @param[in] center         Center of the PSF in the image's PARENT coordinate system
     *                            (i.e. xy0 is used).
     */
    static FitPsfModel apply(
        FitPsfControl const & ctrl,
        afw::image::Image<double> const & image,
        afw::geom::Point2D const & center
    );

    /**
     *  @brief Fit a PSF object evaluated at a point.
     *
     *  This overload accepts the configuration options (inner and outer order) as separate
     *  values, and hence does not require a control object or an Algorithm instance.
     *
     *  @param[in] ctrl           Details of the model to fit.
     *  @param[in] psf            PSF object
     *  @param[in] center         Point at which to evaluate the PSF.
     */
    static FitPsfModel apply(
        FitPsfControl const & ctrl,
        afw::detection::Psf const & psf,
        afw::geom::Point2D const & center
    ) {
        return apply(ctrl, *psf.computeImage(center), center);
    }

private:

    template <typename PixelT>
    void _apply(
        afw::table::SourceRecord & source,
        afw::image::Exposure<PixelT> const & exposure,
        afw::geom::Point2D const & center
    ) const;

    LSST_MEAS_ALGORITHM_PRIVATE_INTERFACE(FitPsfAlgorithm);

    afw::table::Key< afw::table::Array<float> > _innerKey;
    afw::table::Key< afw::table::Array<float> > _outerKey;
    afw::table::Key< afw::table::Moments<float> > _ellipseKey;
    afw::table::Key< afw::table::Flag > _flagKey;
};

inline PTR(FitPsfAlgorithm) FitPsfControl::makeAlgorithm(
    afw::table::Schema & schema,
    PTR(daf::base::PropertyList) const & metadata
) const {
    return boost::static_pointer_cast<FitPsfAlgorithm>(_makeAlgorithm(schema, metadata));
}

}}}} // namespace lsst::meas::extensions::multiShapelet

#endif // !MULTISHAPELET_FitMultiGaussian_h_INCLUDED
