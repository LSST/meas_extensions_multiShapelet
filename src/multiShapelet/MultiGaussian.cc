#include "lsst/meas/extensions/multiShapelet/MultiGaussian.h"

namespace lsst { namespace meas { namespace extensions { namespace multiShapelet {

shapelet::ShapeletFunction GaussianComponent::makeShapelet(
    afw::geom::ellipses::Ellipse const & ellipse, int order
) const {
    static double const FACTOR = 1.0 / (shapelet::NORMALIZATION * shapelet::NORMALIZATION * 2.0);
    shapelet::ShapeletFunction result(order=0, shapelet::HERMITE, ellipse);
    result.getEllipse().getCore().scale(radius);
    result.getCoefficients()[0] = flux * FACTOR / result.getEllipse().getCore().getArea();
    return result;
}

void GaussianComponent::readShapeletAmplitude(
    double coeff0, afw::geom::ellipses::BaseCore const & ellipse
) {
    static double const FACTOR = shapelet::NORMALIZATION * shapelet::NORMALIZATION * 2.0;
    double area = ellipse.getArea() * radius * radius;
    flux *= coeff0 * area * FACTOR;
}


double MultiGaussian::integrate() const {
    double flux = 0.0;
    for (List::const_iterator i = _components.begin(); i != _components.end(); ++i) {
        flux += i->flux;
    }
    return flux;
}

afw::geom::ellipses::Quadrupole MultiGaussian::deconvolve(
    afw::geom::ellipses::Quadrupole const & fullMoments,
    afw::geom::ellipses::Quadrupole const & psfMoments,
    MultiGaussian const & psfComponents
) const {
    static const double EPS = std::numeric_limits<double>::epsilon();
    //
    // We treat fullMoments as unweighted, infinite, zero-noise moments,
    // and match that to the analytic unweighted, infinite, zero-noise moments
    // of the profile convolved with the PSF.  This is analogous to just
    // subtracting Gaussian moments, but with some weighted sums involved.
    //
    // We might be able to do even better by matching the profile's moments
    // as weighted by a Gaussian with the decorrected measured moments, but
    // that's a nonlinear problem that would require a numerical root finder.
    // But it may still be worth considering if it saves us fitting
    // iterations later.
    //
    double lhs = 0.0;
    Eigen::Matrix2d rhs = Eigen::Matrix2d::Zero();
    Eigen::Matrix2d w = fullMoments.getMatrix();
    Eigen::Matrix2d p = psfMoments.getMatrix();
    for (List::const_iterator i = psfComponents.begin(); i != psfComponents.end(); ++i) {
        for (List::const_iterator j = this->begin(); j != this->end(); ++j) {
            double ab = i->flux * j->flux;
            rhs += ab * (w - p * i->radius * i->radius);
            lhs += ab * j->radius * j->radius;
        }
    }
    afw::geom::ellipses::Quadrupole::Matrix q(rhs / lhs);
    if (q(0,0) <= EPS || q(1,1) <= EPS || q.determinant() <= EPS) q = Eigen::Matrix2d::Identity() * EPS;
    return afw::geom::ellipses::Quadrupole(q);
}


}}}} // namespace lsst::meas::extensions::multiShapelet
