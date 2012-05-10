#include "lsst/meas/extensions/multiShapelet/MultiGaussian.h"

namespace lsst { namespace meas { namespace extensions { namespace multiShapelet {

shapelet::ShapeletFunction MultiGaussianComponent::makeShapelet(
    afw::geom::ellipses::Ellipse const & ellipse, int order
) const {
    static double const FACTOR = 1.0 / (shapelet::NORMALIZATION * shapelet::NORMALIZATION * 2.0);
    shapelet::ShapeletFunction result(order=0, shapelet::HERMITE, ellipse);
    result.getEllipse().getCore().scale(radius);
    result.getCoefficients()[0] = amplitude * FACTOR / result.getEllipse().getCore().getArea();
    return result;
}

void MultiGaussianComponent::readShapeletAmplitude(
    double coeff0, afw::geom::ellipses::BaseCore const & ellipse
) {
    static double const FACTOR = shapelet::NORMALIZATION * shapelet::NORMALIZATION * 2.0;
    double area = ellipse.getArea() * radius * radius;
    amplitude *= coeff0 * area * FACTOR;
}


double MultiGaussianComponent::integrate(List const & components) {
    double flux = 0.0;
    for (List::const_iterator i = components.begin(); i != components.end(); ++i) {
        flux += i->amplitude;
    }
    return flux;
}

}}}} // namespace lsst::meas::extensions::multiShapelet
