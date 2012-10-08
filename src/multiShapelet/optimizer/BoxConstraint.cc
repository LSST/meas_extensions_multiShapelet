#include "lsst/meas/extensions/multiShapelet/optimizer/BoxConstraint.h"
#include "lsst/utils/ieee.h"

#include <vector>

namespace lsst { namespace meas { namespace extension { namespace multiShapelet { namespace optimizer {

BoxConstraint::BoxConstraint(int dim) : _box(dim) {}

void BoxConstraint::set(int n, double min, double max) {
    assert(n < _box.dim());
    _box.min()[n] = min;
    _box.max()[n] = max;
}

BoxConstraint BoxConstraint::intersect(BoxConstraint const & other) const {
    return BoxConstraint(this->_box.intersection(other._box));
}
   
void BoxConstraint::project(Eigen::VectorXd & x) const {
    assert(x.size() == _box.dim());
    x = x.array().max(_box.min().array()).matrix();
    x = x.array().min(_box.max().array()).matrix();
}

namespace {

struct ZeroCompare {
    bool operator()(double a) const {
        return a < std::numeric_limits<double>::epsilon();
    }
};

struct EqualCompare {
    bool operator()(double a, double b) const {
        return std::abs(a - b) < std::numeric_limits<double>::epsilon();
    }
};

std::vector<double> computeBreaks(Eigen::VectorXd & t) {
    std::replace_if(t.data(), t.data() + t.size(), ZeroCompare(), 0.0);
    std::vector<double> r(t.data(), t.data() + t.size());
    std::sort(r.begin(), r.end());
    if (r.size() > 0u && r.front() == 0.0) { // exact == should be safe because of the replace_if
        r.erase(r.begin());
    }
    std::vector<double>::iterator end = std::unique(r.begin(), r.end(), EqualCompare());
    r.resize(end - r.begin());
    return r;
}

} // anonymous

void BoxConstraint::findCauchyPoint(
    Eigen::VectorXd & x, Eigen::MatrixXd const & G, Eigen::VectorXd const & c
) const {
    Eigen::VectorXd g = G * x + c;
    Eigen::VectorXd t = Eigen::VectorXd::Constant(x.size(), std::numeric_limits<double>::infinity());
    for (int i = 0; i < x.size(); ++i) {
        if (g[i] < 0.0 && !lsst::utils::isinf(_box.max()[i])) {
            t[i] = (x[i] - _box.max()[i]) / g[i];
        } else if (g[i] > 0.0 && !lsst::utils::isinf(-_box.min()[i])) {
            t[i] = (x[i] - _box.min()[i]) / g[i];
        }
    }
    std::vector<double> breaks = computeBreaks(t); // sort, remove duplicates and zeros
    double t0 = 0.0;
    for (std::size_t j = 0; j < breaks.size(); ++j) {
        Eigen::VectorXd p = (t.array() > t0).select(-g, Eigen::VectorXd::Zero(t.size()));
        Eigen::VectorXd Gx = G * x;
        double f1 = p.dot(c + Gx);
        double f2 = p.dot(G * p);
        double dtMax = breaks[j] - t[0];
        if (f1 > 0.0) return; // minimum is at current value of x.
        double dt = -f1 / f2;
        if (dt < dtMax) {
            x += p * dt;
            return;
        }
        x += dtMax * p;
        t0 = breaks[j];
    }
}

}}}}} // namespace lsst::meas::extension::multiShapelet::optimizer
