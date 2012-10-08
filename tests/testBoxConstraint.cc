#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE BoxConstraint
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wunused-variable"
#include "boost/test/unit_test.hpp"
#pragma clang diagnostic pop
#include "boost/format.hpp"

#include "lsst/meas/extensions/multiShapelet/optimizer/BoxConstraint.h"
#include "Eigen/Core"

namespace opt = lsst::meas::extensions::multiShapelet::optimizer;

namespace {

opt::BoxConstraint makeRandomBox(int n) {
    Eigen::VectorXd a = Eigen::VectorXd::Random(n);
    Eigen::VectorXd b = Eigen::VectorXd::Random(n);
    opt::BoxConstraint r(n);
    for (int i = 0; i < n; ++i) {
        if (a[i] < -0.5) a[i] = -std::numeric_limits<double>::infinity();
        if (b[i] > 0.5) b[i] = std::numeric_limits<double>::infinity();
        r.set(i, std::min(a[i],b[i]), std::max(a[i],b[i]));
    }
    return r;
}

void doProjectTest(opt::BoxConstraint const & a) {
    static int const N_MC_POINTS = 100;
    Eigen::MatrixXd p = Eigen::MatrixXd::Random(a.asEigen().min().size(), N_MC_POINTS);
    for (int i = 0; i < N_MC_POINTS; ++i) {
        Eigen::VectorXd q = p.col(i);
        a.project(q);
        if (a.asEigen().contains(p.col(i))) {
            BOOST_CHECK(p.col(i).isApprox(q));
        } else {
            BOOST_CHECK(p.col(i) != q);
        }
        BOOST_CHECK(a.asEigen().contains(q));
    }
}

void doCauchyPointTest(
    opt::BoxConstraint const & box, Eigen::VectorXd const & x,
    Eigen::MatrixXd const & G, Eigen::VectorXd const & c
) {
    BOOST_CHECK(box.asEigen().contains(x)); 
    Eigen::VectorXd y(x);
    box.findCauchyPoint(y, G, c);
    BOOST_CHECK(box.asEigen().contains(y));
    Eigen::VectorXd g = G * x + c;
    double qx = c.dot(x) + 0.5 * x.dot(G * x);
    double qy = c.dot(y) + 0.5 * y.dot(G * y);
    BOOST_CHECK((qy < qx) || y.isApprox(x));    
}

} // anonymous

BOOST_AUTO_TEST_CASE(BoxConstraint) {
    static int N_BOXES = 100;
    for (int i = 0; i < N_BOXES; ++i) {
        opt::BoxConstraint b = makeRandomBox(4);
        doProjectTest(b);
    }
}
