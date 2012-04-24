#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE EllipseSquaredNorm
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wunused-variable"
#include "boost/test/unit_test.hpp"
#pragma clang diagnostic pop
#include "boost/format.hpp"

#include "lsst/base.h"
#include "lsst/afw/geom/ellipses.h"
#include "lsst/meas/extensions/multiShapelet/EllipseSquaredNorm.h"
#include "Eigen/Core"
#include <vector>

namespace el = lsst::afw::geom::ellipses;
namespace ms = lsst::meas::extensions::multiShapelet;

static double const eps = 1E-4;

template <typename Function>
Eigen::Matrix<double,Function::N,1>
computeDerivative(Function f, Eigen::Matrix<double,Function::N,1> const & initial) {
    Eigen::Matrix<double,Function::N,1> x(initial);
    Eigen::Matrix<double,Function::N,1> result;
    for (int n = 0; n < Function::N; ++n) {
        double h = (x[n] + eps) * eps;
        x[n] += h;
        result[n] = f(x);
        x[n] -= 2.0 * h;
        result[n] -= f(x);
        x[n] += h;
        result[n] /= (2.0 * h);
    }
    return result;
}

class CoordFunc {
public:

    static int const N = 2;

    double operator()(Eigen::Matrix<double,N,1> const & p) const {
        double rx, ry, z;
        esn(p.x(), p.y(), rx, ry, z);
        return z;
    }

    ms::EllipseSquaredNorm esn;
};

class EllipseFunc {
public:

    static int const N = 3;

    double operator()(Eigen::Matrix<double,N,1> const & p) {
        double rx, ry, z;
        core->readParameters(p.data());
        esn.update(*core, false);
        esn(x, y, rx, ry, z);
        return z;
    }

    ms::EllipseSquaredNorm esn;
    PTR(el::BaseCore) core;
    double x;
    double y;
};

typedef std::vector<PTR(el::BaseCore)> EllipseVector;
typedef EllipseVector::const_iterator EllipseIter;

EllipseVector makeTestEllipses() {
    EllipseVector result;
    result.push_back(boost::make_shared<el::Quadrupole>());
    result.push_back(boost::make_shared<el::SeparableDistortionDeterminantRadius>());
    result.push_back(boost::make_shared<el::SeparableConformalShearLogTraceRadius>());
    result.push_back(boost::make_shared<el::SeparableReducedShearTraceRadius>());
    result.push_back(boost::make_shared<el::Quadrupole>(3.0, 4.2, 0.25));
    result.push_back(boost::make_shared<el::SeparableDistortionDeterminantRadius>(0.5, 0.3, 1.5));
    result.push_back(boost::make_shared<el::SeparableConformalShearLogTraceRadius>(0.6, 0.1, 1.7));
    result.push_back(boost::make_shared<el::SeparableReducedShearTraceRadius>(0.0, 1E-4, 2.0));
    return result;
}

BOOST_AUTO_TEST_CASE(Vectorization) {
    // Test that the vectorized versions of functions do the same things as the scalar ones.
    ms::EllipseSquaredNorm esn;
    EllipseVector const ellipses = makeTestEllipses();
    int const n = 50;
    Eigen::VectorXd x = Eigen::VectorXd::Random(n);
    Eigen::VectorXd y = Eigen::VectorXd::Random(n);
    x[0] = 0.0; // want to make sure we don't divide by zero anywhere.
    y[0] = 0.0;
    Eigen::VectorXd rx = Eigen::VectorXd::Zero(n);
    Eigen::VectorXd ry = Eigen::VectorXd::Zero(n);
    Eigen::VectorXd rxc = Eigen::VectorXd::Zero(n);
    Eigen::VectorXd ryc = Eigen::VectorXd::Zero(n);
    Eigen::VectorXd z = Eigen::VectorXd::Zero(n);
    Eigen::VectorXd dz_dx = Eigen::VectorXd::Zero(n);
    Eigen::VectorXd dz_dy = Eigen::VectorXd::Zero(n);
    Eigen::Matrix<double,Eigen::Dynamic,3> dz_de = Eigen::Matrix<double,Eigen::Dynamic,3>::Zero(n, 3);
    for (EllipseIter i = ellipses.begin(); i != ellipses.end(); ++i) {
        el::BaseCore::Jacobian jacobian = esn.update(**i, true);
        esn(x, y, rx, ry, z);
        for (int j = 0; j < n; ++j) {
            double rxj=0., ryj=0., zj=0.;
            esn(x[j], y[j], rxj, ryj, zj);
            BOOST_CHECK_CLOSE(rxj, rx[j], 1E-15);
            BOOST_CHECK_CLOSE(ryj, ry[j], 1E-15);
            BOOST_CHECK_CLOSE(zj, z[j], 1E-15);
        }
        esn.dCoords(x, y, rx, ry, dz_dx, dz_dy);
        for (int j = 0; j < n; ++j) {
            double dz_dxj=0., dz_dyj=0.;
            esn.dCoords(x[j], y[j], rx[j], ry[j], dz_dxj, dz_dyj);
            BOOST_CHECK_CLOSE(dz_dxj, dz_dx[j], 1E-15);
            BOOST_CHECK_CLOSE(dz_dyj, dz_dy[j], 1E-15);
        }
        rxc = rx;
        ryc = ry;
        dz_de.setZero();
        esn.dEllipse(x, y, rxc, ryc, jacobian, dz_de);
        for (int j = 0; j < n; ++j) {
            Eigen::Vector3d dz_dej = Eigen::Vector3d::Zero();
            esn.dEllipse(x[j], y[j], rx[j], ry[j], jacobian, dz_dej);
            for (int k = 0; k < 3; ++k) {
                BOOST_CHECK_CLOSE(dz_dej[k], dz_de(j, k), 1E-15);
            }
        }
    }
}

BOOST_AUTO_TEST_CASE(Functor) {
    ms::EllipseSquaredNorm esn;
    EllipseVector const ellipses = makeTestEllipses();
    int const n = 50;
    Eigen::VectorXd x = Eigen::VectorXd::Random(n);
    Eigen::VectorXd y = Eigen::VectorXd::Random(n);
    x[0] = 0.0; // want to make sure we don't divide by zero anywhere.
    y[0] = 0.0;
    for (EllipseIter i = ellipses.begin(); i != ellipses.end(); ++i) {
        esn.update(**i, false);
        el::Quadrupole quadrupole(**i);
        Eigen::Matrix2d qInv = quadrupole.getMatrix().inverse();
        for (int j = 0; j < n; ++j) {
            double rx=0., ry=0., z=0.;
            esn(x[j], y[j], rx, ry, z);
            Eigen::Vector2d xy(x[j], y[j]);
            double zc = xy.dot(qInv * xy);
            BOOST_CHECK_CLOSE(z, zc, 1E-13);
        }
    }
}

BOOST_AUTO_TEST_CASE(dCoords) {
    ms::EllipseSquaredNorm esn;
    EllipseVector const ellipses = makeTestEllipses();
    int const n = 50;
    Eigen::VectorXd x = Eigen::VectorXd::Random(n);
    Eigen::VectorXd y = Eigen::VectorXd::Random(n);
    x[0] = 0.0; // want to make sure we don't divide by zero anywhere.
    y[0] = 0.0;
    for (EllipseIter i = ellipses.begin(); i != ellipses.end(); ++i) {
        el::BaseCore::Jacobian jacobian = esn.update(**i, true);
        CoordFunc func = { esn };
        for (int j = 0; j < n; ++j) {
            double rx=0., ry=0., z=0.;
            double dz_dx=0., dz_dy=0.;
            esn(x[j], y[j], rx, ry, z);
            esn.dCoords(x[j], y[j], rx, ry, dz_dx, dz_dy);
            Eigen::Vector2d dzn = computeDerivative(func, Eigen::Vector2d(x[j], y[j]));
            BOOST_CHECK_CLOSE(dz_dx, dzn.x(), 1E-4);
            BOOST_CHECK_CLOSE(dz_dy, dzn.y(), 1E-4);
        }
    }
}

BOOST_AUTO_TEST_CASE(dEllipse) {
    ms::EllipseSquaredNorm esn;
    EllipseVector const ellipses = makeTestEllipses();
    int const n = 4;
    Eigen::VectorXd x = Eigen::VectorXd::Random(n);
    Eigen::VectorXd y = Eigen::VectorXd::Random(n);
    x[0] = 0.0; // want to make sure we don't divide by zero anywhere.
    y[0] = 0.0;
    for (EllipseIter i = ellipses.begin(); i != ellipses.end(); ++i) {
        el::BaseCore::Jacobian jacobian = esn.update(**i, true);
        EllipseFunc func = { esn, (**i).clone(), 0.0, 0.0 };
        for (int j = 0; j < n; ++j) {
            func.x = x[j];
            func.y = y[j];
            double rx=0., ry=0., z=0.;
            Eigen::Vector3d dza = Eigen::Vector3d::Zero();
            esn(x[j], y[j], rx, ry, z);
            esn.dEllipse(x[j], y[j], rx, ry, jacobian, dza);
            Eigen::Vector3d dzn = computeDerivative(func, (**i).getParameterVector());
            BOOST_CHECK_CLOSE(dza[0], dzn[0], 1E-4);
            BOOST_CHECK_CLOSE(dza[1], dzn[1], 1E-4);
            BOOST_CHECK_CLOSE(dza[2], dzn[2], 1E-4);
        }
    }
}
