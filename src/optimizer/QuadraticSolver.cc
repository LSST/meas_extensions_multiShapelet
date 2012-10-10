#include "lsst/optimizer/QuadraticSolver.h"

namespace lsst { namespace optimizer {

namespace {

double computeStepScale(Eigen::VectorXd const & y, Eigen::VectorXd const & py, double tau=1.0) {
    double alpha = 1.0;
    for (int i = 0; i < y.size(); ++i) {
        if (py[i] < 0.0) {
            alpha = std::min(alpha, -tau * y[i] / py[i]);
        }
    }
    return alpha;
}

} // anonymous

void QuadraticSolver::solve(
    Eigen::VectorXd & x,
    Eigen::MatrixXd const & G, Eigen::VectorXd const & c,
    Eigen::MatrixXd const & A, Eigen::VectorXd const & b
) const {
    int const n = G.rows();
    int const m = A.rows();
    _sym->factor(G);
    Eigen::VectorXd rhs(-c);
    _sym->solve(x, rhs);
    Eigen::VectorXd s = A * x - b;
    bool done = true;
    for (int i = 0; i < m; ++i) {
        if (s[i] < 0.0)
            done = false;
        s[i] = std::max(1.0, std::abs(s[i]));
    }
    if (done) return; // unconstrained solution obeys constraints; active set at solution is empty
    Eigen::VectorXd z = Eigen::VectorXd::Ones(m);
    Eigen::VectorXd px = Eigen::VectorXd::Zero(n);
    Eigen::VectorXd ps = Eigen::VectorXd::Zero(m);
    Eigen::VectorXd pz = Eigen::VectorXd::Zero(m);
    Eigen::VectorXd r1(n);
    Eigen::VectorXd r2(m);
    Eigen::VectorXd r3(m);
    double tau = _tauInitial;
    for (int k = 0; k < _maxIter; ++k) {
        assert((s.array() >= 0.0).all());
        assert((z.array() >= 0.0).all());
        r1 = A.adjoint() * z - (G * x + c);
        r2 = s - (A * x - b);
        r3 = -(s.array() * z.array()).matrix();
        computeStep(G, A, s, z, px, ps, pz, r1, r2, r3);
        double mu = s.dot(z) / m;
        double alpha_aff_z = computeStepScale(z, pz);
        double alpha_aff_s = computeStepScale(s, ps);
        double mu_aff = (z + alpha_aff_z * pz).dot(s + alpha_aff_s * ps) / m;
        double sigma = (mu_aff * mu_aff * mu_aff) / (mu * mu * mu);
        r3 += Eigen::VectorXd::Constant(m, sigma * mu);
        r3 -= (ps.array() * pz.array()).matrix();
        computeStep(G, A, s, z, px, ps, pz, r1, r2, r3);
        double alpha = std::min(computeStepScale(z, pz, tau), computeStepScale(s, ps, tau));
        tau *= _tauFactor;
        x += alpha * px;
        s += alpha * ps;
        z += alpha * pz;
    }
}

void QuadraticSolver::computeStep(
    Eigen::MatrixXd const & G, Eigen::MatrixXd const & A,
    Eigen::VectorXd const & s, Eigen::VectorXd const & z, 
    Eigen::VectorXd & px, Eigen::VectorXd & ps, Eigen::VectorXd & pz,
    Eigen::VectorXd const & r1, Eigen::VectorXd const & r2, Eigen::VectorXd const & r3
) const {
    if (_solveReduced) {
        Eigen::VectorXd sInv = s.array().inverse().matrix();
        Eigen::MatrixXd H = G + A.adjoint() * (z.array() * sInv.array()).matrix().asDiagonal() * A;
        _sym->factor(H);
        Eigen::VectorXd rhs = r1 + A.adjoint() * sInv.asDiagonal() *
            (r3 + z.asDiagonal() * r2);
        _sym->solve(px, rhs);
        ps = A * px - r2;
        pz = sInv.asDiagonal() * (r3 - z.asDiagonal() * ps);
    } else {
        int const n = G.rows();
        int const m = A.rows();
        Eigen::MatrixXd H(n + 2*m, n + 2*m);
        H <<
            G,              Eigen::MatrixXd::Zero(n,m),               -A.adjoint(),
            A,         -Eigen::MatrixXd::Identity(m,m), Eigen::MatrixXd::Zero(m,m),
            Eigen::MatrixXd::Zero(m,n), z.asDiagonal(),             s.asDiagonal();
        Eigen::VectorXd rhs(n + 2*m);
        rhs <<
            r1,
            r2,
            r3;
        _lu.compute(H);
        rhs = _lu.solve(rhs);
        px = rhs.segment(0, n);
        ps = rhs.segment(n, m);
        pz = rhs.segment(n+m, m);
    }
}

}} // namespace lsst::optimizer
