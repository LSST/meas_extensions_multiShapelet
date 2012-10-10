#include "Eigen/Cholesky"
#include "Eigen/Eigenvalues"
#include "boost/make_shared.hpp"

#include "lsst/optimizer/SymmetricSolver.h"
#include "lsst/pex/exceptions.h"

namespace lsst { namespace optimizer {

namespace {

class EigenSolver : public SymmetricSolver {
public:

    virtual void factor(Eigen::MatrixXd const & matrix) {
        _solver.compute(matrix);
        if (_solver.info() != Eigen::Success) {
            throw LSST_EXCEPT(
                pex::exceptions::RuntimeErrorException,
                "Could not compute Eigensystem of matrix"
            );
        }
        _diag = _solver.eigenvalues().array().inverse();
        double cond = _threshold * _solver.eigenvalues()[_solver.eigenvalues().size() - 1];
        int n = 0;
        while (n < _diag.size() && _solver.eigenvalues()[n] < cond) ++n;
        _rank = _diag.size() - n;
        _diag.head(n).setZero();
    }

    virtual void solve(Eigen::VectorXd & x, Eigen::VectorXd const & rhs) const {
        _workspace.head(_rank) = _solver.eigenvectors().rightCols(_rank).adjoint() * rhs;
        _workspace.head(_rank).array() *= _diag.tail(_rank);
        x = _solver.eigenvectors().rightCols(_rank) * _workspace.head(_rank);
    }

    EigenSolver(int n, double threshold) :
        _rank(n), _threshold(threshold), _solver(n), _diag(n), _workspace(n)
        {}

private:
    int _rank;
    double _threshold;
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> _solver;
    Eigen::ArrayXd _diag;
    mutable Eigen::VectorXd _workspace;
};

class CholeskySolver : public SymmetricSolver {
public:

    virtual void factor(Eigen::MatrixXd const & matrix) {
        _solver.compute(matrix);
        if (_solver.info() != Eigen::Success) {
            throw LSST_EXCEPT(
                pex::exceptions::RuntimeErrorException,
                "Could not compute robust Cholesky decomposition of matrix"
            );
        }
    }

    virtual void solve(Eigen::VectorXd & x, Eigen::VectorXd const & rhs) const {
        x = _solver.solve(rhs);
    }

    explicit CholeskySolver(int n) : _solver(n) {}

private:
    Eigen::LDLT<Eigen::MatrixXd> _solver;
};

} // anonymous

double const SymmetricSolver::DEFAULT_THRESHOLD = std::sqrt(std::numeric_limits<double>::epsilon());

PTR(SymmetricSolver) SymmetricSolver::makeEigenSolver(int n, double threshold) {
    return boost::make_shared<EigenSolver>(n, threshold);
}

PTR(SymmetricSolver) SymmetricSolver::makeCholeskySolver(int n) {
    return boost::make_shared<CholeskySolver>(n);
}

}} // namespace lsst::optimizer
