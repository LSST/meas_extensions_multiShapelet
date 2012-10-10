// -*- lsst-c++ -*-
#ifndef LSST_OPTIMIZER_SymmetricSolver_h_INCLUDED
#define LSST_OPTIMIZER_SymmetricSolver_h_INCLUDED

#include "Eigen/Core"
#include "lsst/base.h"

namespace lsst { namespace optimizer {

class SymmetricSolver {
public:

    static double const DEFAULT_THRESHOLD; // == sqrt(epsilon)

    virtual void factor(Eigen::MatrixXd const & matrix) = 0;
    virtual void solve(Eigen::VectorXd & x, Eigen::VectorXd const & rhs) const = 0;
    void solve(Eigen::VectorXd & x) const { solve(x, x); }
    
    static PTR(SymmetricSolver) makeEigenSolver(int n, double threshold=DEFAULT_THRESHOLD);
    static PTR(SymmetricSolver) makeCholeskySolver(int n);
};

}} // namespace lsst::optimizer

#endif // !LSST_OPTIMIZER_SymmetricSolver_h_INCLUDED
