// -*- lsst-c++ -*-
#ifndef LSST_OPTIMIZER_QuadraticSolver_h_INCLUDED
#define LSST_OPTIMIZER_QuadraticSolver_h_INCLUDED

#include "Eigen/Core"
#include "Eigen/LU"
#include "lsst/optimizer/SymmetricSolver.h"

namespace lsst { namespace optimizer {

class QuadraticSolver {
public:

    /**
     *  @brief Minimize @f$ \frac{1}{2}x^T G x + c^T x @f$ subject to @f$ A x - b \ge 0 @f$.
     */
    void solve(
        Eigen::VectorXd & x,
        Eigen::MatrixXd const & G, Eigen::VectorXd const & c,
        Eigen::MatrixXd const & A, Eigen::VectorXd const & b
    ) const;

private:

    void computeStep(
        Eigen::MatrixXd const & G, Eigen::MatrixXd const & A,
        Eigen::VectorXd const & s, Eigen::VectorXd const & z, 
        Eigen::VectorXd & px, Eigen::VectorXd & ps, Eigen::VectorXd & pz,
        Eigen::VectorXd const & r1, Eigen::VectorXd const & r2, Eigen::VectorXd const & r3
    ) const;

    bool _solveReduced;
    int _maxIter;
    double _tauInitial;
    double _tauFactor;
    mutable Eigen::FullPivLU<Eigen::MatrixXd> _lu;
    PTR(SymmetricSolver) _sym;
};

}} // namespace lsst::optimizer

#endif // !LSST_OPTIMIZER_QuadraticSolver_h_INCLUDED
