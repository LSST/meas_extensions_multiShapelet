// -*- lsst-c++ -*-
#ifndef LSST_MEAS_EXT_MULTISHAPELET_OPT_BoxConstraint_h_INCLUDED
#define LSST_MEAS_EXT_MULTISHAPELET_OPT_BoxConstraint_h_INCLUDED

#include "Eigen/Core"
#include "Eigen/Geometry"

namespace lsst { namespace meas { namespace extensions { namespace multiShapelet { namespace optimizer {

typedef Eigen::Array<bool,Eigen::Dynamic,1> ActiveSet;

class BoxConstraint {
public:
    
    typedef Eigen::AlignedBox<double,Eigen::Dynamic> EigenBox;

    explicit BoxConstraint(EigenBox const & b) : _box(b) {}

    explicit BoxConstraint(int dim) : _box(dim) {}

    EigenBox const & asEigen() const { return _box; }

    void set(int n, double min, double max);

    /**
     *  Return a bool array with true values in dimensions where x is on the constraint (within epsilon).
     *  The vector will be modified such that constrained dimensions are exactly on the constraint.
     */
    ActiveSet computeActiveSet(Eigen::VectorXd & x) const;

    /**
     *  Project a point onto the constraint.
     */
    void project(Eigen::VectorXd & x) const;

    /**
     *  Minimize q(x) = x^T G x / 2 + x^T c along the steepest descent direction.
     */
    void findCauchyPoint(Eigen::VectorXd & x, Eigen::MatrixXd const & G, Eigen::VectorXd const & c) const;

private:

    EigenBox _box;
};

}}}}} // namespace lsst::meas::extensions::multiShapelet::optimizer

#endif // !LSST_MEAS_EXT_MULTISHAPELET_OPT_BoxConstraint_h_INCLUDED
