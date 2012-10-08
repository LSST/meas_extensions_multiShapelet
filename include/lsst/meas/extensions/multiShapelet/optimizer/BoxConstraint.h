// -*- lsst-c++ -*-
#ifndef LSST_MEAS_EXT_MULTISHAPELET_OPT_BoxConstraint_h_INCLUDED
#define LSST_MEAS_EXT_MULTISHAPELET_OPT_BoxConstraint_h_INCLUDED

#include "Eigen/Core"
#include "Eigen/Geometry"

namespace lsst { namespace meas { namespace extension { namespace multiShapelet { namespace optimizer {

class BoxConstraint {
public:
    
    explicit BoxConstraint(int dim);

    void set(int n, double min, double max);

    BoxConstraint intersect(BoxConstraint const & other) const;

    /**
     *  Project a point onto the constraint.
     */
    void project(Eigen::VectorXd & x) const;

    /**
     *  Minimize q(x) = x^T G x / 2 + x^T c along the steepest descent direction.
     */
    void findCauchyPoint(Eigen::VectorXd & x, Eigen::MatrixXd const & G, Eigen::VectorXd const & c) const;

private:

    typedef Eigen::AlignedBox<double,Eigen::Dynamic> Box;

    explicit BoxConstraint(Box const & b) : _box(b) {}

    Box _box;
};

}}}}} // namespace lsst::meas::extension::multiShapelet::optimizer

#endif // !LSST_MEAS_EXT_MULTISHAPELET_OPT_BoxConstraint_h_INCLUDED
