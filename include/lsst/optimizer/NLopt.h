// -*- lsst-c++ -*-
#ifndef LSST_OPTIMIZER_NLopt_h_INCLUDED
#define LSST_OPTIMIZER_NLopt_h_INCLUDED

#include "boost/noncopyable.hpp"
#include "nlopt.hp"

namespace lsst { namespace optimizer {

class Objective {
public:

    virtual double evaluate(Eigen::VectorXd const & x) const = 0;
    virtual double evaluate(Eigen::VectorXd const & x, Eigen::VectorXd & gradient) const = 0;

    virtual ~Objective() {}

};

class Constraint {
public:

    int getSize() const { return _size; }

    virtual void evaluate(Eigen::VectorXd const & x, Eigen::VectorXd & c) const = 0;
    virtual void evaluate(Eigen::VectorXd const & x, Eigen::VectorXd & c, Eigen::MatrixXd & gc) const = 0;

    virtual ~Constraint() {}

protected:

    explicit Constraint(int size) : _size(size) {}

    Constraint(Constraint const & other) : _size(other._size) {}

    void operator=(Constraint const & other) { _size = other._size; }

private:
    int _size;
};

class NLopt {
public:
    
    typedef ::nlopt_algorithm AlgorithmEnum;
    typedef ::nlopt_result ResultEnum;

    NLopt(AlgorithmEnum algorithm, int n);

    NLopt(NLopt const & other);

    NLopt & operator=(NLopt const & other);

    ResultEnum optimize(Eigen::VectorXd const & x);

    void setObjective(CONST_PTR(Objective) const & obj);

    void setBounds(Eigen::AlignedBox2d const & box);
    void setBounds(Eigen::VectorXd const & lower, Eigen::VectorXd const & upper);
    void setBounds(double lower, Eigen::VectorXd const & upper);
    void setBounds(Eigen::VectorXd const & lower, double upper);
    void setLowerBounds();
    void setUpperBounds();
    void setLowerBounds(double lower);
    void setUpperBounds(double upper);
    void setLowerBounds(Eigen::VectorXd const & lower);
    void setUpperBounds(Eigen::VectorXd const & upper);

    void setInequalityConstraint(CONST_PTR(Constraint) const & iqc);
    void setEqualityConstraint(CONST_PTR(Constraint) const & eqc);

    ~NLopt();

private:
    CONST_PTR(Objective) _obj;
    CONST_PTR(Constraint) _iqc;
    CONST_PTR(Constraint) _eqc;
    ::nlopt_opt _opt;
};

}} // namespace lsst::optimizer

#endif // !LSST_OPTIMIZER_NLopt_h_INCLUDED
