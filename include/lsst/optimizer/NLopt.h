// -*- lsst-c++ -*-
#ifndef LSST_OPTIMIZER_NLopt_h_INCLUDED
#define LSST_OPTIMIZER_NLopt_h_INCLUDED

#include <vector>

#include "boost/noncopyable.hpp"
#include "nlopt.hp"

namespace lsst { namespace optimizer {

class Objective {
public:

    int getParameterSize() const { return _parameterSize; }

    virtual double evaluate(Eigen::VectorXd const & x) const = 0;
    virtual double evaluate(Eigen::VectorXd const & x, Eigen::VectorXd & gradient) const = 0;

    virtual ~Objective() {}

protected:

    explicit Objective(int parameterSize) : _parameterSize(parameterSize) {}

    Objective(Objective const & other) : _parameterSize(other._parameterSize) {}

private:
    void operator=(Objective const & other); // intentionally undefined
    int _parameterSize;
};

class Constraint {
public:

    int getParameterSize() const { return _parameterSize; }
    int getConstraintSize() const { return _constraintSize; }

    virtual void evaluate(Eigen::VectorXd const & x, Eigen::VectorXd & c) const = 0;
    virtual void evaluate(Eigen::VectorXd const & x, Eigen::VectorXd & c, Eigen::MatrixXd & gc) const = 0;

    virtual ~Constraint() {}

protected:

    explicit Constraint(int parameterSize, int constraintSize) :
        _parameterSize(parameterSize), _constraintSize(constraintSize)
    {}

    Constraint(Constraint const & other) :
        _parameterSize(other._parameterSize), _constraintSize(other._constraintSize)
    {}

private:
    void operator=(Constraint const & other);

    int _parameterSize;
    int _constraintSize;
};

class NLopt {
public:
    
    typedef ::nlopt_algorithm AlgorithmEnum;
    typedef ::nlopt_result ResultEnum;

    NLopt(AlgorithmEnum algorithm, int n);

    NLopt(NLopt const & other);

    NLopt & operator=(NLopt const & other);

    void swap(NLopt & other);

    ~NLopt();

    AlgorithmEnum getAlgorithm() const;

    int getParameterSize() const;

    std::pair<ResultEnum,double> optimize(Eigen::VectorXd & x);

    void setObjective(CONST_PTR(Objective) const & func, bool maximize=false);

    void setBounds(Eigen::AlignedBoxXd const & box);
    void unsetLowerBounds();
    void unsetUpperBounds();
    void setLowerBounds(double lower);
    void setUpperBounds(double upper);
    void setLowerBounds(Eigen::VectorXd const & lower);
    void setUpperBounds(Eigen::VectorXd const & upper);

    void addInequalityConstraint(CONST_PTR(Objective) const & func, double tolerance);
    void addEqualityConstraint(CONST_PTR(Objective) const & func, double tolerance);
    void removeInequalityConstraints();
    void removeEqualityConstraints();

private:
    typedef std::vector< CONST_PTR(Objective)> > ConstraintVector;

    CONST_PTR(Objective) _obj;
    ConstraintVector _iqc;
    ConstraintVector _eqc;
    ::nlopt_opt _opt;
};

}} // namespace lsst::optimizer

#endif // !LSST_OPTIMIZER_NLopt_h_INCLUDED
