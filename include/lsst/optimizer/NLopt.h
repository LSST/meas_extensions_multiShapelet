// -*- lsst-c++ -*-
#ifndef LSST_OPTIMIZER_NLopt_h_INCLUDED
#define LSST_OPTIMIZER_NLopt_h_INCLUDED

#include <vector>

#include "boost/noncopyable.hpp"
#include "nlopt.h"

#include "ndarray.h"
#include "lsst/base.h"
#include "lsst/pex/exceptions.h"

namespace lsst { namespace optimizer {

class OptFunc {
public:
    virtual ~OptFunc() {}
};

class Objective : public OptFunc {
public:

    int getParameterSize() const { return _parameterSize; }

    virtual double evaluate(
        ndarray::Array<double const,1,1> const & x,
        ndarray::Array<double,1,1> const & gradient
    ) const = 0;

protected:

    explicit Objective(int parameterSize) : _parameterSize(parameterSize) {}

    Objective(Objective const & other) : _parameterSize(other._parameterSize) {}

private:
    void operator=(Objective const & other); // intentionally undefined
    int _parameterSize;
};

class Constraint : public OptFunc {
public:

    int getParameterSize() const { return _parameterSize; }
    int getConstraintSize() const { return _constraintSize; }

    virtual void evaluate(
        ndarray::Array<double const,1,1> const & x,
        ndarray::Array<double,1,1> const & c,
        ndarray::Array<double,2,2> const & gc
    ) const = 0;

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

    AlgorithmEnum getAlgorithm() const { return ::nlopt_get_algorithm(_opt); }

    int getParameterSize() const { return ::nlopt_get_dimension(_opt); }

    std::pair<ResultEnum,double> optimize(ndarray::Array<double,1,1> const & x);

    void setObjective(CONST_PTR(Objective) const & func, bool maximize=false);

    void setLowerBounds(double lower);
    void setUpperBounds(double upper);
    void setLowerBounds(ndarray::Array<double const,1,1> const & lower);
    void setUpperBounds(ndarray::Array<double const,1,1> const & upper);
    void unsetLowerBounds();
    void unsetUpperBounds();

    void addInequalityConstraint(CONST_PTR(Objective) const & func, double tolerance);
    void addEqualityConstraint(CONST_PTR(Objective) const & func, double tolerance);
    void addInequalityConstraints(CONST_PTR(Constraint) const & func,
                                  ndarray::Array<double const,1,1> const & tolerance);
    void addEqualityConstraints(CONST_PTR(Constraint) const & func,
                                ndarray::Array<double const,1,1> const & tolerance);
    void removeInequalityConstraints();
    void removeEqualityConstraints();

    /// Stop if the objective function is goes beyond the given value.
    void setStopValue(double value);
    double getStopValue() const { return ::nlopt_get_stopval(_opt); }

    /// Stop if the objective function changes by less than the given relative amount in a step.
    void setObjectiveRelativeTolerance(double value);
    double getObjectiveRelativeTolerance() const { return ::nlopt_get_ftol_rel(_opt); }

    /// Stop if the objective function changes by less than the given absolute amount in a step.
    void setObjectiveAbsoluteTolerance(double value);
    double getObjectiveAbsoluteTolerance() const { return ::nlopt_get_ftol_abs(_opt); }

    /// Stop if the largest parameter change is less than the given relative amount in a step.
    void setParameterRelativeTolerance(double value);
    double getParameterRelativeTolerance() const { return ::nlopt_get_xtol_rel(_opt); }

    /// Stop if the largest parameter changes is less than the given absolute amount(s) in a step.
    void setParameterAbsoluteTolerance(double value);
    void setParameterAbsoluteTolerance(ndarray::Array<double const,1,1> const & values);
    void getParameterRelativeTolerance(ndarray::Array<double,1,1> const & values) const;

    /// Stop after more than the given number of functions evaluations.
    void setMaxEvaluations(int value);
    int getMaxEvaluations() const { return ::nlopt_get_maxeval(_opt); }

    /// Stop after more than the given nubmer of seconds
    void setMaxTime(double value);
    int getMaxTime() const { return ::nlopt_get_maxtime(_opt); }

    /// Set the initial step.
    void setInitialStep(ndarray::Array<double const,1,1> const & dx);
    void setInitialStep(double dx);
    void unsetInitialStep();

private:
    typedef std::vector<CONST_PTR(OptFunc)> ConstraintVector;

    CONST_PTR(Objective) _obj;
    ConstraintVector _iqc;
    ConstraintVector _eqc;
    ::nlopt_opt _opt;
};

}} // namespace lsst::optimizer

#endif // !LSST_OPTIMIZER_NLopt_h_INCLUDED
