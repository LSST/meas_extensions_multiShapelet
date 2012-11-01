#ifdef HAVE_NLOPT
#include "lsst/optimizer/NLopt.h"

namespace lsst { namespace optimizer {

namespace {

// Function whose pointer we pass to NLopt; just casts and calls the Objective in func_data.
double callObjective(
    unsigned n, double const * x, double * gradient, void * func_data
) {
    Objective const * func = reinterpret_cast<Objective*>(func_data);
    int const p = func->getParameterSize();
    assert(int(n) == p);
    ndarray::Array<double const,1,1> xv = ndarray::external(x, ndarray::makeVector(p));
    ndarray::Array<double,1,1> gv = ndarray::external(gradient, ndarray::makeVector(p));
    return func->evaluate(xv, gv);
}

// Function whose pointer we pass to NLopt; just casts and calls the Constraint in func_data.
void callConstraint(
    unsigned m, double * result, unsigned n, double const * x, double * gradient, void * func_data
) {
    Constraint const * func = reinterpret_cast<Constraint*>(func_data);
    int const p = func->getParameterSize();
    int const c = func->getConstraintSize();
    assert(int(n) == p);
    assert(int(m) == c);
    ndarray::Array<double const,1,1> xv = ndarray::external(x, ndarray::makeVector(p));
    ndarray::Array<double,1,1> rv = ndarray::external(result, ndarray::makeVector(c));
    ndarray::Array<double,2,2> gv = ndarray::external(gradient, ndarray::makeVector(c,p));
    func->evaluate(xv, rv, gv);
}

} // anonymous

NLopt::NLopt(AlgorithmEnum algorithm, int n) :
    _obj(), _iqc(), _eqc(),
    _opt(::nlopt_create(algorithm, n))
{
    if (!_opt) {
        throw LSST_EXCEPT(
            pex::exceptions::RuntimeErrorException,
            "Unknown error initializing NLopt optimizer"
        );
    }
}

NLopt::NLopt(NLopt const & other) : 
    _obj(other._obj), _iqc(other._iqc), _eqc(other._eqc), _opt(::nlopt_copy(other._opt))
{
    if (!_opt) {
        throw LSST_EXCEPT(
            pex::exceptions::RuntimeErrorException,
            "Unknown error copying NLopt optimizer"
        );
    }
}

NLopt & NLopt::operator=(NLopt const & other) {
    if (&other != this) {
        NLopt tmp(other);
        swap(tmp);
    }
    return *this;
}

void NLopt::swap(NLopt & other) {
    std::swap(other._opt, _opt);
    _obj.swap(other._obj);
    _iqc.swap(other._iqc);
    _eqc.swap(other._eqc);
}

NLopt::~NLopt() {
    ::nlopt_destroy(_opt);
}

std::pair<NLopt::ResultEnum,double> NLopt::optimize(ndarray::Array<double,1,1> const & x) {
    assert(x.getSize<0>() == getParameterSize());
    double fval = 0;
    ResultEnum result = ::nlopt_optimize(_opt, x.getData(), &fval);
    return std::make_pair(result, fval);
}

void NLopt::setObjective(CONST_PTR(Objective) const & func, bool maximize) {
    assert(func->getParameterSize() == getParameterSize());
    _obj = func;
    if (maximize) {
        ::nlopt_set_max_objective(_opt, callObjective, const_cast<Objective*>(_obj.get()));
    } else {
        ::nlopt_set_min_objective(_opt, callObjective, const_cast<Objective*>(_obj.get()));
    }
}

void NLopt::unsetLowerBounds() {
    ::nlopt_set_lower_bounds1(_opt, -HUGE_VAL);
}
void NLopt::unsetUpperBounds() {
    ::nlopt_set_upper_bounds1(_opt, HUGE_VAL);
}
void NLopt::setLowerBounds(double lower) {
    ::nlopt_set_lower_bounds1(_opt, lower);
}
void NLopt::setUpperBounds(double upper) {
    ::nlopt_set_upper_bounds1(_opt, upper);
}
void NLopt::setLowerBounds(ndarray::Array<double const,1,1> const & lower) {
    assert(lower.size() == getParameterSize());
    ::nlopt_set_lower_bounds(_opt, lower.getData());
}
void NLopt::setUpperBounds(ndarray::Array<double const,1,1> const & upper) {
    assert(upper.size() == getParameterSize());
    ::nlopt_set_upper_bounds(_opt, upper.getData());
}

void NLopt::addInequalityConstraint(CONST_PTR(Objective) const & func, double tolerance) {
    assert(func->getParameterSize() == getParameterSize());
    _iqc.push_back(func);
    ::nlopt_add_inequality_constraint(_opt, callObjective, const_cast<Objective*>(func.get()), tolerance);
}
void NLopt::addEqualityConstraint(CONST_PTR(Objective) const & func, double tolerance) {
    assert(func->getParameterSize() == getParameterSize());
    _eqc.push_back(func);
    ::nlopt_add_equality_constraint(_opt, callObjective, const_cast<Objective*>(func.get()), tolerance);
}
void NLopt::addInequalityConstraints(
    CONST_PTR(Constraint) const & func,
    ndarray::Array<double const,1,1> const & tolerance
) {
    assert(func->getParameterSize() == getParameterSize());
    assert(tolerance.getSize<0>() == getParameterSize());
    _iqc.push_back(func);
    ::nlopt_add_inequality_mconstraint(
        _opt, func->getConstraintSize(),
        callConstraint,
        const_cast<Constraint*>(func.get()),
        tolerance.getData()
    );
}

void NLopt::addEqualityConstraints(
    CONST_PTR(Constraint) const & func,
    ndarray::Array<double const,1,1> const & tolerance
) {
    assert(func->getParameterSize() == getParameterSize());
    assert(tolerance.getSize<0>() == getParameterSize());
    _eqc.push_back(func);
    ::nlopt_add_equality_mconstraint(
        _opt, func->getConstraintSize(),
        callConstraint,
        const_cast<Constraint*>(func.get()),
        tolerance.getData()
    );
}

void NLopt::removeInequalityConstraints() {
    _iqc.clear();
    ::nlopt_remove_inequality_constraints(_opt);
}
void NLopt::removeEqualityConstraints() {
    _eqc.clear();
    ::nlopt_remove_equality_constraints(_opt);
}

void NLopt::setStopValue(double value) {
    if (::nlopt_set_stopval(_opt, value) < 0) {
        throw LSST_EXCEPT(pex::exceptions::RuntimeErrorException, "Error setting NLopt stopvalue");
    }
}

void NLopt::setObjectiveRelativeTolerance(double value) {
    if (::nlopt_set_ftol_rel(_opt, value) < 0) {
        throw LSST_EXCEPT(pex::exceptions::RuntimeErrorException, "Error setting NLopt ftol_rel");
    }
}

void NLopt::setObjectiveAbsoluteTolerance(double value) {
    if (::nlopt_set_ftol_abs(_opt, value) < 0) {
        throw LSST_EXCEPT(pex::exceptions::RuntimeErrorException, "Error setting NLopt ftol_abs");
    }
}

void NLopt::setParameterRelativeTolerance(double value) {
    if (::nlopt_set_xtol_rel(_opt, value) < 0) {
        throw LSST_EXCEPT(pex::exceptions::RuntimeErrorException, "Error setting NLopt xtol_rel");
    }
}

void NLopt::setParameterAbsoluteTolerance(double value) {
    if (::nlopt_set_xtol_abs1(_opt, value) < 0) {
        throw LSST_EXCEPT(pex::exceptions::RuntimeErrorException, "Error setting NLopt xtol_abs1");
    }
}

void NLopt::setParameterAbsoluteTolerance(ndarray::Array<double const,1,1> const & values) {
    assert(values.getSize<0>() == getParameterSize());
    if (::nlopt_set_xtol_abs(_opt, values.getData()) < 0) {
        throw LSST_EXCEPT(pex::exceptions::RuntimeErrorException, "Error setting NLopt xtol_abs1");
    }
}

void NLopt::getParameterRelativeTolerance(ndarray::Array<double,1,1> const & values) const {
    assert(values.getSize<0>() == getParameterSize());
    if (::nlopt_get_xtol_abs(_opt, values.getData()) < 0) {
        throw LSST_EXCEPT(pex::exceptions::RuntimeErrorException, "Error getting NLopt xtol_abs");
    }    
}

void NLopt::setMaxEvaluations(int value) {
    if (::nlopt_set_maxeval(_opt, value) < 0) {
        throw LSST_EXCEPT(pex::exceptions::RuntimeErrorException, "Error setting NLopt maxeval");
    }
}

void NLopt::setMaxTime(double value) {
    if (::nlopt_set_maxtime(_opt, value) < 0) {
        throw LSST_EXCEPT(pex::exceptions::RuntimeErrorException, "Error setting NLopt maxtime");
    }
}

void NLopt::setInitialStep(ndarray::Array<double const,1,1> const & dx) {
    assert(dx.getSize<0>() == getParameterSize());
    if (::nlopt_set_initial_step(_opt, dx.getData()) < 0) {
        throw LSST_EXCEPT(pex::exceptions::RuntimeErrorException, "Error setting NLopt initial step");
    }
}

void NLopt::setInitialStep(double dx) {
    if (::nlopt_set_initial_step1(_opt, dx) < 0) {
        throw LSST_EXCEPT(pex::exceptions::RuntimeErrorException, "Error setting NLopt initial step");
    }
}

void NLopt::unsetInitialStep() {
    if (::nlopt_set_initial_step(_opt, NULL) < 0) {
        throw LSST_EXCEPT(pex::exceptions::RuntimeErrorException, "Error setting NLopt initial step");
    }
}

}} // lsst::optimizer

#endif // HAVE_NLOPT
