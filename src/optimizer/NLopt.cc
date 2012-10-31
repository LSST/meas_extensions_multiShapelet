#include "lsst/optimizer/NLopt.h"

namespace lsst { namespace optimizer {

namespace {

// Function whose pointer we pass to NLopt; just casts and calls the Objective in func_data.
double callObjective(unsigned n, double const * x, double * gradient, void * func_data) {
    Objective const * obj = reinterpret_cast<Objective*>(func_data);
    Eigen::VectorXd xv = Eigen::Map<Eigen::VectorXd const>(x, obj->getParameterSize());
    if (gradient) {
        Eigen::Map<Eigen::VectorXd> gm(gradient, obj->getParameterSize());
        Eigen::VectorXd gv = gm;
        double r = obj->evaluate(xv, gv);
        gm = gv;
        return r;
    } else {
        return obj->evaluate(xv);
    }
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

NLopt::AlgorithmEnum NLopt::getAlgorithm() const { return ::nlopt_get_algorithm(_opt); }

int NLopt::getParameterSize() const { return ::nlopt_get_dimension(_opt); }

std::pair<NLopt::ResultEnum,double> NLopt::optimize(Eigen::VectorXd & x) {
    assert(x.size() == getParameterSize());
    double fval = 0;
    ResultEnum result = ::nlopt_optimize(_opt, x.data(), &fval);
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

void NLopt::setBounds(Eigen::AlignedBoxXd const & box) {
    setLowerBounds(box.min());
    setUpperBounds(box.max());
}
void NLopt::unsetBounds() {
    unsetLowerBounds();
    unsetUpperBounds();
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
void NLopt::setLowerBounds(Eigen::VectorXd const & lower) {
    assert(lower.size() == getParameterSize());
    ::nlopt_set_lower_bounds(_opt, lower.data());
}
void NLopt::setUpperBounds(Eigen::VectorXd const & upper) {
    assert(upper.size() == getParameterSize());
    ::nlopt_set_upper_bounds(_opt, upper.data());
}

void NLopt::addInequalityConstraint(CONST_PTR(Objective) const & iqc, double tolerance) {
    assert(iqc->getParameterSize() == getParameterSize());
    _iqc.push_back(iqc);
    ::nlopt_add_inequality_constraint(_opt, callObjective, const_cast<Objective*>(iqc.get()), tolerance);
}
void NLopt::addEqualityConstraint(CONST_PTR(Objective) const & eqc, double tolerance) {
    assert(eqc->getParameterSize() == getParameterSize());
    _eqc.push_back(eqc);
    ::nlopt_add_equality_constraint(_opt, callObjective, const_cast<Objective*>(eqc.get()), tolerance);
}
void NLopt::removeInequalityConstraints() {
    _iqc.clear();
    ::nlopt_remove_inequality_constraints(_opt);
}
void NLopt::removeEqualityConstraints() {
    _eqc.clear();
    ::nlopt_remove_equality_constraints(_opt);
}

}} // lsst::optimizer
