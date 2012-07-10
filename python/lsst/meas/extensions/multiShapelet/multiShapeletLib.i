// -*- lsst-c++ -*-
%define meas_extensions_multiShapelet_DOCSTRING
"
Measurement algorithms using galaxy models built from multi-scale shapelets.
"
%enddef

%feature("autodoc", "1");
%module(package="lsst.meas.extensions.multiShapelet", 
        docstring=meas_extensions_multiShapelet_DOCSTRING) multiShapeletLib

%{
#include "lsst/afw/geom.h"
#include "lsst/afw/geom/ellipses.h"
#include "lsst/afw/detection.h"
#include "lsst/pex/logging.h"
#include "lsst/meas/algorithms.h"
#include "lsst/shapelet.h"
#include "lsst/meas/extensions/multiShapelet.h"
%}

%include "lsst/p_lsstSwig.i"

%{
#define PY_ARRAY_UNIQUE_SYMBOL LSST_AFW_MATH_SHAPELETS_NUMPY_ARRAY_API
#include "numpy/arrayobject.h"
#include "ndarray/swig.h"
#include "ndarray/swig/eigen.h"
%}

%init %{
    import_array();
%}

%include "ndarray.i"
%include "std_vector.i"
%include "std_pair.i"

%lsst_exceptions()
%import "lsst/afw/geom/geomLib.i"
%import "lsst/afw/geom/ellipses/ellipsesLib.i"
%import "lsst/afw/detection/detectionLib.i"
%import "lsst/shapelet/shapeletLib.i"
%import "lsst/meas/algorithms/algorithmsLib.i"

%template(BoolPair) std::pair<bool,bool>;

%declareNumPyConverters(ndarray::Array<double const,1>);
%declareNumPyConverters(ndarray::Array<double const,1,1>);
%declareNumPyConverters(ndarray::Array<double,1,1>);
%declareNumPyConverters(ndarray::Array<double,2,-1>);
%declareNumPyConverters(ndarray::Array<double,2,-2>);
%declareNumPyConverters(Eigen::Matrix<double,3,Eigen::Dynamic>);

%include "lsst/meas/extensions/multiShapelet/ModelInputHandler.h"

%template(ModelInputHandler) lsst::meas::extensions::multiShapelet::ModelInputHandler::ModelInputHandler<float>;
%template(ModelInputHandler) lsst::meas::extensions::multiShapelet::ModelInputHandler::ModelInputHandler<double>;

%rename(__len__) lsst::meas::extensions::multiShapelet::MultiGaussian::size;
%rename(__getitem__) lsst::meas::extensions::multiShapelet::MultiGaussian::operator[];
%include "lsst/meas/extensions/multiShapelet/MultiGaussian.h"
%extend lsst::meas::extensions::multiShapelet::MultiGaussian {
%pythoncode %{
def __iter__(self):
    for i in xrange(len(self)):
        yield self[i]
%}
}

%include "lsst/meas/extensions/multiShapelet/GaussianModelBuilder.h"

%shared_ptr(lsst::meas::extensions::multiShapelet::Objective);
%include "lsst/meas/extensions/multiShapelet/HybridOptimizer.h"

%include "lsst/meas/extensions/multiShapelet/MultiGaussianRegistry.h"

%shared_ptr(lsst::meas::extensions::multiShapelet::MultiGaussianObjective);
%returnCopy(lsst::meas::extensions::multiShapelet::MultiGaussianObjective::getInputs);
%include "lsst/meas/extensions/multiShapelet/MultiGaussianObjective.h"

namespace lsst { namespace meas { namespace extensions { namespace multiShapelet {

%extend Objective {
%pythoncode %{
    def cast(self, type): return type._cast(self)
%}
}

%pythoncode %{
    import lsst.afw.geom.ellipses
%}

%extend MultiGaussianObjective {
    static PTR(MultiGaussianObjective) _cast(PTR(Objective) const & input) {
        return boost::dynamic_pointer_cast<lsst::meas::extensions::multiShapelet::MultiGaussianObjective>(input);
    }
    %pythoncode %{
        EllipseCore = lsst.afw.geom.ellipses.SeparableConformalShearLogTraceRadius
    %}
}

}}}} // namespace lsst::meas::extensions::multiShapelet

%shared_ptr(lsst::meas::extensions::multiShapelet::FitPsfControl);
%shared_ptr(lsst::meas::extensions::multiShapelet::FitPsfAlgorithm);
%include "lsst/meas/extensions/multiShapelet/FitPsf.h"

%shared_ptr(lsst::meas::extensions::multiShapelet::FitProfileControl);
%shared_ptr(lsst::meas::extensions::multiShapelet::FitProfileAlgorithm);
%include "lsst/meas/extensions/multiShapelet/FitProfile.h"

%template(adjustInputs) lsst::meas::extensions::multiShapelet::FitProfileAlgorithm::adjustInputs<float>;
%template(adjustInputs) lsst::meas::extensions::multiShapelet::FitProfileAlgorithm::adjustInputs<double>;

%shared_ptr(lsst::meas::extensions::multiShapelet::FitComboControl);
%shared_ptr(lsst::meas::extensions::multiShapelet::FitComboAlgorithm);
%include "lsst/meas/extensions/multiShapelet/FitCombo.h"

%template(adjustInputs) lsst::meas::extensions::multiShapelet::FitComboAlgorithm::adjustInputs<float>;
%template(adjustInputs) lsst::meas::extensions::multiShapelet::FitComboAlgorithm::adjustInputs<double>;
