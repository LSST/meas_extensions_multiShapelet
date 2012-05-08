// -*- lsst-c++ -*-
/* 
 * LSST Data Management System
 * Copyright 2012 LSST Corporation.
 * 
 * This product includes software developed by the
 * LSST Project (http://www.lsst.org/).
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the LSST License Statement and 
 * the GNU General Public License along with this program.  If not, 
 * see <http://www.lsstcorp.org/LegalNotices/>.
 */
#ifndef MULTISHAPELET_MultiGaussianObjective_h_INCLUDED
#define MULTISHAPELET_MultiGaussianObjective_h_INCLUDED

#include "lsst/meas/extensions/multiShapelet/MultiGaussian.h"
#include "lsst/meas/extensions/multiShapelet/HybridOptimizer.h"
#include "lsst/meas/extensions/multiShapelet/ModelInputHandler.h"
#include "lsst/meas/extensions/multiShapelet/GaussianModelBuilder.h"

namespace lsst { namespace meas { namespace extensions { namespace multiShapelet {

class MultiGaussianObjective : public Objective {
public:

    typedef afw::geom::ellipses::SeparableConformalShearTraceRadius EllipseCore;

    virtual void computeFunction(
        ndarray::Array<double const,1,1> const & parameters, 
        ndarray::Array<double,1,1> const & function
    );

    virtual void computeDerivative(
        ndarray::Array<double const,1,1> const & parameters, 
        ndarray::Array<double const,1,1> const & function,
        ndarray::Array<double,2,-2> const & derivative
    );

    double getAmplitude() const { return _amplitude; }
    
    ModelInputHandler const & getInputs() const { return _inputs; }

#ifndef SWIG // Don't need to create this class from Python, just use it.

    MultiGaussianObjective(
        ModelInputHandler const & inputs,
        MultiGaussianList const & components
    );

    MultiGaussianObjective(
        ModelInputHandler const & inputs,
        MultiGaussianList const & components,
        MultiGaussianList const & psfComponents,
        afw::geom::ellipses::Quadrupole const & psfEllipse
    );

#endif

private:

    typedef std::vector<GaussianModelBuilder> BuilderList;

    double _amplitude;
    double _modelSquaredNorm;
    EllipseCore _ellipse;
    EllipseCore _psfEllipse;
    MultiGaussianList _components;
    MultiGaussianList _psfComponents;
    ModelInputHandler _inputs;
    BuilderList _builders;
    ndarray::Array<double,1,1> _model;
};

}}}} // namespace lsst::meas::extensions::multiShapelet

#endif // !MULTISHAPELET_MultiGaussianObjective_h_INCLUDED
