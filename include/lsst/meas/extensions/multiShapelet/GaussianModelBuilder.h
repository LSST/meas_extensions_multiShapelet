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
#ifndef MULTISHAPELET_GaussianModelBuilder_h_INCLUDED
#define MULTISHAPELET_GaussianModelBuilder_h_INCLUDED

#include "ndarray/eigen.h"
#include "lsst/afw/geom/ellipses.h"
#include "lsst/afw/detection/Footprint.h"

#include "lsst/meas/extensions/multiShapelet/EllipseSquaredNorm.h"

namespace lsst { namespace meas { namespace extensions { namespace multiShapelet {

class GaussianModelBuilder {
public:

    GaussianModelBuilder(
        ndarray::Array<double const,1,1> const & x,
        ndarray::Array<double const,1,1> const & y
    );
    
    GaussianModelBuilder(GaussianModelBuilder const & other);

    GaussianModelBuilder & operator=(GaussianModelBuilder const & other);

    int getSize() const { return _x.size(); }

    void update(afw::geom::ellipses::BaseCore const & core);

    ndarray::Array<double const,1,1> computeModel();

    ndarray::Array<double const,1,1> getModel() const { return _model; }

    void computeDerivative(
        ndarray::Array<double,2,-1> const & output,
        Eigen::Matrix<double,3,Eigen::Dynamic> const & jacobian,
        bool add = false
    );

    void setOutput(ndarray::Array<double,1,1> const & array);

private:

    EllipseSquaredNorm _esn;
    Eigen::Matrix<double,3,3> _esnJacobian;
    ndarray::EigenView<double const,1,1> _x;
    ndarray::EigenView<double const,1,1> _y;
    Eigen::VectorXd _rx;
    Eigen::VectorXd _ry;
    ndarray::Array<double,1,1> _model;
};


}}}} // namespace lsst::meas::extensions::multiShapelet

#endif // !MULTISHAPELET_GaussianModelBuilder_h_INCLUDED
