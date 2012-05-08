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

#include "lsst/meas/extensions/multiShapelet/GaussianModelBuilder.h"

namespace lsst { namespace meas { namespace extensions { namespace multiShapelet {

GaussianModelBuilder::GaussianModelBuilder(
    ndarray::Array<double const,1,1> const & x,
    ndarray::Array<double const,1,1> const & y
) :
    _x(x), _y(y), _rx(x.size()), _ry(y.size())
{
    if (_x.size() != _y.size()) {
        throw LSST_EXCEPT(
            pex::exceptions::LengthErrorException,
            (boost::format("x coordinate array size (%d) does not match y coordinate array size (%d)")
             % _x.size() % _y.size()).str()
        );
    }
}

GaussianModelBuilder::GaussianModelBuilder(GaussianModelBuilder const & other) :
    _x(other._x), _y(other._y), _rx(other._rx), _ry(other._ry)
{}

GaussianModelBuilder & GaussianModelBuilder::operator=(GaussianModelBuilder const & other) {
    if (&other != this) {
        _x.reset(other._x.shallow());
        _y.reset(other._y.shallow());
        _rx = other._rx;
        _ry = other._ry;
        _esn = other._esn;
        _esnJacobian = other._esnJacobian;
        if (!other._model.isEmpty()) _model = ndarray::copy(other._model);
    }
    return *this;
}

void GaussianModelBuilder::update(afw::geom::ellipses::BaseCore const & core) {
    _esnJacobian = _esn.update(core);
}

ndarray::Array<double const,1,1> GaussianModelBuilder::computeModel() {
    if (_model.isEmpty()) {
        _model = ndarray::allocate(_x.size());
    }
    ndarray::EigenView<double,1,1> z(_model);
    _esn(_x, _y, _rx, _ry, z);
    z.array() = std::exp(-0.5 * z.array());
    return _model;
}

void GaussianModelBuilder::computeDerivative(
    ndarray::Array<double,2,-1> const & output,
    Eigen::Matrix<double,3,Eigen::Dynamic> const & jacobian,
    bool add
) {
    if (_model.isEmpty()) {
        throw LSST_EXCEPT(
            pex::exceptions::LogicErrorException,
            "computeDerivative called before computeModel"
        );
    }
    if (output.getSize<0>() != _x.size()) {
        throw LSST_EXCEPT(
            pex::exceptions::InvalidParameterException,
            (boost::format("Incorrect number of rows for array: got %d, expected %d")
             % output.getSize<0>() % _x.size()).str()
        );
    }
    if (output.getSize<0>() != _x.size()) {
        throw LSST_EXCEPT(
            pex::exceptions::InvalidParameterException,
            (boost::format("Mismatch between array (%d) and Jacobian dimensions (%d)")
             % output.getSize<1>() % jacobian.cols()).str()
        );
    }
    ndarray::EigenView<double,2,-1> out(output);
    if (!add) out.setZero();
    // Compute derivative wrt ellipse core
    Eigen::MatrixXd dz_de = Eigen::MatrixXd::Zero(_x.size(), jacobian.cols());
    Eigen::Matrix<double,3,Eigen::Dynamic> coreJacobian = _esnJacobian * jacobian;
    _esn.dEllipse(_x, _y, _rx, _ry, coreJacobian, dz_de);
    for (int n = 0; n < jacobian.cols(); ++n) {
        dz_de.col(n).array() *= -0.5 * _model.asEigen<Eigen::ArrayXpr>();
    }
    out += dz_de;
}

void GaussianModelBuilder::setOutput(ndarray::Array<double,1,1> const & array) {
    if (array.getSize<0>() != _x.size()) {
        throw LSST_EXCEPT(
            pex::exceptions::InvalidParameterException,
            (boost::format("Incorrect size for array: got %d, expected %d")
             % array.getSize<0>() % _x.size()).str()
        );
    }
}

}}}} // namespace lsst::meas::extensions::multiShapelet
