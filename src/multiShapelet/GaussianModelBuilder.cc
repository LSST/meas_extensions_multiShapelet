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

GaussianModelBuilder::GaussianModelBuilder(afw::detection::Footprint const & region) :
    _x(region.getArea()), _y(region.getArea()),
    _tx(region.getArea()), _ty(region.getArea()),
    _rx(region.getArea()), _ry(region.getArea())
{
    int n = 0;
    for (
        afw::detection::Footprint::SpanList::const_iterator i = region.getSpans().begin();
        i != region.getSpans().end();
        ++i
    ) {
        for (int x = (**i).getX0(); x <= (**i).getX1(); ++x, ++n) {
            _x[n] = x;
            _y[n] = (**i).getY();
        }
    }
}

GaussianModelBuilder::GaussianModelBuilder(afw::geom::Box2I const & region) :
    _x(region.getArea()), _y(region.getArea()),
    _tx(region.getArea()), _ty(region.getArea()),
    _rx(region.getArea()), _ry(region.getArea())
{
    int n = 0;
    afw::geom::Point2I const llc = region.getMin();
    afw::geom::Point2I const urc = region.getMax();
    for (int y = llc.getY(); y <= urc.getY(); ++y) {
        for (int x = llc.getX(); x <= urc.getX(); ++x, ++n) {
            _x[n] = x;
            _y[n] = y;
        }
    }
}

GaussianModelBuilder::GaussianModelBuilder(GaussianModelBuilder const & other) :
    _x(other._x), _y(other._y), _tx(other._tx), _ty(other._ty), _rx(other._rx), _ry(other._ry)
{}

GaussianModelBuilder & GaussianModelBuilder::operator=(GaussianModelBuilder const & other) {
    if (&other != this) {
        _x = other._x;
        _y = other._y;
        _tx = other._tx;
        _ty = other._ty;
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

void GaussianModelBuilder::update(afw::geom::ellipses::Ellipse const & ellipse) {
    update(ellipse.getCore());
    update(ellipse.getCenter());
}

void GaussianModelBuilder::update(afw::geom::Point2D const & center) {
    _tx.array() = _x.array() - center.getX();
    _ty.array() = _y.array() - center.getY();
}

ndarray::Array<double const,1,1> GaussianModelBuilder::computeModel() {
    if (_model.isEmpty()) {
        _model = ndarray::allocate(_x.size());
    }
    ndarray::EigenView<double,1,1> z(_model);
    _esn(_tx, _ty, _rx, _ry, z);
    z.array() = std::exp(-0.5 * z.array());
    return _model;
}

void GaussianModelBuilder::computeDerivative(
    ndarray::Array<double,2,-1> const & output,
    Eigen::Matrix<double,5,Eigen::Dynamic> const & jacobian,
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
    double const eps = std::numeric_limits<double>::epsilon() * jacobian.lpNorm<Eigen::Infinity>();
    if ((jacobian.bottomRows<2>().array().abs() > eps).any()) {
        // Compute derivative wrt center point
        Eigen::VectorXd dz_dx = Eigen::VectorXd::Zero(_x.size());
        Eigen::VectorXd dz_dy = Eigen::VectorXd::Zero(_y.size());
        _esn.dCoords(_tx, _ty, _rx, _ry, dz_dx, dz_dy);
        dz_dx.array() *= -0.5 * _model.asEigen<Eigen::ArrayXpr>();
        dz_dy.array() *= -0.5 * _model.asEigen<Eigen::ArrayXpr>();
        out -= dz_dx * jacobian.row(3);
        out -= dz_dy * jacobian.row(4);
    }
    if ((jacobian.topRows<3>().array().abs() > eps).any()) {
        // Compute derivative wrt ellipse core
        Eigen::MatrixXd dz_de = Eigen::MatrixXd::Zero(_x.size(), jacobian.cols());
        Eigen::Matrix<double,3,Eigen::Dynamic> coreJacobian = _esnJacobian * jacobian.topRows<3>();
        _esn.dEllipse(_tx, _ty, _rx, _ry, coreJacobian, dz_de);
        for (int n = 0; n < jacobian.cols(); ++n) {
            dz_de.col(n).array() *= -0.5 * _model.asEigen<Eigen::ArrayXpr>();
        }
        out += dz_de;
    }
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
