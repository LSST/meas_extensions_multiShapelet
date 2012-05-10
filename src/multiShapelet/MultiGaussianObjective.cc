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

#include "lsst/meas/extensions/multiShapelet/MultiGaussianObjective.h"

namespace lsst { namespace meas { namespace extensions { namespace multiShapelet {

MultiGaussianObjective::MultiGaussianObjective(
    ModelInputHandler const & inputs,
    MultiGaussianList const & components
) : Objective(inputs.getSize(), 3), _amplitude(1.0), _modelSquaredNorm(1.0),
    _ellipse(), _inputs(inputs), _model(ndarray::allocate(inputs.getSize()))
{
    _builders.reserve(components.size());
    for (MultiGaussianList::const_iterator i = components.begin(); i != components.end(); ++i) {
        _builders.push_back(GaussianModelBuilder(_inputs.getX(), _inputs.getY(), i->flux, i->radius));
    }
}

MultiGaussianObjective::MultiGaussianObjective(
    ModelInputHandler const & inputs,
    MultiGaussianList const & components,
    MultiGaussianList const & psfComponents,
    afw::geom::ellipses::Quadrupole const & psfEllipse
) : Objective(inputs.getSize(), 3), _amplitude(1.0), _modelSquaredNorm(1.0),
    _ellipse(), _inputs(inputs), _model(ndarray::allocate(inputs.getSize()))
{
    _builders.reserve(components.size() * psfComponents.size());
    for (MultiGaussianList::const_iterator j = psfComponents.begin(); j != psfComponents.end(); ++j) {
        afw::geom::ellipses::Quadrupole psfComponentEllipse(psfEllipse);
        psfComponentEllipse.scale(j->radius);
        for (MultiGaussianList::const_iterator i = components.begin(); i != components.end(); ++i) {
            _builders.push_back(
                GaussianModelBuilder(
                    _inputs.getX(), _inputs.getY(), i->flux, i->radius,
                    psfComponentEllipse, j->flux
                )
            );
        }
    }
}

void MultiGaussianObjective::computeFunction(
    ndarray::Array<double const,1,1> const & parameters, 
    ndarray::Array<double,1,1> const & function
) {
    _ellipse.readParameters(parameters.getData());
    ndarray::EigenView<double,1,1> model(_model);
    model.setZero();
    for (std::size_t n = 0; n < _builders.size(); ++n) {
        _builders[n].update(_ellipse);
        model += _builders[n].getModel().asEigen();
    }
    if (!_inputs.getWeights().isEmpty()) {
        model.array() *= _inputs.getWeights().asEigen<Eigen::ArrayXpr>();
    }
    _modelSquaredNorm = model.squaredNorm();
    _amplitude = model.dot(_inputs.getData().asEigen()) / _modelSquaredNorm;
    function.asEigen() = _amplitude * model - _inputs.getData().asEigen();
}

void MultiGaussianObjective::computeDerivative(
    ndarray::Array<double const,1,1> const & parameters, 
    ndarray::Array<double const,1,1> const & function,
    ndarray::Array<double,2,-2> const & derivative
) {
    derivative.asEigen().setZero();
    for (std::size_t n = 0; n < _builders.size(); ++n) {
        _builders[n].computeDerivative(derivative, true);
    }
    if (!_inputs.getWeights().isEmpty()) {
        derivative.asEigen<Eigen::ArrayXpr>() 
            *= (_inputs.getWeights().asEigen() * Eigen::RowVectorXd::Ones(parameters.getSize<0>())).array();
    }
    // Right now, 'derivative' is the partial derivative w.r.t. the objective parameters
    // with flux held fixed at 1.  However, the parameters also affect the flux, so we need
    // to compute the partial derivative of that.
    Eigen::VectorXd tmp = _inputs.getData().asEigen() - 2.0 * _amplitude * _model.asEigen();
    Eigen::VectorXd dAmplitude = (derivative.asEigen().adjoint() * tmp) / _modelSquaredNorm;
    // Now we update 'derivative' so it becomes the complete derivative rather than the partial.
    derivative.asEigen() *= _amplitude;
    derivative.asEigen() += _model.asEigen() * dAmplitude.transpose();
}

}}}} // namespace lsst::meas::extensions::multiShapelet
