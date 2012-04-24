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
#ifndef MULTISHAPELET_EllipseSquaredNorm_h_INCLUDED
#define MULTISHAPELET_EllipseSquaredNorm_h_INCLUDED

#include "lsst/afw/geom/ellipses.h"

namespace lsst { namespace meas { namespace extensions { namespace multiShapelet {

/**
 *  @brief A low-level functor to aid in evaluating elliptically-symmetric functions.
 *
 *  The "ellipse squared norm" this class computes is
 *  @f[
 *     z = \left[\begin{array}{ c c }
 *            x & y
 *         \end{array}\right]
 *         \left[\begin{array}{ c c }
 *            Q_{xx} & Q_{xy} \\
 *            Q_{yx} & Q_{yy}
 *         \end{array}\right]
 *         \left[\begin{array}{ c }
 *           x \\
 *           y
 *         \end{array}\right]
 *  @f]
 *  where @f$Q@f$ is the quadrupole matrix definition of an ellipse.
 */
class EllipseSquaredNorm {
public:

    //@{
    /**
     *  @brief Compute the squared ellipse norm at x and y and return it in z.
     *
     *  x,y ---- Input coordinates.
     *  rx,ry -- Intermediate results that can be passed to the derivative functions.
     *  z ------ Output squared ellipse norm (see class docs).
     */
    void operator()(
        double const x, double const y,
        double & rx, double & ry,
        double & z
    ) const {
        rx = x * _r11;
        ry = x * _r12 + y * _r22;
        z = rx * rx + ry * ry;
    }
    template <typename D1, typename D2, typename D3>
    void operator()(
        Eigen::MatrixBase<D1> const & x, Eigen::MatrixBase<D1> const & y,
        Eigen::MatrixBase<D2> & rx, Eigen::MatrixBase<D2> & ry,
        Eigen::MatrixBase<D3> & z
    ) const {
        rx = x * _r11;
        ry = x * _r12 + y * _r22;
        z.array() = rx.array().square() + ry.array().square();
    }
    //@}

    //@{
    /**
     *  @brief Compute the derivative of the squared ellipse norm wrt the input coordinates.
     *
     *  x,y ------------- Input coordinates (actually ignored; present for consistency).
     *  rx,ry ----------- Intermediate results from operator().
     *  dz_dx, dz_dy ---- Output derivative of z wrt x and y.
     */
    void dCoords(
        double const x, double const y,
        double const rx, double const ry,
        double & dz_dx, double & dz_dy
    ) const {
        dz_dx = 2.0 * (_r11 * rx + _r12 * ry);
        dz_dy = 2.0 * _r22 * ry;
    }
    template <typename D1, typename D2, typename D3>
    void dCoords(
        Eigen::MatrixBase<D1> const & x, Eigen::MatrixBase<D1> const & y,
        Eigen::MatrixBase<D2> const & rx, Eigen::MatrixBase<D2> const & ry,
        Eigen::MatrixBase<D3> & dz_dx, Eigen::MatrixBase<D3> & dz_dy
    ) const {
        dz_dx = (2.0 * _r11) * rx + (2.0 * _r12) * ry;
        dz_dy = (2.0 * _r22) * ry;
    }
    //@}

    //@{
    /**
     *  @brief Compute the derivative of the squared ellipse norm wrt the ellipse parameters.
     *
     *  x,y ---------- Input coordinates.
     *  rx,ry -------- Intermediate results from operator() (destroyed on output).
     *  jacobian ----- Derivative of internal representation wrt ellipse parameters; see update().
     *  dz ----------- Output derivative of z wrt ellipse parameters (should be initialized to zero).
     *
     *  In the vectorized form, the x, y, rx, and ry arguments must be column vectors.
     *  The number of columns of the Jacobian matrix should equal the number of columns of dz.
     *  The Jacobian matrix must have exactly 3 rows.  The output dz should be a row vector
     *  in the scalar case and a matrix with number of rows equal to the number of coordinates
     *  in the vectorized case.
     */
    template <typename D1, typename D2>
    void dEllipse(
        double const x, double const y,          // scalar
        double & rx, double & ry,                // scalar
        Eigen::MatrixBase<D1> const & jacobian,  // 3xM elements
        Eigen::MatrixBase<D2> & dz               // 1xM elements (row vector),
    ) const {
        rx *= 2.0 * x;
        dz += rx * jacobian.row(0);
        rx = 2.0 * ry * x;  // use rx as temporary for 2*ry*x 
        dz += rx * jacobian.row(1);
        ry *= 2.0 * y;
        dz += ry * jacobian.row(2);
    }
    template <typename D1, typename D2, typename D3, typename D4>
    void dEllipse(
        Eigen::MatrixBase<D1> const & x, Eigen::MatrixBase<D1> const & y, // Nx1 elements (col vector)
        Eigen::MatrixBase<D2> & rx, Eigen::MatrixBase<D2> & ry,           // Nx1 elements (col vector)
        Eigen::MatrixBase<D3> const & jacobian,                         // 3xM elements
        Eigen::MatrixBase<D4> & dz                                      // NxM elements
    ) const {
        rx.array() *= 2.0 * x.array();
        dz += rx * jacobian.row(0);
        rx.array() = 2.0 * ry.array() * x.array();  // use rx as temporary for 2*ry*x 
        dz += rx * jacobian.row(1);
        ry.array() *= 2.0 * y.array();
        dz += ry * jacobian.row(2);
    }
    //@}

    /**
     *  @brief Update the ellipse and optionally return the derivative of the internal representation
     *         wrt the ellipse parameters.
     *
     *  The returned matrix can be passed directly as the 'jacobian' argument of dEllipse() to compute
     *  the derivative wrt the ellipse parameters, or it can be multiplied on the right by another matrix
     *  to compute the derivative wrt some function of the ellipse parameters.
     */
    afw::geom::ellipses::BaseCore::Jacobian update(
        afw::geom::ellipses::BaseCore const & ellipse,
        bool computeJacobian=true
    );

    // @brief Construct the squared norm functor with a unit circle.
    explicit EllipseSquaredNorm() : _r11(1.0), _r12(0.0), _r22(1.0) {}

private:
    double _r11;
    double _r12;
    double _r22;
};


}}}} // namespace lsst::meas::extensions::multiShapelet

#endif // !MULTISHAPELET_EllipseSquaredNorm_h_INCLUDED
