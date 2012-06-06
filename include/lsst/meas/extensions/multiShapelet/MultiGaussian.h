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
#ifndef MULTISHAPELET_MultiGaussian_h_INCLUDED
#define MULTISHAPELET_MultiGaussian_h_INCLUDED

#include <vector>

#include "lsst/afw/geom/Angle.h"
#include "lsst/shapelet/ShapeletFunction.h"

namespace lsst { namespace meas { namespace extensions { namespace multiShapelet {

struct GaussianComponent {
    double flux;
    double radius;

    shapelet::ShapeletFunction makeShapelet(afw::geom::ellipses::Ellipse const & ellipse, int order=0) const;

    void readShapeletAmplitude(double coeff0, afw::geom::ellipses::BaseCore const & ellipse);

    explicit GaussianComponent(double flux_=1.0, double radius_=1.0) 
        : flux(flux_), radius(radius_) {}
};

class MultiGaussian {
    typedef std::vector<GaussianComponent> List;
public:

    typedef List::iterator iterator;
    typedef List::const_iterator const_iterator;

    double integrate() const;

    afw::geom::ellipses::Quadrupole deconvolve(
        afw::geom::ellipses::Quadrupole const & fullMoments,
        afw::geom::ellipses::Quadrupole const & psfMoments,
        MultiGaussian const & psfMultiGaussian
    ) const;

    iterator begin() { return _components.begin(); }
    const_iterator begin() const { return _components.begin(); }

    iterator end() { return _components.end(); }
    const_iterator end() const { return _components.end(); }

    std::size_t size() const { return _components.size(); }

    GaussianComponent & add(GaussianComponent const & component) {
        _components.push_back(component);
       return _components.back();
    }

    GaussianComponent & operator[](std::size_t n) { return _components[n]; }
    GaussianComponent const & operator[](std::size_t n) const { return _components[n]; }

private:
    List _components;
};



}}}} // namespace lsst::meas::extensions::multiShapelet

#endif // !MULTISHAPELET_MultiGaussian_h_INCLUDED
