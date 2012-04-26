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
#ifndef MULTISHAPELET_MultiGaussianRegistry_h_INCLUDED
#define MULTISHAPELET_MultiGaussianRegistry_h_INCLUDED

#include "ndarray.h"

#include "lsst/base.h"
#include "lsst/meas/extensions/multiShapelet/MultiGaussian.h"

namespace lsst { namespace meas { namespace extensions { namespace multiShapelet {

/**
 *  @brief Namespace/hidden singleton class for registering and retrieving multi-Gaussian
 *         profiles by name.
 *
 *  Lookups are linear in the number of elements, but the most-recently used item is always
 *  checked first.
 *
 *  This is populated at import time in the Python module.
 */
class MultiGaussianRegistry {
public:

#ifndef SWIG

    /// @brief Retrieve the MultiGaussianList with the given name or throw NotFoundException.
    static MultiGaussianList const & lookup(std::string const & name);

    /// @brief Insert a new MultiGaussianList (replaces if name is already present).
    static void insert(std::string const & name, MultiGaussianList const & components);

#endif // !SWIG

    /// @brief Insert a new MultiGaussianList by passing in separate amplitude and radius arrays.
    static void insert(
        std::string const & name,
        ndarray::Array<double const,1> const & amplitudes,
        ndarray::Array<double const,1> const & radii
    );

};

}}}} // namespace lsst::meas::extensions::multiShapelet

#endif // !MULTISHAPELET_MultiGaussianRegistry_h_INCLUDED
