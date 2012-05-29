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
#ifndef MULTISHAPELET_ModelInputHandler_h_INCLUDED
#define MULTISHAPELET_ModelInputHandler_h_INCLUDED

#include "lsst/afw/geom.h"
#include "lsst/afw/geom/ellipses.h"
#include "lsst/afw/image/Image.h"
#include "lsst/afw/image/MaskedImage.h"
#include "lsst/afw/detection/Footprint.h"

namespace lsst { namespace meas { namespace extensions { namespace multiShapelet {

class ModelInputHandler {
public:

    /// @brief Center-subtracted X coordinates for the model's pixels.
    ndarray::Array<double const,1,1> getX() const { return _x; }

    /// @brief Center-subtracted Y coordinates for the model's pixels.
    ndarray::Array<double const,1,1> getY() const { return _y; }

    /// @brief Weighted data array.
    ndarray::Array<double const,1,1> getData() const  { return _data; }

    /// @brief Per-pixel weights; will be empty if initialized with an Image rather than a MaskedImage.
    ndarray::Array<double const,1,1> getWeights() const  { return _weights; }

    /// @brief The number of pixels in the model.
    int getSize() const { return _x.size(); }
    
    /// @brief Return the footprint actually used to flatten the inputs.
    PTR(afw::detection::Footprint) getFootprint() const { return _footprint; }

    template <typename PixelT>
    ModelInputHandler(afw::image::Image<PixelT> const & image, afw::geom::Point2D const & center, 
                      afw::geom::Box2I const & region);

    template <typename PixelT>
    ModelInputHandler(afw::image::Image<PixelT> const & image, afw::geom::Point2D const & center, 
                      afw::detection::Footprint const & region, int growFootprint=0);

    template <typename PixelT>
    ModelInputHandler(afw::image::Image<PixelT> const & image, afw::geom::ellipses::Ellipse const & ellipse,
                      afw::detection::Footprint const & region, int growFootprint=0);
    
    template <typename PixelT>
    ModelInputHandler(
        afw::image::MaskedImage<PixelT> const & image, afw::geom::Point2D const & center, 
        afw::geom::Box2I const & region, afw::image::MaskPixel badPixelMask=0x0, bool usePixelWeights=false
    );

    template <typename PixelT>
    ModelInputHandler(
        afw::image::MaskedImage<PixelT> const & image, afw::geom::Point2D const & center, 
        afw::detection::Footprint const & region, int growFootprint=0,
        afw::image::MaskPixel badPixelMask=0x0, bool usePixelWeights=false
    );

    template <typename PixelT>
    ModelInputHandler(
        afw::image::MaskedImage<PixelT> const & image, afw::geom::ellipses::Ellipse const & ellipse, 
        afw::detection::Footprint const & region, int growFootprint=0,
        afw::image::MaskPixel badPixelMask=0x0, bool usePixelWeights=false
    );

private:
    ndarray::Array<double,1,1> _x;
    ndarray::Array<double,1,1> _y;
    ndarray::Array<double,1,1> _data;
    ndarray::Array<double,1,1> _weights;
    PTR(afw::detection::Footprint) _footprint;
};

}}}} // namespace lsst::meas::extensions::multiShapelet

#endif // !MULTISHAPELET_ModelInputHandler_h_INCLUDED
