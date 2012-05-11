#!/usr/bin/env python

# 
# LSST Data Management System
# Copyright 2008, 2009, 2010 LSST Corporation.
# 
# This product includes software developed by the
# LSST Project (http://www.lsst.org/).
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the LSST License Statement and 
# the GNU General Public License along with this program.  If not, 
# see <http://www.lsstcorp.org/LegalNotices/>.
#

"""
Tests for FitPsf

Run with:
   ./testFitPsf.py
or
   python
   >>> import testFitPsf; testFitPsf.run()
"""

import unittest
import numpy

import lsst.utils.tests as utilsTests
import lsst.pex.exceptions
import lsst.afw.geom as geom
import lsst.afw.geom.ellipses as ellipses
import lsst.afw.image
import lsst.afw.detection
import lsst.meas.extensions.multiShapelet as ms

numpy.random.seed(5)
numpy.set_printoptions(linewidth=120)

from matplotlib import pyplot

#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

class FitProfileTestMixin(object):

    def assertClose(self, a, b, rtol=1E-5, atol=1E-8):
        self.assert_(numpy.allclose(a, b, rtol=rtol, atol=atol), "\n%s\n!=\n%s" % (a, b))

    def setUp(self):
        self.center = geom.Point2D(50.1, 60.7)
        self.ellipse = ellipses.Ellipse(ellipses.Axes(20.0, 15.0, 0.23), self.center)
        self.footprint = lsst.afw.detection.Footprint(self.ellipse)
        self.bbox = geom.Box2I(geom.Point2I(5, 6), geom.Extent2I(80, 75))
        shape = (self.bbox.getHeight(), self.bbox.getWidth())
        self.mi = lsst.afw.image.MaskedImageF(self.bbox)
        self.mi.getImage().getArray()[:,:] += numpy.random.randn(*shape)
        self.mi.getVariance().getArray()[:,:] = 1.0
        bad = lsst.afw.image.MaskU.getPlaneBitMask("BAD")
        mask = numpy.random.randn(*shape) > 2.5
        assert mask.sum() > 0  # not a real test, just guard against unlikely random numbers
        self.mi.getMask().getArray()[mask] = bad
        self.mi.getImage().getArray()[mask] = float("NaN")
        gt = self.ellipse.getGridTransform()
        x, y = numpy.meshgrid(
            numpy.arange(self.bbox.getBeginX(), self.bbox.getEndX()),
            numpy.arange(self.bbox.getBeginY(), self.bbox.getEndY())
            )
        xt = gt[geom.AffineTransform.XX] * x + gt[geom.AffineTransform.XY] * y + gt[geom.AffineTransform.X]
        yt = gt[geom.AffineTransform.YX] * x + gt[geom.AffineTransform.YY] * y + gt[geom.AffineTransform.Y]
        self.x = x.astype(float).ravel() - self.center.getX()
        self.y = y.astype(float).ravel() - self.center.getY()
        self.mi.getImage().getArray()[:,:] += 20.0 * numpy.exp(0.5 * xt**2 + yt**2)
        self.ctrl = self.config.makeControl()
        self.inputs = ms.ModelInputHandler(self.mi, self.center,
                                           self.footprint, self.ctrl.growFootprint, bad, False)

    def testModel(self):
        components = ms.MultiGaussianRegistry.lookup(self.ctrl.profile)
        self.assertClose(components.integrate(), 1.0)
        parameters = numpy.array([[0.0, 0.0, 1.0], [0.0, 0.0, -16], [0.2,-0.8,2.3]])
        sImage = lsst.afw.image.ImageD(self.bbox)
        builders = [
            ms.GaussianModelBuilder(self.x, self.y, c.flux, c.radius)
            for c in components
            ]
        for row in parameters:
            model = ms.FitProfileModel(self.ctrl, 1.0, row)
            msf = model.asMultiShapelet(self.center)
            self.assertClose(msf.evaluate().integrate(), 1.0)
            z0 = sImage.getArray()
            z0[:,:] = 0.0
            msf.evaluate().addToImage(sImage)
            z1 = numpy.zeros(z0.shape, dtype=float)
            for b in builders:
                b.update(ms.MultiGaussianObjective.readParameters(row))
                z1 += b.getModel().reshape(*z1.shape)
            self.assertClose(z0, z1)

    def testConvolvedModel(self):
        psfModel = ms.FitPsfModel(ms.FitPsfControl(), 1.0, numpy.array([0.1, -0.05, 1.0]))
        psfComponents = psfModel.getComponents()
        components = ms.MultiGaussianRegistry.lookup(self.ctrl.profile)
        parameters = numpy.array([[0.0, 0.0, 0.0], [0.0, 0.0, -16], [0.2,-0.8,2.3]])
        sImage = lsst.afw.image.ImageD(self.bbox)
        builders = []
        for p in psfComponents:
            psfEllipse = psfModel.ellipse.clone()
            psfEllipse.scale(p.radius)
            builders.extend(
                ms.GaussianModelBuilder(self.x, self.y, c.flux, c.radius, psfEllipse, p.flux)
                for c in components
                )
        psfShapelets = psfModel.asMultiShapelet()
        for row in parameters:
            model = ms.FitProfileModel(self.ctrl, 1.0, row)
            msf = model.asMultiShapelet(self.center).convolve(psfShapelets)
            z0 = sImage.getArray()
            z0[:,:] = 0.0
            msf.evaluate().addToImage(sImage)
            z1 = numpy.zeros(z0.shape, dtype=float)
            for b in builders:
                b.update(ms.MultiGaussianObjective.readParameters(row))
                z1 += b.getModel().reshape(*z1.shape)
            self.assertClose(z0, z1)

    def tearDown(self):
        del self.ellipse
        del self.footprint
        del self.bbox
        del self.mi
        del self.config
        del self.ctrl
        del self.inputs

class FitExponentialTestCase(unittest.TestCase,FitProfileTestMixin):

    def setUp(self):
        self.config = ms.FitExponentialConfig()
        FitProfileTestMixin.setUp(self)

    def tearDown(self):
        FitProfileTestMixin.tearDown(self)

#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

def suite():
    """Returns a suite containing all the test cases in this module."""

    utilsTests.init()

    suites = []
    suites += unittest.makeSuite(FitExponentialTestCase)
    suites += unittest.makeSuite(utilsTests.MemoryTestCase)
    return unittest.TestSuite(suites)

def run(shouldExit=False):
    """Run the tests"""
    utilsTests.run(suite(), shouldExit)

if __name__ == "__main__":
    run(True)
