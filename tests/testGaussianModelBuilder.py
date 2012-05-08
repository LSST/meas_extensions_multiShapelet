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
Tests for GaussianModelBuilder

Run with:
   ./testGaussianModelBuilder.py
or
   python
   >>> import testGaussianModelBuilder; testGaussianModelBuilder.run()
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
numpy.set_printoptions(linewidth=110)

#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

class GaussianModelBuilderTestCase(unittest.TestCase):

    def assertClose(self, a, b, rtol=1E-5, atol=1E-8):
        self.assert_(numpy.allclose(a, b, rtol=rtol, atol=atol), "\n%s\n!=\n%s" % (a, b))

    def buildModel(self, ellipse):
        LT = geom.LinearTransform
        gt = ellipse.getGridTransform()
        xt = gt[LT.XX] * self.xg + gt[LT.XY] * self.yg
        yt = gt[LT.YX] * self.xg + gt[LT.YY] * self.yg
        model = numpy.exp(-0.5 * (yt**2 + xt**2))
        return model.ravel()

    def buildNumericalDerivative(self, builder, parameters, makeEllipse):
        eps = 1E-6
        derivative = numpy.zeros((len(parameters), builder.getSize()), dtype=float).transpose()
        for i in range(len(parameters)):
            parameters[i] += eps
            ellipse = makeEllipse(parameters)
            builder.update(ellipse)
            derivative[:,i] = builder.getModel()
            parameters[i] -= 2.0 * eps
            ellipse = makeEllipse(parameters)
            builder.update(ellipse)
            derivative[:,i] -= builder.getModel()
            derivative[:,i] /= 2.0 * eps
        return derivative

    def setUp(self):
        self.ellipse = ellipses.Axes(10, 7, 0.3)
        self.xg, self.yg = numpy.meshgrid(numpy.linspace(-20, 20, 101), numpy.linspace(-15, 25, 95))
        self.x = self.xg.ravel()
        self.y = self.yg.ravel()
        self.model = self.buildModel(self.ellipse)

    def tearDown(self):
        del self.ellipse

    def testModel(self):
        builder = ms.GaussianModelBuilder(self.x, self.y)
        builder.update(self.ellipse)
        self.assertClose(builder.getModel(), self.model)

    def testDerivative1(self):
        """test derivative with no reparameterization"""
        builder = ms.GaussianModelBuilder(self.x, self.y)
        a = numpy.zeros((3, builder.getSize()), dtype=float).transpose()
        jac = numpy.identity(3, dtype=float)
        builder.update(self.ellipse)
        builder.computeDerivative(a, jac)
        def makeAxesEllipse(p):
            return ellipses.Axes(*p)
        n = self.buildNumericalDerivative(builder, self.ellipse.getParameterVector(), makeAxesEllipse)
        # no hard requirement for tolerances here, but I've dialed them to the max to avoid regressions
        self.assertClose(a, n, rtol=1E-4, atol=1E-6)

    def testDerivative2(self):
        """test derivative with nontrivial reparameterization (derivative wrt different core)"""
        builder = ms.GaussianModelBuilder(self.x, self.y)
        builder.update(self.ellipse)
        quad = ellipses.Quadrupole(self.ellipse)
        jac = self.ellipse.dAssign(quad)
        a = numpy.zeros((3, builder.getSize()),dtype=float).transpose()
        builder.update(self.ellipse)
        builder.computeDerivative(a, jac)
        def makeQuadrupole(p):
            return ellipses.Quadrupole(*p)
        n = self.buildNumericalDerivative(builder, quad.getParameterVector(), makeQuadrupole)
        # no hard requirement for tolerances here, but I've dialed them to the max to avoid regressions
        self.assertClose(a, n, rtol=1E-4, atol=1E-6)


#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

def suite():
    """Returns a suite containing all the test cases in this module."""

    utilsTests.init()

    suites = []
    suites += unittest.makeSuite(GaussianModelBuilderTestCase)
    suites += unittest.makeSuite(utilsTests.MemoryTestCase)
    return unittest.TestSuite(suites)

def run(shouldExit=False):
    """Run the tests"""
    utilsTests.run(suite(), shouldExit)

if __name__ == "__main__":
    run(True)
