#!/usr/bin/env python

import numpy
import scipy.special
from matplotlib import pyplot

class ExactSersic(object):

    def __init__(self, n):
        self.n = n
        self.kappa = scipy.special.gammaincinv(2.0 * n, 0.5)

    def __call__(self, r):
        return numpy.exp(-self.kappa * (r**(1.0 / self.n) - 1.0))

exactDeVaucouleur = ExactSersic(4.0)
exactExponential = ExactSersic(1.0)

def sdssDeVaucouleur(r):
    """Truncated/softened de Vaucouleur - copied from SDSS"""
    DEFAC = -7.66925
    DEVOUT = 8.0
    DEVCUT = 7.0
    p = numpy.exp(DEFAC * ((r**2 + 0.0004)**0.125 - 1.0))
    big = r > DEVCUT
    scr = (r[big] - DEVCUT) / (DEVOUT - DEVCUT)
    scr = 1.0 - scr**2
    p[big] *= scr*scr
    p[r > DEVOUT] = 0.0
    return p

def sdssExponential(r):
    """Truncated/softened exponential - copied from SDSS"""
    EXPFAC = -1.67835
    EXPOUT = 4.0
    EXPCUT = 3.0
    p = numpy.exp(EXPFAC * (r - 1.0))
    big = r > EXPCUT
    scr = (r[big] - EXPCUT) / (EXPOUT - EXPCUT);
    scr = 1.0 - scr**2
    p[big] *= scr * scr
    p[r > EXPOUT] = 0.0
    return p

# Special values below from Dustin Lang and David Hogg.

class TractorMultiGaussian(object):

    def __init__(self, amplitude, variance):
        self.sigma = variance**0.5
        # Tractor uses integrated amplitude of each component, but I use peak value here.
        self.amplitude = amplitude / (2.0 * numpy.pi * self.sigma**2)

    def __call__(self, r):
        return (numpy.exp(-0.5 * (numpy.divide.outer(r, self.sigma))**2) * self.amplitude).sum(axis=1)

tractorExponential = TractorMultiGaussian(
    amplitude = numpy.array([3.31636565e-05, 1.06478564e-03, 1.33260624e-02, 1.06217866e-01,
                             6.09924868e-01, 2.43600369e+00, 5.34553250e+00, 3.41379672e+00]),
    variance = numpy.array([4.91509189e-05, 7.91283025e-04, 5.06909854e-03, 2.30018199e-02,
                            8.50523831e-02, 2.73398885e-01, 7.93675135e-01, 2.17065603e+00]),
)
#tractorExponential.sigma /= EXP_HALF_LIGHT_RADIUS # want radius parameter to be half-light radius

tractorDeVaucouleur = TractorMultiGaussian(
    amplitude = numpy.array([1.36305372e-02, 1.08889599e-01, 3.68235229e-01, 9.26647361e-01,
                             2.00002437e+00, 3.77384257e+00, 6.01053703e+00, 7.22968202e+00]),
    variance = numpy.array([1.34654444e-04, 5.04128747e-04, 1.88518562e-03, 7.20439754e-03,
                            2.87959626e-02, 1.25118904e-01, 6.38235086e-01, 4.76437813e+00]),
)

def integral(x, profile):
    dx = x[1] - x[0]
    z = profile(x) * x * 2.0 * numpy.pi * dx
    return numpy.cumsum(z)

def plotSingle(x, profile, color, label, marker, stride):
    y = profile(x)
    #z = integral(x, profile)
    pyplot.plot(x, y, color, alpha=0.75)
    pyplot.plot(x[::stride], y[::stride], color + marker, alpha=0.5, markeredgewidth=0, label=label)
    #pyplot.plot(x, z, color, alpha=0.75)
    #pyplot.plot(x[::stride], z[::stride], color + marker, alpha=0.5, markeredgewidth=0)

def plotFull(exact, sdss, tractor, yMax):
    x1 = numpy.linspace(0.0, 10.0, 1001)
    x2 = numpy.linspace(0.0, 2.0, 1001)
    pyplot.figure()
    pyplot.subplot(2, 1, 1)
    plotSingle(x1, exact, "r", "exact", "o", 50)
    plotSingle(x1, sdss, "g", "sdss", "^", 50)
    plotSingle(x1, tractor, "b", "tractor", "s", 50)
    pyplot.ylim(1E-5, 1E2)
    pyplot.semilogy()
    pyplot.subplot(2, 1, 2)
    plotSingle(x2, exact, "r", "exact", "o", 125)
    plotSingle(x2, sdss, "g", "sdss", "^", 125)
    plotSingle(x2, tractor, "b", "tractor", "s", 125)
    pyplot.ylim(0, yMax)
    pyplot.legend()

def main():
    plotFull(exactExponential, sdssExponential, tractorExponential, yMax=12)
    plotFull(exactDeVaucouleur, sdssDeVaucouleur, tractorDeVaucouleur, yMax=100)
    pyplot.show()

if __name__ == "__main__":
    main()
