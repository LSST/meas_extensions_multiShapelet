#!/usr/bin/env python

import numpy
from matplotlib import pyplot

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
        # Tractor uses integrated amplitude of each component, we use peak value.
        self.amplitude = amplitude / (2.0 * numpy.pi * self.sigma**2)

    def __call__(self, r):
        return (numpy.exp(-0.5 * (numpy.divide.outer(r, self.sigma))) * self.amplitude).sum(axis=1)

tractorExponential = TractorMultiGaussian(
    amplitude = numpy.array([3.31636565e-05, 1.06478564e-03, 1.33260624e-02, 1.06217866e-01,
                             6.09924868e-01, 2.43600369e+00, 5.34553250e+00, 3.41379672e+00]),
    variance = numpy.array([4.91509189e-05, 7.91283025e-04, 5.06909854e-03, 2.30018199e-02,
                            8.50523831e-02, 2.73398885e-01, 7.93675135e-01, 2.17065603e+00]),
)

tractorDeVaucouleur = TractorMultiGaussian(
    amplitude = numpy.array([1.36305372e-02, 1.08889599e-01, 3.68235229e-01, 9.26647361e-01,
                             2.00002437e+00, 3.77384257e+00, 6.01053703e+00, 7.22968202e+00]),
    variance = numpy.array([1.34654444e-04, 5.04128747e-04, 1.88518562e-03, 7.20439754e-03,
                            2.87959626e-02, 1.25118904e-01, 6.38235086e-01, 4.76437813e+00]),
)

