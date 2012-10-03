import numpy
import scipy.optimize
import scipy.special
from matplotlib import pyplot

class ExactSersic(object):

    def __init__(self, n):
        self.n = n
        self.kappa = scipy.special.gammaincinv(2.0 * n, 0.5)

    def __call__(self, r):
        return numpy.exp(-self.kappa * (r**(1.0 / self.n) - 1.0))

# Special values below from Dustin Lang and David Hogg.

class TractorMultiGaussian(object):

    def __init__(self, flux, variance):
        self.sigma = variance**0.5
        self.flux = flux

    def __call__(self, r):
        return self.decompose(r).sum(axis=1)

    def decompose(self, r):
        return numpy.exp(-0.5 * (numpy.divide.outer(r, self.sigma))**2) * self.flux \
            / (2.0 * numpy.pi * self.sigma**2)

tractorExponential = TractorMultiGaussian(
    flux = numpy.array([3.31636565e-05, 1.06478564e-03, 1.33260624e-02, 1.06217866e-01,
                        6.09924868e-01, 2.43600369e+00, 5.34553250e+00, 3.41379672e+00]),
    variance = numpy.array([4.91509189e-05, 7.91283025e-04, 5.06909854e-03, 2.30018199e-02,
                            8.50523831e-02, 2.73398885e-01, 7.93675135e-01, 2.17065603e+00]),
)

tractorDeVaucouleur = TractorMultiGaussian(
    flux = numpy.array([1.36305372e-02, 1.08889599e-01, 3.68235229e-01, 9.26647361e-01,
                        2.00002437e+00, 3.77384257e+00, 6.01053703e+00, 7.22968202e+00]),
    variance = numpy.array([1.34654444e-04, 5.04128747e-04, 1.88518562e-03, 7.20439754e-03,
                            2.87959626e-02, 1.25118904e-01, 6.38235086e-01, 4.76437813e+00]),
)

def safe_nnls(a, b):
    try:
        return scipy.optimize.nnls(a, b)
    except RuntimeError:
        pass
    def func(x):
        d = numpy.dot(a, x)
        d -= b
        f = numpy.dot(d, d)
        g = numpy.dot(a.transpose(), d)
        return f, g
    x0, residuals, rank, sv = numpy.linalg.lstsq(a, b)
    x0[x0 < 0] = 0.0
    bounds = [(0, float("inf"))] * a.shape[1]
    x, neval, rc = scipy.optimize.fmin_tnc(func, x0, bounds=bounds, disp=0)
    f, g = func(x)
    return x, f

class SersicApproximator(object):

    sigmaMin = 1e-2
    sigmaMax = 10.0
    rMin = 1E-2
    rMax = 10.0
    rPoints = 500
    nMin = 1.0
    nMax = 4.0
    nPoints = 5
    fitTolerance = 1E-8

    def __init__(self, weights):
        self.r = numpy.linspace(self.rMin, self.rMax, self.rPoints)
        self.n = numpy.linspace(self.nMin, self.nMax, self.nPoints)
        self.rGrid, self.nGrid = numpy.meshgrid(self.r, self.n)
        self.sersicGrid = self.buildSersicGrid()
        self.z = -0.5 * self.r**2
        self.setWeights(weights)
        self.beta = numpy.ones(self.n.shape, dtype=float)
        self.alpha = None

    def setWeights(self, weights):
        self.weightGrid = weights(self.rGrid, self.nGrid)
        self.wSersicGrid = self.weightGrid * self.sersicGrid

    def buildSersicGrid(self):
        grid = numpy.zeros(self.n.shape + self.r.shape, dtype=float)
        for i, n in enumerate(self.n):
            f = ExactSersic(n)
            grid[i,:] = f(self.r)
        return grid

    def buildGaussianMatrix(self, sigma):
        return numpy.exp(numpy.multiply.outer(sigma**(-2), self.z))

    def fit(self, sigma, i, fitAlpha=False):
        n = self.n[i]
        def func(param):
            matrix = self.buildGaussianMatrix(sigma * param)
            matrix *= self.weightGrid[i,:]
            # NNLS: linear least squares, requiring solution >= 0
            alpha, residual = safe_nnls(matrix.transpose(), self.wSersicGrid[i,:])
            return residual
        if numpy.abs(self.beta[i] - 1.0) < 0.5:
            p1 = 0.5
            p2 = 2.0
        if self.beta[i] > 1.0:
            p1 = 0.5
            p2 = self.beta[i]
        else:
            p1 = self.beta[i]
            p2 = 2.0
        f1 = func(p1)
        f2 = func(p2)
        if f1 < f2:
            bracket = (0.0, p1, p2)
        else:
            bracket = (p1, p2)
        beta, residual, _, _ = scipy.optimize.brent(func, brack=bracket, tol=self.fitTolerance,
                                                    full_output=True)
        if fitAlpha:
            matrix = self.buildGaussianMatrix(sigma * beta)
            matrix *= self.weightGrid[i,:]
            alpha, _ = safe_nnls(matrix.transpose(), self.wSersicGrid[i,:])
        else:
            alpha = None
        return residual, beta, alpha

    def evaluate(self, sigma, fitAlpha=False, plot=False):
        self.alpha = numpy.zeros(self.n.shape + sigma.shape, dtype=float) if (fitAlpha or plot) else None
        self.beta[:] = 1.0
        fullResidual = 0.0
        for i, n in enumerate(self.n):
            residual, beta, alpha = self.fit(sigma, i, fitAlpha=(fitAlpha or plot))
            self.beta[i] = beta
            if fitAlpha or plot:
                self.alpha[i,:] = alpha
            if plot:
                def plotComparison(n, xlim, ylim):
                    pyplot.subplot(3,2,n)
                    pyplot.semilogy(self.r, self.sersicGrid[i,:], 'k', alpha=0.5)
                    pyplot.semilogy(self.r, approx, 'r')
                    pyplot.xlim(*xlim)
                    pyplot.ylim(*ylim)
                def plotDifference(n, xlim, ylim):
                    pyplot.subplot(3,2,n)
                    pyplot.axhline(0.0, color='k', alpha=0.5)
                    pyplot.plot(self.r, self.sersicGrid[i,:] - approx, 'r')
                    pyplot.xlim(*xlim)
                    pyplot.ylim(*ylim)
                def plotRelDifference(n, xlim, ylim):
                    pyplot.subplot(3,2,n)
                    pyplot.axhline(0.0, color='k', alpha=0.5)
                    pyplot.plot(self.r, (self.sersicGrid[i,:] - approx) / self.sersicGrid[i,:], 'r')
                    pyplot.xlim(*xlim)
                    pyplot.ylim(*ylim)
                pyplot.clf()
                pyplot.suptitle("n = %f (%d of %d)" % (n, i+1, self.n.size))
                matrix = self.buildGaussianMatrix(sigma * beta)
                approx = numpy.dot(matrix.transpose(), alpha)
                plotComparison(1, xlim=(self.rMin, self.rMax), ylim=(1E-8, numpy.max(self.sersicGrid)))
                plotComparison(2, xlim=(0,1), ylim=(1E-1, numpy.max(self.sersicGrid)))
                plotDifference(3, xlim=(self.rMin, self.rMax), ylim=(-1, 1))
                plotDifference(4, xlim=(0,1), ylim=(-1, 1))
                plotRelDifference(5, xlim=(self.rMin, self.rMax), ylim=(-1, 1))
                plotRelDifference(6, xlim=(0,1), ylim=(-1, 1))
                print "beta", beta
                print "alpha", alpha
                raw_input("Press Enter to continue...")
            fullResidual += residual
        return fullResidual

    def optimize(self, sigma):
        def iecons(p):
            result = numpy.zeros(sigma.size, dtype=float)
            result[0] = p[0]
            result[1:] = p[1:] - p[:-1]
            return result
        def d_iecons(p):
            result = numpy.zeros((sigma.size, sigma.size), dtype=float)
            result[0,0] = 1.0
            for i in range(1, sigma.size):
                result[i,i] = 1.0
                result[i,i-1] = -1.0
            return result
        sigma, f, niter, imode, smode = scipy.optimize.fmin_slsqp(self.evaluate, sigma,
                                                                  f_ieqcons=iecons, fprime_ieqcons=d_iecons,
                                                                  full_output=True)
        sigma = numpy.array(sigma, dtype=float)
        betaMid = numpy.median(self.beta)
        self.beta /= betaMid
        sigma *= betaMid
        return sigma

class weights(object):

    @staticmethod
    def uniform1d(r, n):
        return 1.0

    @staticmethod
    def uniform2d(r, n):
        return r**0.5

    class Smooth(object):
        def __init__(self, radius=2.0, power=0.4):
            self.radius = radius
            self.power = power
        def __call__(self, r, n):
            return (numpy.exp(-(r**2)/(self.radius**2))*r)**(self.power*0.5)
