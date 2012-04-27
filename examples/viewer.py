import lsst.meas.extensions.multiShapelet as ms
import lsst.afw.image
import numpy
from matplotlib import pyplot

class FitPsfViewer(object):

    class Iteration(object):

        def __init__(self, opt, viewer):
            obj = opt.getObjective().cast(ms.MultiGaussianObjective)
            self.state = opt.getState()
            if self.state & ms.HybridOptimizer.STEP_ACCEPTED:
                self.parameters = opt.getParameters().copy()
                self.chisq = opt.getChiSq()
            else:
                self.parameters = opt.getTrialParameters().copy()
                self.chisq = opt.getTrialChiSq()
            self.fNorm = opt.getFunctionInfNorm()
            self.gNorm = opt.getGradientInfNorm()
            self.amplitude = obj.getAmplitude()
            self.method = opt.getMethod()
            self.model = ms.FitPsfModel(viewer.ctrl, self.amplitude, self.parameters)

        def __str__(self):
            return "state=0x%1.2x, method=%s, chisq=%f, fNorm=%g, gNorm=%g, amplitude=%s, parameters=%s" % (
                self.state, "LM" if self.method==ms.HybridOptimizer.LM else "BFGS",
                self.chisq, self.fNorm, self.gNorm, self.amplitude, self.parameters
                )

        __repr__ = __str__

    def __init__(self, source, psf, ctrl=None):
        if ctrl is None:
            ctrl = ms.FitPsfControl()
        self.ctrl = ctrl
        self.saved = ms.FitPsfModel(self.ctrl, source)
        self.center = source.getCentroid()
        self.image = psf.computeImage(self.center)
        opt = ms.FitPsfAlgorithm.makeOptimizer(self.ctrl, self.image, self.center)
        maxIter = opt.getControl().maxIter
        self.iterations = [self.Iteration(opt, self)]
        for self.iterCount in range(maxIter):
            opt.step()
            self.iterations.append(self.Iteration(opt, self))
            if opt.getState() & ms.HybridOptimizer.FINISHED:
                break
        self.model = ms.FitPsfModel(self.iterations[-1].model)
        ms.FitPsfAlgorithm.fitShapeletTerms(self.ctrl, self.image, self.center, self.model)

    @staticmethod
    def _plotImage(image, title=None, ellipses=()):
        bbox = image.getBBox(lsst.afw.image.PARENT)
        array = image.getArray()
        pyplot.imshow(array, interpolation='nearest', origin='lower',
                      extent=(bbox.getMinX()-0.5, bbox.getMaxX()+0.5, bbox.getMinY()-0.5, bbox.getMaxY()+0.5)
                      )
        if title is not None:
            pyplot.title(title)
        for ellipse in ellipses:
            ellipse.plot(fill=False, rescale=False)
        ticks = [array.min(), 0.5 * (array.min() + array.max()), array.max()]
        pyplot.colorbar(orientation="horizontal", format="%.2g", ticks=ticks)

    def plot(self, model=None):
        if model is None:
            model = self.model
        elif model == "saved":
            model = self.saved
        elif not isinstance(model, ms.FitPsfModel):
            model = self.iterations[model].model
        image = lsst.afw.image.ImageD(self.image.getBBox(lsst.afw.image.PARENT))
        func = model.asMultiShapelet(self.center)
        func.evaluate().addToImage(image)
        residuals = lsst.afw.image.ImageD(self.image, True)
        residuals -= image
        outerEllipse = lsst.afw.geom.ellipses.Quadrupole(model.ellipse)
        outerEllipse.scale(model.radiusRatio)
        ellipses = [
            lsst.afw.geom.ellipses.Ellipse(model.ellipse, self.center),
            lsst.afw.geom.ellipses.Ellipse(outerEllipse, self.center),
            ]
        pyplot.clf()
        pyplot.subplot(1,3,1)
        self._plotImage(self.image, title="original PSF model image", ellipses=ellipses)
        pyplot.subplot(1,3,2)
        self._plotImage(image, title="shapelet fit", ellipses=ellipses)
        pyplot.subplot(1,3,3)
        self._plotImage(residuals, title="original-shapelets", ellipses=ellipses)

