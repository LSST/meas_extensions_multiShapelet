import lsst.meas.extensions.multiShapelet as ms
import lsst.afw.image
import numpy
from matplotlib import pyplot

class ViewerBase(object):

    class Iteration(object):

        def __init__(self, opt, viewer):
            obj = opt.getObjective().cast(ms.MultiGaussianObjective)
            self.state = opt.getState()
            if self.state & ms.HybridOptimizer.STEP_ACCEPTED:
                self.parameters = opt.getParameters().copy()
                self.func = opt.getFunction().copy()
                self.chisq = opt.getChiSq() / (self.func.size - self.parameters.size)
            else:
                self.parameters = opt.getTrialParameters().copy()
                self.func = opt.getTrialFunction().copy()
                self.chisq = opt.getTrialChiSq() / (self.func.size - self.parameters.size)
            self.fNorm = opt.getFunctionInfNorm()
            self.gNorm = opt.getGradientInfNorm()
            self.amplitude = obj.getAmplitude()
            self.method = opt.getMethod()
            self.model = viewer.Model(viewer.ctrl, self.amplitude, self.parameters)

        def __str__(self):
            return "state=0x%1.2x, method=%s, chisq=%f, fNorm=%g, gNorm=%g, amplitude=%s, parameters=%s" % (
                self.state, "LM" if self.method==ms.HybridOptimizer.LM else "BFGS",
                self.chisq, self.fNorm, self.gNorm, self.amplitude, self.parameters
                )

        __repr__ = __str__

    @staticmethod
    def _plotImage(image, title=None, ellipses=(), vmin=None, vmax=None):
        bbox = image.getBBox(lsst.afw.image.PARENT)
        array = image.getArray()
        if vmin is None or vmax is None:
            valid = array[numpy.isfinite(array)]
            if vmin is None:
                vmin = valid.min() - 0.1 * (valid.max() - valid.min())
            if vmax is None:
                vmax = valid.max() + 0.1 * (valid.max() - valid.min())
        pyplot.imshow(array, interpolation='nearest', origin='lower', vmin=vmin, vmax=vmax,
                      extent=(bbox.getMinX()-0.5, bbox.getMaxX()+0.5, bbox.getMinY()-0.5, bbox.getMaxY()+0.5)
                      )
        if title is not None:
            pyplot.title(title)
        for ellipse in ellipses:
            ellipse.plot(fill=False, rescale=False)
        return vmin, vmax

class FitPsfViewer(ViewerBase):

    Model = ms.FitPsfModel
    Control = ms.FitPsfControl
    Algorithm = ms.FitPsfAlgorithm

    def __init__(self, source, psf, ctrl=None):
        if ctrl is None:
            ctrl = self.Control()
        self.ctrl = ctrl
        self.saved = self.Model(self.ctrl, source)
        self.center = source.getCentroid()
        self.image = psf.computeImage(self.center)
        self.inputs = ms.ModelInputHandler(self.image, self.center, self.image.getBBox(lsst.afw.image.PARENT))
        opt = self.Algorithm.makeOptimizer(self.ctrl, self.inputs)
        maxIter = opt.getControl().maxIter
        self.iterations = [self.Iteration(opt, self)]
        for self.iterCount in range(maxIter):
            opt.step()
            self.iterations.append(self.Iteration(opt, self))
            if opt.getState() & ms.HybridOptimizer.FINISHED:
                break
        self.model = self.Model(self.iterations[-1].model)
        self.Algorithm.fitShapeletTerms(self.ctrl, self.inputs, self.model)


    def plot(self, model=None):
        if model is None:
            model = self.model
        elif model == "saved":
            model = self.saved
        elif not isinstance(model, self.Model):
            model = self.iterations[model].model
        data = self.image
        fit = lsst.afw.image.ImageD(data.getBBox(lsst.afw.image.PARENT))
        func = model.asMultiShapelet(self.center)
        func.evaluate().addToImage(fit)
        residuals = lsst.afw.image.ImageD(data, True)
        residuals -= fit
        residuals.getArray()[:,:] *= 10.0
        outerEllipse = lsst.afw.geom.ellipses.Quadrupole(model.ellipse)
        outerEllipse.scale(model.radiusRatio)
        ellipses = [
            lsst.afw.geom.ellipses.Ellipse(model.ellipse, self.center),
            lsst.afw.geom.ellipses.Ellipse(outerEllipse, self.center),
            ]
        pyplot.clf()
        pyplot.subplot(1,3,1)
        vmin, vmax = self._plotImage(data, title="original PSF model image", ellipses=ellipses)
        pyplot.subplot(1,3,2)
        self._plotImage(fit, title="shapelet fit", ellipses=ellipses, vmin=vmin, vmax=vmax)
        pyplot.subplot(1,3,3)
        self._plotImage(residuals, title="10*(original-shapelets)", ellipses=ellipses, vmin=vmin, vmax=vmax)
        cax = pyplot.axes([0.1, 0.05, 0.8, 0.025])
        pyplot.colorbar(orientation="horizontal", cax=cax)

class FitProfileViewer(ViewerBase):

    Model = ms.FitProfileModel
    Control = ms.FitProfileControl
    Algorithm = ms.FitProfileAlgorithm

    def __init__(self, source, exposure, name="exp", psfCtrl=None):
        if name == "exp":
            config = ms.FitExponentialConfig()
        else:
            config = ms.FitDeVaucouleurConfig()           
        if psfCtrl is None:
            psfCtrl = ms.FitPsfControl()
        self.ctrl = config.makeControl()
        self.ctrl.name = "multishapelet." + name
        self.psfModel = ms.FitPsfModel(psfCtrl, source)
        self.saved = self.Model(self.ctrl, source)
        self.center = source.getCentroid()
        self.shape = source.getShape()
        badPixelMask = lsst.afw.image.MaskU.getPlaneBitMask(self.ctrl.badMaskPlanes)
        self.inputs = ms.ModelInputHandler(exposure.getMaskedImage(), self.center, 
                                           source.getFootprint(), self.ctrl.growFootprint,
                                           badPixelMask, self.ctrl.usePixelWeights)
        self.footprint = self.inputs.getFootprint()
        self.bbox = self.footprint.getBBox()
        self.image = lsst.afw.image.ImageD(self.bbox)
        self.image.getArray()[:,:] = exposure.getMaskedImage().getImage().getArray()[self.bbox.getSlices()]
        opt = self.Algorithm.makeOptimizer(self.ctrl, self.psfModel, self.shape, self.inputs)
        maxIter = opt.getControl().maxIter
        self.iterations = [self.Iteration(opt, self)]
        for self.iterCount in range(maxIter):
            opt.step()
            self.iterations.append(self.Iteration(opt, self))
            if opt.getState() & ms.HybridOptimizer.FINISHED:
                break
        self.model = self.Model(self.iterations[-1].model)
        if self.ctrl.usePsfShapeletTerms:
            self.Algorithm.fitShapeletTerms(self.ctrl, self.psfModel, self.inputs, self.model)


    def plot(self, n=-1):
        iteration = self.iterations[n]
        model = iteration.model
        data = self.image
        fit = lsst.afw.image.ImageD(self.bbox)
        func = model.asMultiShapelet(self.center).convolve(self.psfModel.asMultiShapelet())
        func.evaluate().addToImage(fit)
        residuals = lsst.afw.image.ImageD(data, True)
        residuals -= fit
        ellipses = [
            lsst.afw.geom.ellipses.Ellipse(model.ellipse, self.center),
            ]
        wData = lsst.afw.image.ImageD(self.bbox)
        wData.getArray()[:,:] = float("NaN")
        lsst.afw.detection.expandArray(self.footprint, self.inputs.getData(), wData.getArray(),
                                       self.bbox.getBegin())
        wResiduals = lsst.afw.image.ImageD(self.bbox)
        wResiduals.getArray()[:,:] = float("NaN")
        lsst.afw.detection.expandArray(self.footprint, iteration.func, wResiduals.getArray(),
                                       self.bbox.getBegin())
        wFit = lsst.afw.image.ImageD(wResiduals, True)
        wFit += wData
        wResiduals *= -1.0 # want to plot (data - fit); objective computes (fit - data)
        pyplot.clf()
        pyplot.subplot(2,3,1)
        vmin, vmax = self._plotImage(data, title="data", ellipses=ellipses)
        pyplot.subplot(2,3,2)
        self._plotImage(fit, title="fit", ellipses=ellipses, vmin=vmin, vmax=vmax)
        pyplot.subplot(2,3,3)
        self._plotImage(residuals, title="data-fit", ellipses=ellipses, vmin=vmin, vmax=vmax)
        cax = pyplot.axes([0.1, 0.5, 0.8, 0.025])
        pyplot.colorbar(orientation="horizontal", cax=cax)

        pyplot.subplot(2,3,4)
        vmin, vmax = self._plotImage(wData, ellipses=ellipses)
        pyplot.subplot(2,3,5)
        self._plotImage(wFit, ellipses=ellipses, vmin=vmin, vmax=vmax)
        pyplot.subplot(2,3,6)
        self._plotImage(wResiduals, ellipses=ellipses, vmin=vmin, vmax=vmax)

        cax = pyplot.axes([0.1, 0.05, 0.8, 0.025])
        pyplot.colorbar(orientation="horizontal", cax=cax)
