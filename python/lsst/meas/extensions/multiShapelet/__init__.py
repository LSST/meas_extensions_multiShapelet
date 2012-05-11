from .multiShapeletLib import *   # generated by SWIG
from .version import *   # generated by sconsUtils unless you tell it not to

import lsst.pex.config
import lsst.meas.algorithms

@lsst.pex.config.wrap(HybridOptimizerControl)
class HybridOptimizerConfig(lsst.pex.config.Config):
    pass

@lsst.pex.config.wrap(FitProfileControl)
class FitProfileConfig(lsst.pex.config.Config):
    pass

class FitExponentialConfig(FitProfileConfig):
    def setDefaults(self):
        self.profile = "tractor-exponential"

class FitDeVaucouleurConfig(FitProfileConfig):
    def setDefaults(self):
        self.profile = "tractor-devaucouleur"

lsst.meas.algorithms.AlgorithmRegistry.register("multishapelet.psf", FitPsfControl)
lsst.meas.algorithms.AlgorithmRegistry.register("multishapelet.exp", FitProfileControl, FitExponentialConfig)
lsst.meas.algorithms.AlgorithmRegistry.register("multishapelet.dev", FitProfileControl, FitDeVaucouleurConfig)

def loadProfiles():
    import os
    import cPickle
    root = os.environ["MEAS_EXTENSIONS_MULTISHAPELET_DIR"]
    for filename in ("tractor.p",):
        with open(os.path.join(root, "data", filename), "r") as f:
            d = cPickle.load(f)
            for k, v in d.iteritems():
                flux, radius = v
                MultiGaussianRegistry.insert(k, flux, radius, True)
loadProfiles()

# cleanup namespace
del loadProfiles
del lsst
