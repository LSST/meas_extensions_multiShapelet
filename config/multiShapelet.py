# This is a config override, for use with ProcessCcdTask, that enables the MultiShapelet galaxy
# photometry algorithms with their default settings, and sets the ModelFlux source slot to link
# to the flux from the exponential + de Vaucouleur linear combination fit.

import lsst.meas.extensions.multiShapelet

root.measurement.algorithms.names |= lsst.meas.extensions.multiShapelet.algorithms
root.measurement.slots.modelFlux = "multishapelet.combo.flux"
