
import os

from terrain import Terrain
from extendedsurfacenetwork import ExtendedSurfaceNetwork
from iomodule import writeGpkgNode, writeGpkgRidge, writeGpkgThalweg, writeGpkgMultiPoint

directory = "path to the DTM"
filename = directory + "name of the DTM"

terrain = Terrain(filename)

network = ExtendedSurfaceNetwork()
network.buildESNFromTerrain(terrain, sd = True, smooth = False)

network.assignStrahlerOrder(11000)

writeGpkgNode(network.nodedict, terrain, directory, 'node')
writeGpkgRidge(network.ridgedict, terrain, directory, 'ridge')
writeGpkgMultiPoint(network.puddledict, terrain, network.nodedict, directory, 'puddle')
writeGpkgThalweg(network.thalwegdict, terrain, directory, 'thalweg')
