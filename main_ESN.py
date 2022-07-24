
import os

from terrain import Terrain
from extendedsurfacenetwork import ExtendedSurfaceNetwork
#from iomodule import writeShpCandidateSaddle, writeShpSaddle
from iomodule import writeGpkgNode, writeGpkgRidge, writeGpkgThalweg, writeGpkgMultiPoint

directory = "path to the DTM"
filename = directory + "name of the DTM"

terrain = Terrain(filename)

network = ExtendedSurfaceNetwork()
network.buildESNFromTerrain(terrain, sd = True, smooth = False)

network.assignStrahlerOrder(11000)

writeGpkgNode(network.nodedict, directory, 'node', terrain)
writeGpkgRidge(network.ridgedict, directory, 'ridge', terrain)
writeGpkgMultiPoint(network.puddledict, network.nodedict, directory, 'puddle', terrain)
writeGpkgThalweg(network.thalwegdict, directory, 'thalweg', terrain)
