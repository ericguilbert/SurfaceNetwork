
import os

from terrain import Terrain
from extendedsurfacenetwork import ExtendedSurfaceNetwork
from iomodule import writeGpkgNode, writeGpkgRidge, writeGpkgThalweg, writeGpkgPuddle

directory = "the directory where the DTM is located"
filename = directory + 'name of the DTM file'

directory = "/home/eric/Terrain/Berev_brule/"
filename = directory + "MNT_BEREV_BRULE.tif"

terrain = Terrain(filename)

network = ExtendedSurfaceNetwork()
network.buildESNFromTerrain(terrain, sd = True, smooth = True)

network.assignStrahlerOrder(11000)

#writeGpkgNode(network.nodedict, terrain, directory, 'node2')
#writeGpkgRidge(network.ridgedict, terrain, directory, 'ridge2')
#writeGpkgPuddle(network.puddledict, network.nodedict, terrain, directory, 'puddle2')
#writeGpkgThalweg(network.thalwegdict, terrain, directory, 'thalweg2')
