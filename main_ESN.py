
import os

from terrain import Terrain
from extendedsurfacenetwork import ExtendedSurfaceNetwork
from iomodule import writeGpkgNode, writeGpkgRidge, writeGpkgThalweg, writeGpkgPuddle

directory = "C:/Chezmoi/Valley/"
filename = directory + "test2.tif"

terrain = Terrain(filename)

network = ExtendedSurfaceNetwork()
network.buildESNFromTerrain(terrain, sd = True, smooth = True)

network.assignStrahlerOrder(11000)

writeGpkgNode(network.nodedict, terrain, directory, 'node')
writeGpkgRidge(network.ridgedict, terrain, directory, 'ridge')
writeGpkgPuddle(network.puddledict, terrain, network.nodedict, directory, 'puddle')
writeGpkgThalweg(network.thalwegdict, terrain, directory, 'thalweg')
