import os

from terrain import Terrain
from surfacenetwork import SurfaceNetwork
from iomodule import writeGpkgNode, writeGpkgLine

directory = "Pathtofolder/"
filename = directory + 'rasterDTM.tif'

terrain = Terrain(filename)

network = SurfaceNetwork()
network.buildFromTerrain(terrain, divide = False)

writeGpkgNode(network.nodedict, terrain, directory, 'node')
writeGpkgLine(network.linedict, terrain, directory, 'line')
