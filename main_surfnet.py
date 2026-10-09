import os

from terrain import Terrain
from surfacenetwork import SurfaceNetwork
#from iomodule import writeShpCandidateSaddle, writeShpSaddle
from iomodule import writeGpkgNode, writeGpkgLine

directory = "Pathtofolder/"
filename = directory + 'rasterDTM.tif'

terrain = Terrain(filename)

network = SurfaceNetwork()
network.buildFromTerrain(terrain, divide = True)

network.saveNetwork(directory, 'sn')
writeGpkgNode(network.nodedict, terrain, directory, 'node')
writeGpkgLine(network.linedict, terrain, directory, 'line')
