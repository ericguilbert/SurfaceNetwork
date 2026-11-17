
from terrain import Terrain
from thalwegnetwork import ThalwegNetwork
from iomodule import writeGpkgNode, writeGpkgLine

directory = "folder name"
filename = directory + 'file name'

terrain = Terrain(filename)

network = ThalwegNetwork()
network.buildThalwegsFromTerrain(terrain, False)

writeGpkgNode(network.nodedict, terrain, directory, 'node')
writeGpkgLine(network.linedict, terrain, directory, 'line')
