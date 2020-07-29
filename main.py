
from terrain import Terrain
from surfacenetwork import SurfaceNetwork
from iomodule import writeShpNode, writeShpLine


directory = "the directory where the DTM is located"
filename = directory + 'name of the DTM file'

print(filename)
terrain = Terrain(filename)

network = SurfaceNetwork()
network.buildFromTerrain(terrain)

writeShpNode(network.nodedict, directory+'node', terrain)
writeShpLine(network.thalwegdict, directory+'thalweg', terrain)
writeShpLine(network.ridgedict, directory+'ridge', terrain)

