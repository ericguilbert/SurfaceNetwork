
from terrain import Terrain
from surfacenetwork import SurfaceNetwork
from iomodule import writeGpkgNode, writeGpkgRidge, writeGpkgThalweg


directory = "the directory where the DTM is located"
filename = directory + 'name of the DTM file'

print(filename)
terrain = Terrain(filename)

network = SurfaceNetwork()
network.buildFromTerrain(terrain)

writeGpkgNode(network.nodedict, terrain, directory, 'node')
writeGpkgThalweg(network.thalwegdict, terrain, directory, 'thalweg')
writeGpkgRidge(network.ridgedict, terrain, directory, 'ridge')

