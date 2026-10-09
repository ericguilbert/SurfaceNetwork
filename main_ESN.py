
import os
from time import perf_counter_ns

from terrain import Terrain
from extendedsurfacenetwork import ExtendedSurfaceNetwork
from depressiongraph import DepressionGraph
from iomodule import writeGpkgNode, writeGpkgLine, writeGpkgFlowLine

os.environ['GDAL_DATA'] = os.environ['CONDA_PREFIX'] + r'\Library\share\gdal'
os.environ['PROJ_LIB'] = os.environ['CONDA_PREFIX'] + r'\Library\share'

directory = "Pathtofolder/"
filename = directory + 'dtm.tif'

# read the DTM
terrain = Terrain(filename)

bigstart = perf_counter_ns()
network = ExtendedSurfaceNetwork()
# build the ESN with all thalwegs and divides
network.buildESNFromTerrain(terrain, sd = True, smooth = False)
network.addThalwegSlopeInDegree()
network.addDivideSlopeInDegree()

start = perf_counter_ns()
# initiate depression graph with all depressions of the ESN
depgraph = DepressionGraph(network)
# merge those which share the same outlet
depgraph.buildPuddleGraph()
# get the thalwegs below the outlets
puddledict = depgraph.getThalwegsInsidePuddles()
end = perf_counter_ns()
print("Puddle computation time:", end - start)
start = perf_counter_ns()
# these thalwegs flow towards the outlet
network.assignFlowDirection(puddledict)
network.accumulateFlow()
end = perf_counter_ns()
print("Flow direction and accumulation computation time:", end - start)

start = perf_counter_ns()
network.assignStrahlerOrder(11000)
network.assignHortonOrder()
end = perf_counter_ns()
print("Strahler and Horton order computation time:", end - start)
print("Total time = %.3f seconds" % ((end - bigstart)*10**-9))
network.saveNetwork(directory, 'network')

writeGpkgNode(network.nodedict, terrain, directory, 'node')
writeGpkgLine(network.linedict, terrain, directory, 'line')
writeGpkgFlowLine(network.linedict, terrain, directory, 'flowline')
