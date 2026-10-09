# -*- coding: utf-8 -*-
"""
Created on Mon Jul 17 17:45:22 2023

@author: Yassmine Zada
"""

from time import perf_counter_ns
from terrain import Terrain

from iomodule import writeGpkgLine, writeGpkgFlowLine
from iomodule import writeGpkgNode
from topohydronetwork import TopoHydroNetwork, readBDculverts
from depressiongraph import DepressionGraph

start = perf_counter_ns()


directory = "./Data/"
filename = directory + 'MNT_BEREV_p.tif'
print("Loading the DTM")
terrain = Terrain(filename)
terrain.crs = {'init' :'epsg:6622'}

bigstart = perf_counter_ns()
# read the culvert file
culvert_filename = directory + "ponceaux_berev_p.gpkg"
culvert = readBDculverts(culvert_filename, terrain)
length = 4
culvert = terrain.pixelCulverts(culvert, length)

network = TopoHydroNetwork()
network.buildTopoHydroFromTerrain(terrain, culvert)
start = perf_counter_ns()
depgraph = DepressionGraph(network, True)
depgraph.buildPuddleGraph()
puddledict = depgraph.getThalwegsInsidePuddles()
end = perf_counter_ns()
print("Puddle computation time:", end - start)
start = perf_counter_ns()
network.assignFlowDirection(puddledict)
network.accumulateFlow()
end = perf_counter_ns()
print("Flow direction and accumulation computation time:", end - start)
start = perf_counter_ns()
network.assignStrahlerOrder(11000)
network.assignHortonOrder()
end = perf_counter_ns()
print("Strahler and Horton order computation time:", end - start)

elapsed_time = (end - bigstart)/1_000_000_000

print(f"Algorithm execution time: {elapsed_time:.5f} secondes")
print("Saving files")
writeGpkgNode(network.nodedict, terrain, directory + 'output/', 'nodes')
writeGpkgLine(network.linedict, terrain, directory + 'output/', 'lines')
writeGpkgFlowLine(network.linedict, terrain, directory + 'output/', 'flowline')

