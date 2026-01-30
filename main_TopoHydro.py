# -*- coding: utf-8 -*-
"""
Created on Mon Jul 17 17:45:22 2023

@author: Yassmine Zada
"""

from time import perf_counter_ns
from terrain import Terrain

from iomodule import writeGpkgLine
from thalwegnetwork import ThalwegNetwork
#from iomodule import writeGpkgMultiPoint, writeGpkgLine, writeGpkgThalweg, writeGpkgNode, writeGpkgBasinsC, writeGpkgPuddle, writeGpkgRidge
from iomodule import writeGpkgMultiPoint, writeGpkgLine, writeGpkgThalweg, writeGpkgNode, writeGpkgPuddle, writeGpkgRidge
from topohydronetwork import TopoHydroNetwork, readBDculverts

start = perf_counter_ns()


directory = "./Data/"
filename = directory + 'MNT_BEREV_p.tif'
print("Loading the DTM")
terrain = Terrain(filename)

culvert_filename = directory + "ponceaux_berev_p.gpkg"
culvert = readBDculverts(culvert_filename, terrain)

topohydro_network = TopoHydroNetwork()
topohydro_network.buildTopoHydroFromTerrain(terrain, culvert)

writeGpkgNode(topohydro_network.nodedict, terrain, directory + 'output/', 'nodes')
writeGpkgThalweg(topohydro_network.thalwegdict, terrain, directory + 'output/', 'thalwegs')
writeGpkgThalweg(topohydro_network.culvertdict, terrain, directory + 'output/', 'culverts')
writeGpkgThalweg(topohydro_network.flowdict, terrain, directory + 'output/', 'flows')
writeGpkgRidge(topohydro_network.ridgedict, terrain, directory + 'output/', 'ridges')
# writeGpkgPuddle(topohydro_network.puddledict, topohydro_network.nodedict, terrain, directory + 'output/', 'puddles')

end = perf_counter_ns()
elapsed_time = (end - start)/1_000_000_000

print(f"Algorithm execution time: {elapsed_time:.5f} secondes")

