
from terrain import Terrain
from thalwegnetwork import ThalwegNetwork

directory = "folder name"
filename = directory + 'file name'

print(filename)
terrain = Terrain(filename)

network = ThalwegNetwork()
network.buildThalwegsFromTerrain(terrain, False)
