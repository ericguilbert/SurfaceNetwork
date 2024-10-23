# SurfaceNetwork
Construction of a surface network from a raster DTM
This code is an implementation of the algorithm presented in 
https://www.doi.org/10.5311/JOSIS.2021.22.681

The code includes the following files:
main.py the main file
terrain.py defines a DTM class that can be triangulated by adding a diagonal in grid cells
thalwegnetwork.py, a class that stores a set of critical points, a network of thalwegs and a set of hills computed from the DTM
surfacenetwork.py, a class that inherits from thalwegnetwork and adds the ridge network and the dales to build a complete surface network
deffunction.py, some functions used by previous files
iomodule.py, a few more functions for reading a DTM and saving the surface network in shapefiles

The present code can handle any connected set but holes are considered as low areas belonging to the set.

It requires several libraries to run:
numpy
time (for testing only)
shapely.geometry
functools (only the reduce function)
collections (only the deque class)
math (for some basic math operations)
os, ogr, osr, gdal for all i/o operations

# Extended Surface Network
Construction of an extended surface network from a raster DTM. The code is the implementation of the algorithm presented in
https://doi.org/10.5311/JOSIS.2023.26.240

Two files are added to build a new network that includes the surface network and the drainage network
streamnetwork.py, a class that inherits from thalwegnetwork, that computes puddles (places where the water accumulates and the flow is interrupted) and adds a flow direction to each thalweg. The flow direction is directed towards the lowest point except in puddles where it is directed towards the outlet.
extendedsurfacenetwork.py, a class that inherits from surfacenetwork and streamnetwork. It adds ridges starting from confluences and pits. Dales are redefined in a one-to-one relationship to thalwegs. Flow direction and dales are used to compute a flow accumulation for each thalweg. Drainage networks can be extracted from the extended surface network by defining an accumulation threshold. Methods are added to compute drainage basins and the Strahler order.
