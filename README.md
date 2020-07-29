# SurfaceNetwork
Construction of a surface network from a raster DTM

The code includes the following files:
main.py the main file
terrain.py defines a DTM class that can be triangulated by adding a diagonal in grid cells
thalwegnetwork.py, a class that stores a set of critical points and a network of thalwegs computed a DTM
surfacenetwork.py, a class that inherits from thalwegnetwork and adds a ridge network to compose a complete surface network
deffunction.py, some functions used by previous files
iomodule.py, a few more functions for reading a DTM and saving the surface network in shapefiles

It requires several libraries to run:
numpy
time (for testing only)
shapely.geometry
functools (only the reduce function)
collections (only the deque class)
math (for some basic math operations)
os, ogr, osr, gdal for all i/o operations
