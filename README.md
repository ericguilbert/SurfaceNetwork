# SurfaceNetwork
Construction of a surface network from a raster DTM
This code is an implementation of the algorithm presented in 
https://josis.org/index.php/josis/article/view/681

The code includes the following files:
main.py the main file
terrain.py defines a DTM class that can be triangulated by adding a diagonal in grid cells
thalwegnetwork.py, a class that stores a set of critical points and a network of thalwegs computed a DTM
surfacenetwork.py, a class that inherits from thalwegnetwork and adds a ridge network to build a complete surface network
deffunction.py, some functions used by previous files
iomodule.py, a few more functions for reading a DTM and saving the surface network in shapefiles

The present code can handle convex and star domains. Other domains and non-connex DTM may lead to an error when organising thalwegs around the virtual pit.

It requires several libraries to run:
numpy
time (for testing only)
shapely.geometry
functools (only the reduce function)
collections (only the deque class)
math (for some basic math operations)
os, ogr, osr, gdal for all i/o operations
