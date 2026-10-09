# Drainage network mapping
Extraction of the drainage network from a raster digital terrain model (DTM). The drainage network is computed directly from the DTM without drainage enforcement (no filling or breaching required).

# Surface network
The construction of the drainage network relies first on the extraction of the surface network (Pfaff, 1976), a Morse-Smale complex where points are the maximums, minimums and saddles of the DTM and lines are ridges (lines of diverging flow) and thalwegs (lines of converging flow). Each maximum is surrounded by a set of thalwegs forming a hill and each minimum is surrounded by a set of ridges forming a depression. 

The code is based on the algorithm presented in https://www.doi.org/10.5311/JOSIS.2021.22.681. The DTM is triangulated to form a piecewise linear surface. This can be done either by maximising the number of saddles or by favouring downhill flow. The first case maximises the number of points, retaining more noise (good to test the robustness of the algorithm but that’s it). The second case reduces the number of pits which are always seen as noise in drainage mapping.

Thalwegs are computed by starting from saddles and following a steepest descent towards a pit. Two paths are possible: we compute the steepest path to one of the neighbours or we compute the steepest descent and move to the closest neighbour. In that case, the offset is recorded and a correction is applied when the offset is larger than a pixel. This is more relevant in terrains with regular slopes to reduce the staircase effect. In high-resolution DTMs, like one-metre resolution DTM obtained from lidar, it does not have much impact.

Ridges are computed second and their position is constrained by thalwegs: a ridge must remain in its hill to avoid crossings. It is still possible that two lines overlap but they are topologically distinct: a ridge can geometrically overlap a thalweg but always remains on the same side of the thalweg. An offset can also be considered when computing the steepest ascent. Hills and depressions are explicitly defined by ordered lists of directed thalwegs or ridges forming a linear ring. 

# Extended surface network
The extended surface network (ESN) adds water divides starting at junctions between thalwegs. It creates a new cell complex where each thalweg is located inside a dale bounded by water divides. The dale corresponds to the direct catchment area of the thalweg. Dales of the ESN partition the depressions of the SN. A depression is formed by all the dales whose thalwegs points to the same pit.

The ESN also adds a flow direction to each thalweg and computes a flow accumulation. The flow inside the depression is directed towards the outlet, the lowest saddle on its border. If two depressions share the same outlet, the depressions are merged into a larger one. In the original algorithm, this was done by defining puddles formed by groups of points below the outlet (https://doi.org/10.5311/JOSIS.2023.26.240). In this current implementation, depressions form a depression graph and are merged until each depression has its own outlet (https://doi.org/10.4230/LIPIcs.COSIT.2024.2). This approach is faster and more robust. The flow accumulation is computed for each thalweg by summing up the area of all dales located upstream. The drainage network can then be extracted by setting an accumulation threshold.

The algorithm can also compute the Horton and the Strahler orders. Thalwegs with an accumulation below the threshold have an order of 0. Drainage basins can be extracted on demand for any thalweg. They are defined by the union of all dales above the thalweg.

# Culvert insertion
The data structure was further extended to culverts into the drainage network (https://doi.org/10.1080/15230406.2026.2633277). The structure follows a new conceptual model where watercourses are separated into surface and underground watercourses. A thalweg is thus a kind of surface watercourse and a culvert is an underground watercourse. Culverts are provided as a table of segments. They are first rasterised and then connected to existing thalwegs. This is done by snapping the culverts to a neighbouring thalweg. A parameter is required to set a snapping distance. The inlet of the culvert is always snapped to a thalweg (located in a gutter along the road). If no thalweg is found at the outlet, a surface watercourse is added and runs down to the next watercourse.

# Class definition
The code defines the following classes:
- Terrain, a DTM class that can be triangulated by adding a diagonal in grid cells
- CellType, inherits from IntFlag. CellType defines the different types of cells that compose the ESN. 
- ThalwegNetwork, a class that stores points (saddles, maximums, minimums and thalweg junctions), a network of thalwegs and a set of hills computed from the DTM
- SurfaceNetwork, the class that defines the surface network as a Morse-Smale complex. It inherits from ThalwegNetwork and adds the ridge network and the depressions to build a complete surface network
- StreamNetwork, a class that inherits from SurfaceNetwork and adds some functions to compute flow directions and accumulation. It can be seen as an abstract class as it mainly adds some functions that are used by the ESN.
- ExtendedSurfaceNetwork, the class that defines the ESN, including the drainage network characteristics.
- DepressionGraph, a graph of depressions used to merge puddles and direct the flow in the ESN.
- TopoHydroNetwork, a class that inherits from ExtendedSurfaceNetwork and adds the culverts. New functions rasterise the culverts and connect them to thalwegs. It also adds surface watercourses (that are not thalwegs).

Over modules are
- deffunction.py, with some functions handling some geometry or structures and used by previous files
- iomodule.py, a few more functions for reading a DTM and saving the surface network in geopackages

It requires several libraries to run:
numpy
time (for testing only)
shapely.geometry
functools (only the reduce function)
collections (only the deque class)
math (for some basic math operations)
os, ogr, osr, gdal, pickle, json, geopandas for i/o operations
