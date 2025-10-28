# coding=utf-8
"""
iomodule.py:
    Several functions used to load a raster DTM and save the different dictionaries composing the surface network into geopackages
    
    author: Eric Guilbert
"""
import os
import numpy as np
import pickle
import shapely
import json
os.environ["USE_PATH_FOR_GDAL_PYTHON"] = "YES"
from osgeo import ogr, osr, gdal, gdal_array


# read a raster DTM from a file
def readRasterDTM(filename, terrain):
    raster = gdal.Open(filename)
    
    if not raster:
        print("Layer failed to load!")
    terrain.crs = raster.GetProjection()

    # need to test if the provider is a raster file
    terrain.n = raster.RasterXSize  # number of columns
    terrain.m = raster.RasterYSize  # number of rows
    rastersize = raster.GetGeoTransform()
    terrain.upper_left_x = rastersize[0]
    terrain.upper_left_y = rastersize[3]
    terrain.x_size = rastersize[1]
    terrain.y_size = rastersize[5] # sign to be checked
    
    # get nodata value or set it below min
    band = raster.GetRasterBand(1)
    minvalue = band.GetMinimum()
    maxvalue = band.GetMaximum()
    nodata = band.GetNoDataValue()
    if minvalue is None or maxvalue is None:
        band.ComputeStatistics(0)
        minvalue = band.GetMinimum()
        maxvalue = band.GetMaximum()
        nodata = band.GetNoDataValue()
    print("in readRasterDTM", nodata, minvalue)
    if not nodata or np.isnan(nodata) or nodata > minvalue:
        nodata = minvalue - 1
    else:
        if minvalue < nodata:
            print("in readRasterDTM: warning! minz < nodata")
    print(terrain.m, 'rows', terrain.n, 'columns')
    print("minvalue =", minvalue, 'nodata =', nodata)

    terrain.nodata = nodata
    terrain.dtm = raster.ReadAsArray()

# read a polygon from a geopackage
def readGpkgPolygon(directory, filename):
    driver = ogr.GetDriverByName("GPKG")
    filename = directory+filename+".gpkg"
    if os.path.exists(filename):
        polygonfile = driver.Open(filename, 0)
    polygonlayer = polygonfile.GetLayer()
    layerdef = polygonlayer.GetLayerDefn()
    for i in range(0, layerdef.GetFieldCount()):
        fielddef = layerdef.GetFieldDefn(i)
        fieldname = fielddef.GetName()
        print(fieldname)
    geomlist = []
    for feature in polygonlayer:
        tmpgeom = feature.GetGeometryRef()
        geom = tmpgeom.Clone()
        geomlist.append(geom)
        print(type(geom))

    filename = None
    polygonlayer = None
    return geomlist
    
# write a line geopackage from a dict
def writeGpkgLine(linedict, terrain, directory, linename):
    driver = ogr.GetDriverByName("GPKG")
    filename = directory+linename+".gpkg"
    if os.path.exists(filename):
        driver.DeleteDataSource(filename)
    linefile = driver.CreateDataSource(filename)
    crs = osr.SpatialReference(terrain.crs)
    
    linelayer = linefile.CreateLayer(linename, crs, ogr.wkbLineString)

    linelayer.CreateField(ogr.FieldDefn("id", ogr.OFTInteger))
    #linelayer.CreateField(ogr.FieldDefn("start", ogr.OFTInteger))
    #linelayer.CreateField(ogr.FieldDefn("end", ogr.OFTInteger))
    layerdef = linelayer.GetLayerDefn()

    linelayer.StartTransaction()
    for id, line in linedict.items():
        feature = ogr.Feature(layerdef)
        feature.SetField("id", id)
        #feature.SetField("start", line['start'])
        #feature.SetField("end", line['end'])

        polyline = ogr.Geometry(ogr.wkbLineString)
        #for p in line['polyline']:
        for p in line:
            x, y = terrain.fromIndexToCoordinates(p[0],p[1])
            polyline.AddPoint(x, y)

        feature.SetGeometry(polyline)
        
        linelayer.CreateFeature(feature)
        feature = None
    linelayer.CommitTransaction()
    linefile = None
    linelayer = None

def writeGpkgPolygon(pointlist, terrain, directory, polygonname):
    driver = ogr.GetDriverByName("GPKG")
    filename = directory+polygonname+".gpkg"
    if os.path.exists(filename):
        driver.DeleteDataSource(filename)
    polygonfile = driver.CreateDataSource(filename)
    crs = osr.SpatialReference(terrain.crs)
    
    polygonlayer = polygonfile.CreateLayer(polygonname, crs, ogr.wkbPolygon)

    #polygonlayer.CreateField(ogr.FieldDefn("thalweg", ogr.OFTInteger))
    #polygonlayer.CreateField(ogr.FieldDefn("order", ogr.OFTInteger))
    layerdef = polygonlayer.GetLayerDefn()

    polygonlayer.StartTransaction() 
    feature = ogr.Feature(layerdef)
    #feature.SetField("thalweg", thalweg)
    #feature.SetField("order", order)
    ring = ogr.Geometry(ogr.wkbLinearRing)
    for p in pointlist:
        x, y = terrain.fromIndexToCoordinates(p[0],p[1])
        ring.AddPoint(x, y)

    polygon = ogr.Geometry(ogr.wkbPolygon)
    polygon.AddGeometry(ring)        
    feature.SetGeometry(polygon)
    
    polygonlayer.CreateFeature(feature)
    feature = None
    polygonlayer.CommitTransaction()
    polygonfile = None
    polygonlayer = None

def writeGpkgMultiPoint(pointlist, terrain, directory, multipointname):
    driver = ogr.GetDriverByName("GPKG")
    filename = directory+multipointname+".gpkg"
    if os.path.exists(filename):
        driver.DeleteDataSource(filename)
    multipointfile = driver.CreateDataSource(filename)
    crs = osr.SpatialReference(terrain.crs)
    
    multipointlayer = multipointfile.CreateLayer(multipointname, crs, ogr.wkbMultiPoint)

    layerdef = multipointlayer.GetLayerDefn()

    multipointlayer.StartTransaction() 
    feature = ogr.Feature(layerdef)
    multipoint = ogr.Geometry(ogr.wkbMultiPoint)
    for p in pointlist:
        x, y = terrain.fromIndexToCoordinates(p[0],p[1])
        pt = ogr.Geometry(ogr.wkbPoint)
        pt.AddPoint(x,y)
        multipoint.AddGeometry(pt)

    feature.SetGeometry(multipoint)
    
    multipointlayer.CreateFeature(feature)
    feature = None
    multipointlayer.CommitTransaction()
    multipointfile = None
    multipointlayer = None

def writeGpkgMultiPolygonFromMultiPoint(pointlistdict, terrain, directory, multipolygonname):
    driver = ogr.GetDriverByName("GPKG")
    filename = directory+multipolygonname+".gpkg"
    if os.path.exists(filename):
        driver.DeleteDataSource(filename)
    multipolygonfile = driver.CreateDataSource(filename)
    crs = osr.SpatialReference(terrain.crs)
    
    multipolygonlayer = multipolygonfile.CreateLayer(multipolygonname, crs, ogr.wkbMultiPolygon)
    multipolygonlayer.CreateField(ogr.FieldDefn("id", ogr.OFTInteger))

    layerdef = multipolygonlayer.GetLayerDefn()
    res = 0.5*terrain.x_size

    multipolygonlayer.StartTransaction()
    for i, pointlist in pointlistdict.items():
        feature = ogr.Feature(layerdef)
        feature.SetField("id", i)

        multipoint = [terrain.fromIndexToCoordinates(p[0],p[1]) for p in pointlist]
        multipoint = shapely.geometry.MultiPoint(multipoint)
        plgn = shapely.buffer(multipoint, res, cap_style = "square")   
        plgn = ogr.CreateGeometryFromWkb(plgn.wkb)
        feature.SetGeometry(plgn)
        
        multipolygonlayer.CreateFeature(feature)
        feature = None
    multipolygonlayer.CommitTransaction()
    multipolygonfile = None
    multipolygonlayer = None

def writeGpkgNode(nodedict, terrain, directory, nodename):
    driver = ogr.GetDriverByName("GPKG")
    filename = directory+nodename+".gpkg"
    if os.path.exists(filename):
        driver.DeleteDataSource(filename)
    nodefile = driver.CreateDataSource(filename)
    crs = osr.SpatialReference(terrain.crs)
    
    nodelayer = nodefile.CreateLayer(nodename, crs, ogr.wkbPoint)
    
    nodelayer.CreateField(ogr.FieldDefn("id", ogr.OFTInteger))
    nodelayer.CreateField(ogr.FieldDefn("i", ogr.OFTReal))
    nodelayer.CreateField(ogr.FieldDefn("j", ogr.OFTReal))
    nodelayer.CreateField(ogr.FieldDefn("z", ogr.OFTReal))
    fieldt = ogr.FieldDefn("thalweg", ogr.OFTString)
    fieldt.SetWidth(128)
    nodelayer.CreateField(fieldt)
    fieldr = ogr.FieldDefn("ridge", ogr.OFTString)
    fieldr.SetWidth(128)
    nodelayer.CreateField(fieldr)
    nodelayer.CreateField(ogr.FieldDefn("type", ogr.OFTInteger))
    layerdef = nodelayer.GetLayerDefn()
    
    virtualpit = nodedict.pop(-1, None)
    nodelayer.StartTransaction()
    for id, pss in nodedict.items():
        (i, j) = pss['ij']
        x, y = terrain.fromIndexToCoordinates(i, j)
        thalwegtxt = ', '.join(str(x) for x in pss['thalweg'])
        ridgetxt = ', '.join(str(x) for x in pss['ridge'])

        feature = ogr.Feature(layerdef)
        feature.SetField("id", id)
        feature.SetField("i", i)
        feature.SetField("j", j)
        feature.SetField("z", float(pss["z"]))
        feature.SetField("type", pss['type'])
        feature.SetField("thalweg", thalwegtxt)
        feature.SetField("ridge", ridgetxt)
        
        wkt = "POINT(%f %f)" %(x,y)
        point = ogr.CreateGeometryFromWkt(wkt)
        feature.SetGeometry(point)

        nodelayer.CreateFeature(feature)
        feature = None
    nodelayer.CommitTransaction()
    if virtualpit:
        nodedict[-1] = virtualpit
    nodefile = None
    nodelayer = None

# write shapefiles with saddles, ridge and thalweg initial segments
def writeShpCandidateSaddle(saddledict, name, terrain):
    """
    Records the saddles computed by terrain.computeSaddles() in two files:
        a shapefile for the nodes and a shapefile for thalweg and ridge segments

    Parameters
    ----------
    saddledict : dictionary of saddles
        Saddles are described by their (i,j) coordinates, their z and their
        lines. Lines are segments defined by two points (first is the saddle), 
        a type (thalweg or ridge) and the z of the second point
    name : string
        the name of the shapefile.
    terrain : a terrain class instance recording the dtm
        Used only for the spatial reference system.

    Returns
    -------
    None.

    """
    driver = ogr.GetDriverByName("ESRI Shapefile")

    saddlename = name + 'node'
    filename = saddlename+".shp"
    if os.path.exists(filename):
        driver.DeleteDataSource(filename)
    saddlefile = driver.CreateDataSource(filename)
    crs = osr.SpatialReference(terrain.crs)

    print(saddlename, filename)
    saddlelayer = saddlefile.CreateLayer(saddlename, crs, ogr.wkbPoint)
    
    saddlelayer.CreateField(ogr.FieldDefn("id", ogr.OFTInteger))
    saddlelayer.CreateField(ogr.FieldDefn("i", ogr.OFTReal))
    saddlelayer.CreateField(ogr.FieldDefn("j", ogr.OFTReal))
    #saddlelayer.CreateField(ogr.FieldDefn("z", ogr.OFTReal))

    for id in saddledict:
        node = saddledict[id]
        x, y = terrain.fromIndexToCoordinates(
            saddledict[id]['ij'][0], saddledict[id]['ij'][1])
        feature = ogr.Feature(saddlelayer.GetLayerDefn())
        feature.SetField("id", id)
        feature.SetField("i", node['ij'][0])
        feature.SetField("j", node['ij'][1])
        #feature.SetField("z", node['z'])

        wkt = "POINT(%f %f)" %(x,y)
        point = ogr.CreateGeometryFromWkt(wkt)
        feature.SetGeometry(point)

        saddlelayer.CreateFeature(feature)
        feature = None
    saddlefile = None
    saddlelayer = None

    ridgename = name + 'line'
    filename = ridgename+".shp"
    if os.path.exists(filename):
        driver.DeleteDataSource(filename)
    ridgefile = driver.CreateDataSource(filename)
    crs = osr.SpatialReference(terrain.crs)

    ridgelayer = ridgefile.CreateLayer(saddlename, crs, ogr.wkbLineString)

    ridgelayer.CreateField(ogr.FieldDefn("i", ogr.OFTInteger))
    ridgelayer.CreateField(ogr.FieldDefn("j", ogr.OFTInteger))
    #ridgelayer.CreateField(ogr.FieldDefn("z", ogr.OFTReal))
    ridgelayer.CreateField(ogr.FieldDefn("type", ogr.OFTInteger))

    for id in saddledict:
        node = saddledict[id]
        for l in node['line']:
            feature = ogr.Feature(ridgelayer.GetLayerDefn())
            feature.SetField("i", node['ij'][0])
            feature.SetField("j", node['ij'][1])
            #feature.SetField("z", l[3])
            feature.SetField('type', l[2])
            polyline = ogr.Geometry(ogr.wkbLineString)
            x, y = terrain.fromIndexToCoordinates(l[0][0], l[0][1])
            polyline.AddPoint(x, y)
            x, y = terrain.fromIndexToCoordinates(l[1][0], l[1][1])
            polyline.AddPoint(x, y)
            feature.SetGeometry(polyline)
            
            ridgelayer.CreateFeature(feature)
    ridgefile = None
    ridgelayer = None

def writeGpkgPuddle(puddledict, terrain, nodedict, directory, puddlename):
    """
    Parameters
    ----------
    puddledict : TYPE
        DESCRIPTION.
    puddlename : TYPE
        DESCRIPTION.
    terrain : TYPE
        DESCRIPTION.

    Returns
    -------
    None.

    """
    driver = ogr.GetDriverByName("GPKG")
    filename = directory+puddlename+".gpkg"
    if os.path.exists(filename):
        driver.DeleteDataSource(filename)
    puddlefile = driver.CreateDataSource(filename)
    crs = osr.SpatialReference(terrain.crs)
    
    puddlelayer = puddlefile.CreateLayer(puddlename, crs, ogr.wkbMultiPoint)

    puddlelayer.CreateField(ogr.FieldDefn("id", ogr.OFTInteger))
    puddlelayer.CreateField(ogr.FieldDefn("outlet", ogr.OFTInteger))
#    puddlelayer.CreateField(ogr.FieldDefn("outflow", ogr.OFTString))
    layerdef = puddlelayer.GetLayerDefn()

    puddlelayer.StartTransaction()
    for id in puddledict:
        puddle = puddledict[id]
#        outlettext = ', '.join(str(x) for x in puddle['outlet'])
#        outflowtext = ', '.join(str(x) for x in puddle['outflow'])
        
        feature = ogr.Feature(layerdef)
        feature.SetField("id", id)
        feature.SetField("outlet", puddle['outlet'])
#        feature.SetField("outflow", outflowtext)
        multipoint = ogr.Geometry(ogr.wkbMultiPoint)
        for ip in puddle['nodes']:
            if ip != -1:
                p = nodedict[ip]['ij']
                x, y = terrain.fromIndexToCoordinates(p[0],p[1])
                pt = ogr.Geometry(ogr.wkbPoint)
                pt.AddPoint(x,y)
                multipoint.AddGeometry(pt)

        feature.SetGeometry(multipoint)
        
        puddlelayer.CreateFeature(feature)
        feature = None
    puddlelayer.CommitTransaction()
    puddlefile = None
    puddlelayer = None

def writeGpkgThalweg(linedict, terrain, directory, linename):
    driver = ogr.GetDriverByName("GPKG")
    filename = directory+linename+".gpkg"
    if os.path.exists(filename):
        driver.DeleteDataSource(filename)
    linefile = driver.CreateDataSource(filename)
    crs = osr.SpatialReference(terrain.crs)
    
    linelayer = linefile.CreateLayer(linename, crs, ogr.wkbLineString)

    linelayer.CreateField(ogr.FieldDefn("id", ogr.OFTInteger))
    linelayer.CreateField(ogr.FieldDefn("start", ogr.OFTInteger))
    linelayer.CreateField(ogr.FieldDefn("end", ogr.OFTInteger))
    accumulationfield = False
    orderfield = False
    hortonfield = False
    slopefield = False
    key = list(linedict)[0]
    if 'accumulation' in linedict[key]:
        linelayer.CreateField(ogr.FieldDefn("accumulation", ogr.OFTReal))
        accumulationfield = True
    if 'order' in linedict[key]:
        linelayer.CreateField(ogr.FieldDefn("order", ogr.OFTInteger))
        orderfield = True
    if 'horton' in linedict[key]:
        linelayer.CreateField(ogr.FieldDefn("horton", ogr.OFTInteger))
        hortonfield = True
    if 'slope' in linedict[key]:
        linelayer.CreateField(ogr.FieldDefn("slope", ogr.OFTReal))
        slopefield = True
    layerdef = linelayer.GetLayerDefn()

    linelayer.StartTransaction()    
    for id in linedict:
        line = linedict[id]
        feature = ogr.Feature(layerdef)
        feature.SetField("id", id)
        feature.SetField("start", line['start'])
        feature.SetField("end", line['end'])
        if accumulationfield:
            feature.SetField("accumulation", float(line['accumulation']))
        if orderfield:
            feature.SetField("order", line['order'])
        if hortonfield:
            feature.SetField("horton", line['horton'])
        if slopefield:
            feature.SetField("slope", line['slope'])
        directedline = list(line['polyline'])
        if 'flowdirection' in line and not line['flowdirection']:
            directedline = reversed(directedline)
        polyline = ogr.Geometry(ogr.wkbLineString)
        for p in directedline:
            x, y = terrain.fromIndexToCoordinates(p[0],p[1])
            polyline.AddPoint(x, y)

        feature.SetGeometry(polyline)
        
        linelayer.CreateFeature(feature)
        feature = None
    linelayer.CommitTransaction()
    linefile = None
    linelayer = None

# write a line Shapefile from a dict
def writeGpkgRidge(linedict, terrain, directory, linename):
    driver = ogr.GetDriverByName("GPKG")
    filename = directory+linename+".gpkg"
    if os.path.exists(filename):
        driver.DeleteDataSource(filename)
    linefile = driver.CreateDataSource(filename)
    crs = osr.SpatialReference(terrain.crs)
    
    linelayer = linefile.CreateLayer(linename, crs, ogr.wkbLineString)

    slopefield = False
    key = list(linedict)[0]
    linelayer.CreateField(ogr.FieldDefn("id", ogr.OFTInteger))
    linelayer.CreateField(ogr.FieldDefn("start", ogr.OFTInteger))
    linelayer.CreateField(ogr.FieldDefn("end", ogr.OFTInteger))
    if 'slope' in linedict[key]:
        linelayer.CreateField(ogr.FieldDefn("slope", ogr.OFTReal))
        slopefield = True
    layerdef = linelayer.GetLayerDefn()

    linelayer.StartTransaction()
    for id in linedict:
        line = linedict[id]
        feature = ogr.Feature(layerdef)
        feature.SetField("id", id)
        feature.SetField("start", line['start'])
        feature.SetField("end", line['end'])
        if slopefield:
            feature.SetField("slope", line['slope'])

        polyline = ogr.Geometry(ogr.wkbLineString)
        for p in line['polyline']:
            x, y = terrain.fromIndexToCoordinates(p[0],p[1])
            polyline.AddPoint(x, y)

        feature.SetGeometry(polyline)
        
        linelayer.CreateFeature(feature)
        feature = None
    linelayer.CommitTransaction()
    linefile = None
    linelayer = None

# write a geopackage from a thalweg dict
def writeGpkgThalwegESN(linedict, terrain, directory, linename):
    driver = ogr.GetDriverByName("GPKG")
    filename = directory+linename+".gpkg"
    if os.path.exists(filename):
        driver.DeleteDataSource(filename)
    linefile = driver.CreateDataSource(filename)
    crs = osr.SpatialReference(terrain.crs)
    
    linelayer = linefile.CreateLayer(linename, crs, ogr.wkbLineString)

    linelayer.CreateField(ogr.FieldDefn("id", ogr.OFTInteger))
    linelayer.CreateField(ogr.FieldDefn("start", ogr.OFTInteger))
    linelayer.CreateField(ogr.FieldDefn("end", ogr.OFTInteger))
    linelayer.CreateField(ogr.FieldDefn("accumulation", ogr.OFTReal))
    linelayer.CreateField(ogr.FieldDefn("order", ogr.OFTInteger))
    linelayer.CreateField(ogr.FieldDefn("reversedflow", ogr.OFTInteger))
    linelayer.CreateField(ogr.FieldDefn("dale", ogr.OFTInteger))
    linelayer.CreateField(ogr.FieldDefn("lefthill", ogr.OFTInteger))
    linelayer.CreateField(ogr.FieldDefn("righthill", ogr.OFTInteger))
    hortonfield = False
    slopefield = False
    key = list(linedict)[0]
    if 'horton' in linedict[key]:
        linelayer.CreateField(ogr.FieldDefn("horton", ogr.OFTInteger))
        hortonfield = True
    if 'slope' in linedict[key]:
        linelayer.CreateField(ogr.FieldDefn("slope", ogr.OFTReal))
        slopefield = True
    layerdef = linelayer.GetLayerDefn()

    linelayer.StartTransaction()    
    for id in linedict:
        line = linedict[id]
        feature = ogr.Feature(layerdef)
        feature.SetField("id", id)
        feature.SetField("start", line['start'])
        feature.SetField("end", line['end'])
        feature.SetField("accumulation", float(line['accumulation']))
        feature.SetField("order", line['order'])
        feature.SetField("reversedflow", line['flowdirection'])
        feature.SetField("dale", line['dale'])
        feature.SetField("lefthill", line['lefthill'])
        feature.SetField("righthill", line['righthill'])
        if hortonfield:
            feature.SetField("horton", line['horton'])
        if slopefield:
            feature.SetField("slope", line['slope'])
        directedline = list(line['polyline'])
        polyline = ogr.Geometry(ogr.wkbLineString)
        for p in directedline:
            x, y = terrain.fromIndexToCoordinates(p[0],p[1])
            polyline.AddPoint(x, y)

        feature.SetGeometry(polyline)
        
        linelayer.CreateFeature(feature)
        feature = None
    linelayer.CommitTransaction()
    linefile = None
    linelayer = None
    
# write a geopackage from a ridge dict
def writeGpkgRidgeESN(linedict, terrain, directory, linename):
    driver = ogr.GetDriverByName("GPKG")
    filename = directory+linename+".gpkg"
    if os.path.exists(filename):
        driver.DeleteDataSource(filename)
    linefile = driver.CreateDataSource(filename)
    crs = osr.SpatialReference(terrain.crs)
    
    linelayer = linefile.CreateLayer(linename, crs, ogr.wkbLineString)

    slopefield = False
    key = list(linedict)[0]
    linelayer.CreateField(ogr.FieldDefn("id", ogr.OFTInteger))
    linelayer.CreateField(ogr.FieldDefn("start", ogr.OFTInteger))
    linelayer.CreateField(ogr.FieldDefn("end", ogr.OFTInteger))
    linelayer.CreateField(ogr.FieldDefn("hill", ogr.OFTInteger))
    linelayer.CreateField(ogr.FieldDefn("leftdale", ogr.OFTInteger))
    linelayer.CreateField(ogr.FieldDefn("rightdale", ogr.OFTInteger))
    if 'slope' in linedict[key]:
        linelayer.CreateField(ogr.FieldDefn("slope", ogr.OFTReal))
        slopefield = True
    layerdef = linelayer.GetLayerDefn()

    linelayer.StartTransaction()
    for id in linedict:
        line = linedict[id]
        feature = ogr.Feature(layerdef)
        feature.SetField("id", id)
        feature.SetField("start", line['start'])
        feature.SetField("end", line['end'])
        feature.SetField("hill", line['hill'])
        feature.SetField("leftdale", line['leftdale'])
        feature.SetField("rightdale", line['rightdale'])
        if slopefield:
            feature.SetField("slope", line['slope'])

        polyline = ogr.Geometry(ogr.wkbLineString)
        for p in line['polyline']:
            x, y = terrain.fromIndexToCoordinates(p[0],p[1])
            polyline.AddPoint(x, y)

        feature.SetGeometry(polyline)
        
        linelayer.CreateFeature(feature)
        feature = None
    linelayer.CommitTransaction()
    linefile = None
    linelayer = None

def writeJsonHills(hilldict, directory, hillname):
    filename = directory+hillname+".json"
    with open(filename, 'w') as outfile:
        json.dump(hilldict, outfile)
        
def writeJsonDales(daledict, directory, dalename):
    filename = directory+dalename+".json"
    with open(filename, 'w') as outfile:
        json.dump(daledict, outfile)
        
def writeGpkgBasins(thalweglist, basinlist, network, directory, basinname):
    driver = ogr.GetDriverByName("GPKG")
    filename = directory+basinname+".gpkg"
    if os.path.exists(filename):
        driver.DeleteDataSource(filename)
    polygonfile = driver.CreateDataSource(filename)
    terrain = network.terrain    
    crs = osr.SpatialReference(terrain.crs)
    
    polygonlayer = polygonfile.CreateLayer(basinname, crs, ogr.wkbPolygon)

    polygonlayer.CreateField(ogr.FieldDefn("thalweg", ogr.OFTInteger))
    polygonlayer.CreateField(ogr.FieldDefn("order", ogr.OFTInteger))
    layerdef = polygonlayer.GetLayerDefn()

    polygonlayer.StartTransaction() 
    n = len(thalweglist)
    for i in range(n):
        thalweg = thalweglist[i]
        order = network.thalwegdict[thalweg]['order']
        feature = ogr.Feature(layerdef)
        feature.SetField("thalweg", thalweg)
        feature.SetField("order", order)
        basin = basinlist[i]
        ring = ogr.Geometry(ogr.wkbLinearRing)
        if not basin:
            print('basin', i, "is empty")
            polygonlayer.CreateFeature(feature)
            feature = None
            continue
        for ir in basin:
            if ir > 0:
                directedline = network.ridgedict[ir]['polyline']
            else:
                directedline = network.ridgedict[-ir]['polyline']
                directedline = reversed(directedline)
            for p in directedline:
                x, y = terrain.fromIndexToCoordinates(p[0],p[1])
                ring.AddPoint(x, y)

        polygon = ogr.Geometry(ogr.wkbPolygon)
        polygon.AddGeometry(ring)        
        feature.SetGeometry(polygon)
        
        polygonlayer.CreateFeature(feature)
        feature = None
    polygonlayer.CommitTransaction()
    polygonfile = None
    polygonlayer = None
        
def loadNetwork(directory, name, terrain = None):
    with open(directory + name + '.snw', 'rb') as networkfile:
        network = pickle.load(networkfile)
        network.terrain = terrain
        return network        


