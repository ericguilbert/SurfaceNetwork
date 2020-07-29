# coding=utf-8
"""
iomodule.py:
    Several functions used to load a raster DTM and save the different dictionaries composing the surface network into shapefiles
    
    author: Eric Guilbert
"""
import os, ogr, osr, gdal

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
    if not nodata or nodata > minvalue:
        nodata = minvalue - 1
    else:
        if minvalue < nodata:
            print("in readRasterDTM: warning! minz < nodata")
    print(terrain.m, 'rows', terrain.n, 'columns')
    print("minvalue =", minvalue, 'nodata =', nodata)

    terrain.nodata = nodata
    terrain.dtm = raster.ReadAsArray()

# write a line Shapefile from a dict
def writeShpLine(linedict, linename, terrain):
    driver = ogr.GetDriverByName("ESRI Shapefile")
    filename = linename+".shp"
    if os.path.exists(filename):
        driver.DeleteDataSource(filename)
    linefile = driver.CreateDataSource(filename)
    crs = osr.SpatialReference(terrain.crs)
    
    linelayer = linefile.CreateLayer(linename, crs, ogr.wkbLineString)

    linelayer.CreateField(ogr.FieldDefn("id", ogr.OFTInteger))
    linelayer.CreateField(ogr.FieldDefn("start", ogr.OFTInteger))
    linelayer.CreateField(ogr.FieldDefn("end", ogr.OFTInteger))

    for id in linedict:
        line = linedict[id]
        feature = ogr.Feature(linelayer.GetLayerDefn())
        feature.SetField("id", id)
        feature.SetField("start", line['start'])
        feature.SetField("end", line['end'])
#        polyline = "LINESTRING("
        polyline = ogr.Geometry(ogr.wkbLineString)
        for p in line['polyline']:
            x, y = terrain.fromIndexToCoordinates(p[0],p[1])
            polyline.AddPoint(x, y)
#            polyline += "%f %f," %(x, y)
#        polyline = polyline[:-1]+")"

        feature.SetGeometry(polyline)
        
        linelayer.CreateFeature(feature)
        feature = None
    linefile = None
    linelayer = None

def writeShpNode(nodedict, nodename, terrain):
    driver = ogr.GetDriverByName("ESRI Shapefile")
    filename = nodename+".shp"
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

    virtualpit = nodedict.pop(-1, None)
    for id, pss in nodedict.items():
        (i, j) = pss['ij']
        x, y = terrain.fromIndexToCoordinates(i, j)
        thalwegtxt = ', '.join(str(x) for x in pss['thalweg'])
        ridgetxt = ', '.join(str(x) for x in pss['ridge'])

        feature = ogr.Feature(nodelayer.GetLayerDefn())
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

def writeShpSaddle(nodedict, name, terrain):
    driver = ogr.GetDriverByName("ESRI Shapefile")
    nodename = name + "node"
    filename = nodename+".shp"
    if os.path.exists(filename):
        driver.DeleteDataSource(filename)
    nodefile = driver.CreateDataSource(filename)
    crs = osr.SpatialReference(terrain.crs)
    
    nodelayer = nodefile.CreateLayer(nodename, crs, ogr.wkbPoint)
    
    nodelayer.CreateField(ogr.FieldDefn("id", ogr.OFTInteger))
    nodelayer.CreateField(ogr.FieldDefn("i", ogr.OFTReal))
    nodelayer.CreateField(ogr.FieldDefn("j", ogr.OFTReal))
    nodelayer.CreateField(ogr.FieldDefn("type", ogr.OFTInteger))

    virtualpit = nodedict.pop(-1, None)
    for id, pss in nodedict.items():
        (i, j) = pss['ij']
        x, y = terrain.fromIndexToCoordinates(i, j)

        feature = ogr.Feature(nodelayer.GetLayerDefn())
        feature.SetField("id", id)
        feature.SetField("i", i)
        feature.SetField("j", j)
        feature.SetField("z", float(pss["z"]))
        feature.SetField("type", pss['type'])
        
        wkt = "POINT(%f %f)" %(x,y)
        point = ogr.CreateGeometryFromWkt(wkt)
        feature.SetGeometry(point)

        nodelayer.CreateFeature(feature)
        feature = None
    if virtualpit:
        nodedict[-1] = virtualpit
    nodefile = None
    nodelayer = None

    ridgename = name + 'line'
    filename = ridgename+".shp"
    if os.path.exists(filename):
        driver.DeleteDataSource(filename)
    ridgefile = driver.CreateDataSource(filename)
    crs = osr.SpatialReference(terrain.crs)

    ridgelayer = ridgefile.CreateLayer(nodename, crs, ogr.wkbLineString)

    ridgelayer.CreateField(ogr.FieldDefn("i", ogr.OFTInteger))
    ridgelayer.CreateField(ogr.FieldDefn("j", ogr.OFTInteger))
    ridgelayer.CreateField(ogr.FieldDefn("type", ogr.OFTInteger))

    for id in nodedict:
        node = nodedict[id]
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

def writeShpMultiPoint(puddledict, nodedict, puddlename, terrain):
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
    driver = ogr.GetDriverByName("ESRI Shapefile")
    filename = puddlename+".shp"
    if os.path.exists(filename):
        driver.DeleteDataSource(filename)
    puddlefile = driver.CreateDataSource(filename)
    crs = osr.SpatialReference(terrain.crs)
    
    puddlelayer = puddlefile.CreateLayer(puddlename, crs, ogr.wkbMultiPoint)

    puddlelayer.CreateField(ogr.FieldDefn("id", ogr.OFTInteger))
    puddlelayer.CreateField(ogr.FieldDefn("outlet", ogr.OFTInteger))
#    puddlelayer.CreateField(ogr.FieldDefn("outflow", ogr.OFTString))

    for id in puddledict:
        puddle = puddledict[id]
#        outlettext = ', '.join(str(x) for x in puddle['outlet'])
#        outflowtext = ', '.join(str(x) for x in puddle['outflow'])
        
        feature = ogr.Feature(puddlelayer.GetLayerDefn())
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
    puddlefile = None
    puddlelayer = None

def writeShpStream(linedict, linename, terrain):
    driver = ogr.GetDriverByName("ESRI Shapefile")
    filename = linename+".shp"
    if os.path.exists(filename):
        driver.DeleteDataSource(filename)
    linefile = driver.CreateDataSource(filename)
    crs = osr.SpatialReference(terrain.crs)
    
    linelayer = linefile.CreateLayer(linename, crs, ogr.wkbLineString)

    linelayer.CreateField(ogr.FieldDefn("id", ogr.OFTInteger))
    linelayer.CreateField(ogr.FieldDefn("start", ogr.OFTInteger))
    linelayer.CreateField(ogr.FieldDefn("end", ogr.OFTInteger))
    linelayer.CreateField(ogr.FieldDefn("order", ogr.OFTInteger))

    for id in linedict:
        line = linedict[id]
        feature = ogr.Feature(linelayer.GetLayerDefn())
        feature.SetField("id", id)
        feature.SetField("start", line['start'])
        feature.SetField("end", line['end'])
        feature.SetField("order", line['order'])
        directedline = list(line['polyline'])
        if not line['flowdirection']:
            directedline = reversed(directedline)
        polyline = ogr.Geometry(ogr.wkbLineString)
        for p in directedline:
            x, y = terrain.fromIndexToCoordinates(p[0],p[1])
            polyline.AddPoint(x, y)
#            polyline += "%f %f," %(x, y)
#        polyline = polyline[:-1]+")"

        feature.SetGeometry(polyline)
        
        linelayer.CreateFeature(feature)
        feature = None
    linefile = None
    linelayer = None
