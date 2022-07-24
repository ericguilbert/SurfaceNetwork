# coding=utf-8
# coding=utf-8
"""
surfacenetwork.py:
    Class defining a surface network.
    
    author: Eric Guilbert
"""

#import time

from time import perf_counter_ns
from math import sqrt
from functools import reduce

import shapely.geometry

import deffunction as df
from thalwegnetwork import ThalwegNetwork

class SurfaceNetwork(ThalwegNetwork):
    # a class describing a surface network composed of
    # nodedict: a dictionary of critical points in the network (pit, peak, saddle, confluence, junction)
    # thalwegdict: a dictionary of thalwegs
    # ridgedict: a dictionary of ridges
    # hilldict and daledict: dictionaries of hills (ascending cells) and dales (descending cells)
    
    # the class inherits from ThalwegNetwork.
    # the constructor first builds the thalweg network and the hills and adds the ridges and dales

    
    def buildFromTerrain(self, terrain, divide = False, correction = True):
        """
        Builds a surface network on a given terrain. Start by detecting saddles.
        Follows with thalwegs and ridges. Ridges are constrained to avoid 
        conflicts with thalwegs.

        Parameters
        ----------
        terrain : Terrain
            Instance of the class Terrain on which the network is built.
        divide : Boolean
            Indicate if ridges have to be extracted in view of watershed computation
            If false (default), they are not and we stick to the basic definition
            If true, a ridge is computed between each pair of consecutive thalwegs
            and ridges are split at confluences
        correction : boolean
            True if a correction is done to stay closer to the gradient
            False if the descent direction is not corrected as in D8

        Returns
        -------
        None.

        """
        if correction:
            df.distancefactor = 2/3
        else:
            df.distancefactor = sqrt(2)/2
        self.buildThalwegsFromTerrain(terrain, correction)
        start = perf_counter_ns()
        self.instantiateRidgesAtSaddles(divide)
        self.traceRidges(divide, correction)
        self.orderRidgesAroundNodes()
        end = perf_counter_ns()
        print("Ridge computation time:", end - start)
        start = perf_counter_ns()
        if divide:
            print("Dales are not computed")
        else:
            self.computeDales()
        end = perf_counter_ns()
        print("Dale computation time:", end - start)
        
    def computeRidgeLength(self, ir):
        """
        Compute the length of a ridge

        Parameters
        ----------
        it : integer, ridge index

        Returns
        -------
        length: length of the ridge.

        """
        r = self.ridge[ir]
        return df.length(r['polyline'])

    def addNewRidge(self, ridgekey, ihill, istart, polylist):
        """
        Add a new ridge to the ridge dict and update node and hill dicts

        Parameters
        ----------
        ridgekey : integer
            key of the new ridge.
        ihill : integer
            key of the hill containing the ridge.
        istart : integer
            key of the starting node.
        polylist : list
            list of pixel coordinates.

        Returns
        -------
        None.

        """
        self.nodedict[istart]['ridge'].append(ridgekey)
        self.hilldict[ihill]['ridge'].append(ridgekey)
        self.ridgedict[ridgekey]={
            'start' : istart,
            'polyline' : polylist,
            'end' : None,
            'hill' : ihill
        }
        
    def computeAscent(self, i, j, v, z, ldr, plgn, danglingconfluence, confluencewedge, thalwegside, dtoverlap, danglingmap):
        """
        Compute the direction of ascent, along the edge closest to the gradient

        Parameters
        ----------
        i : integer
            x coordinate of the ridge point.
        j : integer
            y coordinate of the ridge point.
        v : list of floats
            elevation of neighbour points.
        z : float
            point elevation.
        ldr : list of pairs
            coordinates of neighbour points.
        plgn : shapely polygon
            polygon of the hill containing the ridge.
        danglingconfluence : boolean
            true if the pixel is at the confluence of a dangling thalweg.
        confluencewedge : triplet of points
            wedge formed by the confluence point and the points before and after
            on the thalwegs connecting at the wedge.
        thalwegside : -1, 0, 1
            side of the dangling thalweg on which the ridge is located.
        dtoverlap : integer
            id of the overlapping thalweg.
        danglingmap : array
            map of dangling thalwegs.

        Returns
        -------
        validridge : boolean
            true if no error occured.
        kmax : integer
            index in v of the point of ascent
        inlist : list of points
            points of v located in the hill

        """
        dz = df.gradLength(z, v, ldr)

        maximum = dz[0] - 1
        kmax = 0
        maxassigned = False
        inlist = []
        validridge = True
        
        checkcovers = False
        if self.terrain.pixelclass[i,j] == df.THALWEG or self.terrain.pixelclass[i,j] == df.CONFLUENCE:
            checkcovers = True
        for kk,value in enumerate(dz):
            ij = (i+ldr[kk][0], j+ldr[kk][1])
            #We only take the points that are in the polygon defined by the thalwegs
            if not self.terrain.isInside(ij[0], ij[1]):
                continue

            tside = 0
            # if we leave a dangling thalweg
            if danglingconfluence:
                tside = df.sideOfWedge(confluencewedge, ij) * thalwegside
#                                print('xx', confluencewedge, ij, thalwegside, tside)
            elif dtoverlap > 0 and thalwegside != 0:
                tmppoly = self.thalwegdict[dtoverlap]['polyline']
                kt = tmppoly.index((i,j))
                #if danglingmap[ij] <= 0:
                if danglingmap[ij] != dtoverlap:
                    # get the thalweg point
                    #print('oo',ij, dtoverlap, danglingmap[ij], thalwegside)
                    tside = df.sideOfLine(tmppoly, kt, ij) * thalwegside
                elif ij in tmppoly:
                    kij = tmppoly.index(ij)
                    if kij > kt + 1 or kij < kt - 1:
                        tside = df.sideOfLine(tmppoly, kt, ij) * thalwegside
                        
            covers = True
            try:
                if checkcovers:
                    line = shapely.geometry.LineString([(i, j), (ij[0], ij[1])])
                    covers = plgn.relate_pattern(line, '******FF*')
            except:
                print('in computeAscent topology exception from shapely', i,j)
                validridge = False
            #print(ij, covers, tside)
            if covers and tside >= 0:
                inlist.append(ldr[kk])
                if maxassigned == False:
                    maximum = value
                    kmax = kk
                    maxassigned = True
                elif value>maximum:
                    maximum=value
                    kmax=kk
                elif value==maximum:
                    if df.lexicoSign(ldr[kk][0]-ldr[kmax][0], ldr[kk][1]-ldr[kmax][1])>0:
                        maximum=value
                        kmax=kk

        if maxassigned == False:
            print(i, j, "no next point in createRidge")
            validridge = False
        return validridge, kmax, inlist

    def correctAscent(self, v, z, kmax, ldr, inlist, oldldr, oldlag):
        """
        Correction of the steepest ascent direction according to the offset
        observed at each step.

        Parameters
        ----------
        v : list of float
            Elevation of neighbour pixels.
        z : float
            elevation of the pixel.
        kmax : integer
            index of the ascent direction in the ldr list.
        ldr : list of pairs
            vectors to neighbour pixels.
        inlist : list of integer pairs
            neighbour pixels inside the hill.
        oldldr : pair of integer
            former direction of ascent.
        oldlag : float
            former offset.

        Returns
        -------
        newldr: pair of integer
            corrected direction of ascent.
        newlag: float
            new offset.

        """

        # we first check in which triangle the gradient is
        lag = 0
        vleft = 0
        vright = 0
        kleft = kmax - 1
        kright = kmax + 1
        if kright == len(ldr):
            kright = 0
        if abs(ldr[kmax][0]-ldr[kleft][0])+abs(ldr[kmax][1]-ldr[kleft][1]) == 1:
            if ldr[kleft] in inlist:
                vleft = v[kleft] - z
        if abs(ldr[kmax][0]-ldr[kright][0])+abs(ldr[kmax][1]-ldr[kright][1]) == 1:
            if ldr[kright] in inlist:
                vright = v[kright] - z

#                    print(ldr)
#                    print('vleft, vright', vleft, vright)
        if vleft == vright or max(vleft, vright) <= 0: # no lag, it goes along ldr[kmin]
            lag = 0
        else: # either to the left or the right triangle
            k1 = kmax
            if vleft > vright: # if the left point is higher, gradient is on the left
                k1 = kleft
                v1 = vleft
            if vright > vleft:
                k1 = kright
                v1 = vright
            # here we compute the gradient direction by solving a 2x2 system
            # may be simplified noting that one coefficient is zero, others are ones
            x1 = ldr[kmax][0]
            y1 = ldr[kmax][1]
            x2 = ldr[k1][0]
            y2 = ldr[k1][1]
            z1 = v[kmax] - z
            z2 = v1
            det = 1/(x1*y2 - x2*y1)
            a = det*(y2*z1 - y1*z2)
            b = det*(x1*z2 - x2*z1)
            # the gradient direction is given by (a,b)
            # now we compute the lag
            # the lag is the difference between ldr[kmin] and (a,b)
#                        print('kmax', kmax, 'k1', k1)
#                        print('ldr[kmax]', ldr[kmax], 'ldr[k1]', ldr[k1], 'oldldr', oldldr)
            if (ldr[kmax][0] == ldr[k1][0]): # compute b for a = 1
                b = b/a
                a = ldr[kmax][0]
                if abs(b) < abs(a):
                    lag = b - ldr[kmax][1]
                else:
                    lag = 0
            else:
                a = a/b
                b = ldr[kmax][1]
                if abs(a) < abs(b):
                    lag = a - ldr[kmax][0]
                else: 
                    lag = 0
            # if we keep going the same way, correct with the lag
#                        print('gradient', a, b)
#                        print('oldlag', oldlag, 'lag', lag)
            if (ldr[kmax] == oldldr):
                lag += oldlag
                if lag > 0.5:
                    kmax = k1
                    lag -= 1
                if lag < -0.5:
                    kmax = k1
                    lag += 1
        return ldr[kmax], lag
                   
    def instantiateRidgesAtSaddles(self, divide = False):
        """
        Creates the ridge objects. The procedure takes all the saddles and 
        initiates the ridges by using the segments that were computed when 
        detecting the saddles. Each ridge is defined by one segment starting at
        the saddle. Each ridge is also assigned to a hill. A specific check is 
        done for this if the ridge overlaps a thalweg.
        
        Parameter
        ---------
        divide : Boolean
            Indicates if ridges have to be extracted in view of watershed computation
            If false, only ridges already identified are computed. In case a ridge
            overlaps a thalweg, a side has to be chosen.
            If true, a ridge is computed between two consecutive thalwegs even
            if none was identified at the beginning.            

        Returns
        -------
        None.
        
        """
        
        ridgekey = 1
        divideside = None # record the side for divides overlapping dangling thalwegs
        if divide:
            divideside = {}

        # initialises the ridge lists in hilldict
        for i in self.hilldict:
            self.hilldict[i]['ridge'] = []
            
        # get all the saddles
        # we check z>nodata to avoid starting a ridge in a hole (to be argued)
        saddlelist = [i for i in self.nodedict if self.nodedict[i]['type'] == df.SADDLE and self.nodedict[i]['z']>self.terrain.nodata]
        
        #testipt = 103
        for ipt in saddlelist:
            #if ipt!=testipt:
            #    continue
            pt = self.nodedict[ipt]
                    
            # we get the list of thalwegs
            thalweglist = pt['thalweg']
            #print('Node ', ipt, thalweglist)
            # we need at least two thalwegs to create a ridge
            # if we don't have this, we move to the next node
            if len(thalweglist)<2:
            	# there should not be less than two thalwegs
                print('less than two thalwegs', ipt, thalweglist)
                continue
            (i,j) = pt['ij']
            # orderedneighbours is the 8 neighbours of (i,j) in clockwise order
#            orderedneighbours = [(i+di, j+dj) for (di,dj) in df.ldr8]
            orderedneighbours = [(i+di, j+dj) for (di,dj) in self.terrain.neighbour[i,j]]
            # for each thalweg, we start a ridge on its right
            for it, t in enumerate(thalweglist):
                #print('thalweg', it, t)
                # get the hill right of t (left and right depends on the direction of the thalweg)
                if t > 0:
                    ihill = self.thalwegdict[t]['righthill']
                else:
                    ihill = self.thalwegdict[-t]['lefthill']
                hill = self.hilldict[ihill]

                thalwegline = hill['dangling']
                
                nt = len(thalweglist) -1
                # tright is the next thalweg in the list
                # because the list is in clockwise order, it is the thalweg on the right
                if it == nt:
                    tright = thalweglist[0]
                else:
                    tright = thalweglist[it+1]
                # pleft and pright are coordinates of thalwegs' second points
                if t>0:
                    pleft = self.thalwegdict[t]['polyline'][1]
                else:
                    pleft = self.thalwegdict[-t]['polyline'][-2]
                if tright>0:
                    pright = self.thalwegdict[tright]['polyline'][1]
                else:
                    pright = self.thalwegdict[-tright]['polyline'][-2]
                # ileft and iright are the neighbour pixels where the thalwegs pass through
                ileft = orderedneighbours.index(pleft)
                iright = orderedneighbours.index(pright)
                #print('tleft, tright', t, tright, ihill)
                #print('ileft, iright', ileft, iright)
                # slopeneighbours is the list of pixels between the thalwegs
                slopeneighbours = []
                if ileft <= iright:
                    slopeneighbours = orderedneighbours[ileft:iright+1]
                else:
                    slopeneighbours = orderedneighbours[ileft:]+orderedneighbours[:iright+1]
                # the ridge in this hill must go through one of these pixels
                # if there is already such a ridge in 'line', we take it
                ridgestart = [ir for ir in self.nodedict[ipt]['line'] if ir[2]==1 and ir[1] in slopeneighbours]
                #print(slopeneighbours, ridgestart)
                slopeneighbours = [v for v in slopeneighbours if self.terrain.isInside(v[0],v[1])]
                
                polylist = None
                if not ridgestart: # we did not find a ridge
                    if divide: # we need to define one
                        # take the steepest slope in pixels of slopeneighbours
                        z = self.terrain.dtm[i, j]
                        slope = [2*(self.terrain.dtm[v] - z)/(abs(v[0] - i) + abs(v[1] - j) + 1) for v in slopeneighbours]
                        # get the maximum slope
                        maxslope = max(slope)
                        maxcount = slope.count(maxslope)
                        point = None
                        if maxcount > 1: # if there are several pixels with the max slope
                            maxcandidates = [slopeneighbours[i] for i in range(len(slope)) if slope[i] == maxslope]
                            maxcandidates.sort(key = lambda x: df.order8.index((x[0]-i, x[1]-j)))
                            point = maxcandidates[-1]
                        else:
                            maxpos = slope.index(maxslope)
                            point = slopeneighbours[maxpos]
                        polylist = [(i,j), point]
                    else:
                        continue

                else:                
                    #print('ridgestart', ridgestart[0][0:2])
                    polylist = ridgestart[0][0:2]
                # case where the ridge overlaps a dangling thalweg ending at the saddle

                if divide:
                    # if the ridge overlaps a dangling thalweg, the side of the thalweg must be recorded
                    if pleft == polylist[1] and -t in thalwegline:
                        divideside[ridgekey] = -1
                    if pright == polylist[1] and -tright in thalwegline:
                        divideside[ridgekey] = 1
                else:
                    if -t in thalwegline:
                        if pleft == polylist[1]:
                        # if the ridge overlaps a dangling thalweg, the ridge has already been calculated
                            continue
                    else:
                        if pleft == polylist[1] or (pright == polylist[1] and -tright not in thalwegline):
                            p = polylist[0]
                            q = polylist[1]
                            dp = (q[0] - p[0], q[1] - p[1])
                            d = df.ldr8.index(dp)
                            dz = [0, 0]
                            neighbours =[(0,0)] * 6
                            for k in [1,2,3]:
                                ldr1 = df.ldr8[d - k]
                                ldr2 = df.ldr8[d - 8 + k]
                                neighbours[k-1] = (ldr1, ldr2)
                                ldr1 = (2*ldr1[0], 2*ldr1[1])
                                ldr2 = (2*ldr2[0], 2*ldr2[1])
                                neighbours[k+2] = (ldr1, ldr2)
                            k = 0
                            while dz[0] == dz[1] and k < 6:
    #                            dp1 = df.ldr8[d - k]
    #                            dp2 = df.ldr8[d - 8 + k]
                                (dp1, dp2) = neighbours[k]
                                ldr = [dp1, dp2]
                                z1 = self.terrain.dtm[p[0]+dp1[0], p[1]+dp1[1]]
                                z2 = self.terrain.dtm[p[0]+dp2[0], p[1]+dp2[1]]
                                dz = df.gradLength(pt['z'], [z1, z2], ldr)
                                if k > 2:
                                    print(k, ipt, pleft, pright, ldr, [z1, z2], dz)
                                k += 1
                            if not ((dz[1] > dz[0] and pleft == polylist[1]) or (dz[0] > dz[1] and pright == polylist[1])):
                                #print("on n'y va pas")
                                continue

                #print("New ridge ", ridgekey, divideside, ihill)                
                self.addNewRidge(ridgekey, ihill, ipt, polylist) 
                ridgekey += 1
        if divide:
            for ir in divideside:
                self.ridgedict[ir]['startside'] = divideside[ir]
        
    def traceRidges(self, divide = False, correction = True):
        """
        Computation of ridges starting from saddles. Ridges do not overlap or
        go through saddles. If two ridges join, a junction is inserted and a
        new ridge runs from the junction to the next node. The computation is 
        done in a similar way as for thalwegs and gradients are also constrained 
        to vertices with a correction of the offset. The difference is that a
        ridge must always remain within a polygon formed by its surrounding 
        thalwegs to avoid any thalweg-ridge crossing.

        Parameter
        ---------
        divide : Boolean
            Indicates if ridges have to be extracted in view of watershed computation
            If false, ridges do not stop at confluences. 
            If true, ridges are still computed the same way but intersections with
            confluences are recorded and ridges are split at confluences. 
        correction : boolean
            True if a correction is done to stay closer to the gradient
            False if the descent direction is not corrected as in D8
        Returns
        -------
        None.

        """

        # key for new ridges        
        nextridge = max(self.ridgedict) + 1
        # hill polygons and map of dangling thalwegs
        hillpolygon, danglingmap = self.createHillPolygons()

#        testir = [10363]
        for ihill in self.hilldict:
            floorkey={} # coordinate index recording ridge id
            junctionidx = {} # coordinate index to junctions
            ridgeindangling = {} # coordinate index for ridges on a dangling thalweg
            # for each coordinate, records the ridge id and the side of the thalweg to check
            # if another rige comes there, a junction should be inserted
            
            if ihill % 1000 == 0:
                print(ihill)
            hill = self.hilldict[ihill]
            ridgelist = hill['ridge'][:]
            for ir in ridgelist:
#                if ir not in testir:
#                    continue
                r = self.ridgedict[ir]
                polylist = r['polyline']
                ipt = r['start']
                #pt = self.nodedict[ir]
                # computes a hill starting from t at ipt
                #thalwegring, thalwegline = self.getThalwegRing(ipt, t)
                
                # in regular cases, thalwegline is empty and thalwegring contains only one ring, the boundary,
                # and the ring goes clockwise
                
                # if there are several rings, the biggest one is the boundary, others are holes
                
                # singular case: a node has two thalwegs that yield the same ring in different directions
                # this case can occur for two thalwegs from a pass connecting the same pit
                
                # create the shapely polygon from the hill
                plgn = hillpolygon[ihill]
                
                (i,j) = polylist[1]
                dtoverlap = danglingmap[i,j] # thalweg id on which the ridge is
                thalwegside = 0 # -1/0/1 if left/undefined/right side
                # on which side of the dangling line we come from
                if dtoverlap > 0:
                    if divide and 'startside' in self.ridgedict[ir]: 
                        thalwegside = self.ridgedict[ir]['startside']
                    else:
                        tmppoly = self.thalwegdict[dtoverlap]['polyline']
                        kt = tmppoly.index((i,j))
                        if kt == 0:
                            print('kt = 0, thalweg', dtoverlap, 'ridge', ir)
                            thalwegside = 0
                        else:
                            thalwegside = df.sideOfLine(tmppoly, kt, polylist[0])
                danglingconfluence = False
                confluencewedge = ()

                keepgoing = True # false once a critical point is reached
                #outbound = False # true if out of bound (obsolete)
                junction = False # true if a new junction point to be created
                okjoin = True # false if the ridge reaches another ridge on the other side of the thalweg
                endindex = -1
                
                nextrounds = False # only for the first round
                validridge = True # false if an exception occured in ascent computation
                
                oldlag = 0
                oldldr = (0,0)
                while keepgoing:
                    if nextrounds:
                        ldr = self.terrain.neighbour[i,j]
                        v = self.terrain.getNeighbourHeights(i, j)
                        z = self.terrain.dtm[i,j]
    
                        [validridge, kmax, inlist] = self.computeAscent(i, j, v, z, ldr, plgn, danglingconfluence, confluencewedge, thalwegside, dtoverlap, danglingmap)

                        if not validridge:
                            print('validridge', ir, ihill)
                        oldij = (i,j) # used later to check sides
                        
                        # the steepest gradient is for v[kmax]
                        # we check if we need to modify it according to the lag
                        if correction:
                            [oldldr, oldlag] = self.correctAscent(v, z, kmax, ldr, inlist, oldldr, oldlag) 
                        else:
                            oldldr = ldr[kmax]
                        i += oldldr[0]
                        j += oldldr[1]
                        # if (i,j) in polylist:
                        #     print(i, j, "we're turning back")
                        #     print(ipt, ir, ihill, polylist)
                        #     polylist = polylist[:polylist.index((i,j))]
                        #     validridge = False
                        #     break
                        polylist.append((i,j))

                    nextrounds = True                    
                    # check dangling thalwegs
                    okjoin = True
                    danglingconfluence = False
                    oldoverlap = dtoverlap
                    dtoverlap = danglingmap[i,j]
                    # no dangling thalweg under this pixel
                    if dtoverlap <= 0:
                        # if we reach a confluence
                        if self.terrain.pixelclass[i,j] == df.CONFLUENCE:
                            iconf = self.nodeidx[(i,j)]
                            confthalweg = self.nodedict[iconf]['thalweg']
                            # check if first node of dangling thalweg
                            lthalweg = [-t for t in confthalweg if -t in hill['dangling']]
                            if lthalweg:
                                oldtlwg = [t for t in confthalweg if t > 0][0]
                                if len(lthalweg) > 1:
                                    print("trouble ahead: more than one dangling thalwegs", lthalweg)
                                dtoverlap = lthalweg[0]
                                otherthalweg = [-t for t in confthalweg if t < 0 and -t != dtoverlap][0]
                                p3 = self.thalwegdict[oldtlwg]['polyline'][1]
                                p2 = self.thalwegdict[dtoverlap]['polyline'][-2]
                                p1 = self.thalwegdict[otherthalweg]['polyline'][-2]
                                confluencewedge = (p3, (i,j), p2)
                                danglingconfluence = True
                                thalwegside = -df.sideOfWedge(confluencewedge, p1)
                                #print(ir, confluencewedge, p1, thalwegside)
                        else:
                            thalwegside = 0 # we are no longer on a thalweg
                    # if we are on a dangling thalweg
                    else:
                        # if it is a confluence (that's the last pixel of the thalweg)
                        if self.terrain.pixelclass[i,j] == df.CONFLUENCE:
                            iconf = self.nodeidx[(i,j)]
                            lthalweg = self.nodedict[iconf]['thalweg']
                            ithalweg = lthalweg.index(dtoverlap)
                            it1 = df.getSecondPointPosition(lthalweg[ithalweg])
                            p1 = self.thalwegdict[dtoverlap]['polyline'][it1]
                            # if we don't know on which side to stay
                            if thalwegside == 0:
                                jthalweg = ithalweg - 1
                                kthalweg = ithalweg + 1
                                if kthalweg == len(lthalweg):
                                    kthalweg = 0
                                it2 = df.getSecondPointPosition(lthalweg[jthalweg])
                                it3 = df.getSecondPointPosition(lthalweg[kthalweg])
                                p2 = self.thalwegdict[abs(lthalweg[jthalweg])]['polyline'][it2]
                                p3 = self.thalwegdict[abs(lthalweg[kthalweg])]['polyline'][it3]
                                confluencewedge = (p3, (i,j), p2)
                                thalwegside = df.sideOfWedge(confluencewedge, p1)
                            # if we know on which side to stay
                            else:
                                jthalweg = ithalweg + thalwegside
                                if jthalweg == len(lthalweg):
                                    jthalweg = 0
                                #print(ipt, 'reached confluence', iconf, lthalweg[jthalweg], thalwegside)
                                it2 = df.getSecondPointPosition(lthalweg[jthalweg])
                                p2 = self.thalwegdict[abs(lthalweg[jthalweg])]['polyline'][it2]
                                confluencewedge = (p1,(i,j),p2)
                            danglingconfluence = True
                        # if we are not at a confluence (regular thalweg pixel)
                        else:
                            if oldoverlap > 0 and dtoverlap != oldoverlap and self.terrain.pixelclass[i,j] != df.SADDLE:
                                thalwegside = 0
                            # if we don't know on which side to stay, find it out
                            if thalwegside == 0:
                                #print('xx', dtoverlap, (i,j))
                                tmppoly = self.thalwegdict[dtoverlap]['polyline']
                                kt = tmppoly.index((i,j))
                                if kt != 0:
                                    thalwegside = df.sideOfLine(tmppoly, kt, oldij)
                                #print(ir, 'thalwegside = ', thalwegside)
                            # if there is already a ridge here
                            if (i,j) in floorkey:
                                #if len(ridgeindangling[(i,j)]) > 1:
                                    #print("Check this (more than one ridge)", ipt, ir, i, j, ridgeindangling[(i,j)])
                                jnctn = [rd[1] for rd in ridgeindangling[(i,j)]]
                                if not thalwegside in jnctn:
                                    okjoin = False
                                #elif -thalwegside in jnctn:
                                #    print(i,j, ir, thalwegside, ridgeindangling[(i,j)])
                                ridgeindangling[(i,j)].append((ir, thalwegside))
                            else:
                                ridgeindangling[(i,j)] = [(ir, thalwegside)]
                            
                    # check if the point is a node or out (it should not happen)
                    if not self.terrain.isInside(i,j):
                        keepgoing = False
                        #outbound = True
                        print("we're out")
                        
                    # check what kind of pixel we reach and whether we keep going
                    if self.terrain.pixelclass[i,j] == df.SLOPE:
                        # nothing special, keep going
                        keepgoing = True
                        self.terrain.pixelclass[i,j] = df.RIDGE
                        floorkey[i,j] = ir
                    # if divide is false (default) confluences are not considered 
                    # as end points since ridges only start from saddles
                    # if divide is true, ridges can start and end at confluences
                    elif self.terrain.pixelclass[i,j] == df.THALWEG or (self.terrain.pixelclass[i,j] == df.CONFLUENCE and not divide):
                        if (i,j) in floorkey: # we reach a ridge
                            if not okjoin:
                                junction = False
                                keepgoing = True
                                floorkey[i,j] = ir
                            else:
                                if (i,j) not in junctionidx: # we need to create a junction
                                    junction = True
                                keepgoing = False
                        else:
                            keepgoing = True
                            floorkey[i,j] = ir
                    elif divide and self.terrain.pixelclass[i,j] == df.CONFLUENCE:
                        if (i,j) in floorkey:
                            junction = False
                            keepgoing = False
                            if dtoverlap > 0:
                                ridgeindangling[(i,j)].append((ir, thalwegside))
                        else:
                            floorkey[i,j] = ir
                            if dtoverlap > 0:
                                ridgeindangling[(i,j)] = [(ir, thalwegside)]
                            keepgoing = False
                            junction = False
                    elif self.terrain.pixelclass[i,j] == df.RIDGE:
                        keepgoing = False
                        junction = True
                    else: # ridge terminates at a saddle or peak
                        keepgoing = False
                        #outbound = False
                # end of while loop reached when keepgoing is false
            
            
                if validridge:
#                    idend = -1
#                    if self.terrain.pixelclass[i,j] == CONFLUENCE:
#                        print('Ridge reached a confluence', ipt, ridgekey, i,j)
                    if junction:                    
                        # if we reach an existing ridge, we have to create a junction
                        idridge = floorkey[i,j]
                        if (i,j) in ridgeindangling:
                            idridge = [rd[0] for rd in ridgeindangling[(i,j)] if rd[1] == thalwegside][0]
                            #print('idridge = ', idridge, ir, ridgeindangling[(i,j)])
                        # we create the new node
                        endindex = self.nodekey
                        self.nodeidx[(i,j)] = self.nodekey
                        junctionidx[(i,j)] = self.nodekey
                        self.nodedict[self.nodekey] = {
                            'ij' : (i,j),
                            'z': self.terrain.dtm[i,j],
                            'thalweg' : [],
#                            'ridge' : [-idridge, nextridge],
                            'ridge' : [-idridge],
                            'type' : df.JUNCTION
                        }
                        #print(endindex, ir, idridge, nextridge, self.nodedict[endindex])
                        if self.terrain.pixelclass[i,j] == df.RIDGE:
                            self.terrain.pixelclass[i,j] = df.JUNCTION
                        self.nodekey += 1
                        # split the existing ridge in two
                        r = self.ridgedict[idridge]
                        ptlist = r['polyline']
                        iend = r['end']
                        try:
                            icfl = ptlist.index((i,j)) # point where we split
                        except ValueError:
                            print('in createRidges ptlist.index', ipt, ir, idridge, ihill)
                            print('polylist', polylist)
                            print('ij', i, j, r)
                            break
                        ptlist1 = ptlist[:icfl+1]
                        ptlist2 = ptlist[icfl:]
                        # modify the first ridge to stop at the junction
                        r['polyline'] = ptlist1
                        r['end'] = endindex
                        # create a new ridge
                        self.addNewRidge(nextridge, ihill, endindex, ptlist2)
                        #self.nodedict[endindex]['ridge'].append(nextridge)
                        self.ridgedict[nextridge]['end'] = iend
                        #if thalwegside:
                            #self.ridgedict[nextridge]['endside'] = thalwegside
                        if 'endside' in r:
                            self.ridgedict[nextridge]['endside'] = r['endside']
                        
                        # we update the node at the end of the new ridge
                        ndt = self.nodedict[iend]['ridge']
                        ndt.remove(-idridge)
                        ndt.append(-nextridge)
                        # we update floorkey with nextridge
                        
                        for k in ptlist2[1:-1]:
                            floorkey[k] = nextridge
                            if k in ridgeindangling:
                                pos = [i for i, x in enumerate(ridgeindangling[k]) if x[0] == idridge][0]
                                tmp = ridgeindangling[k][pos][1]
                                del ridgeindangling[k][pos]
                                ridgeindangling[k].append((nextridge, tmp))
                        nextridge += 1
                    else:
                        endindex = self.nodeidx[(i,j)]
                        #if self.terrain.pixelclass[i,j] == df.CONFLUENCE:
                        #    print('Ridge reached a confluence', ipt, ir, i,j)

                    # the last point is a critical point
                    # we add its index to the ridgelist
                    self.nodedict[endindex]['ridge'].append(-ir)
                    # self.ridgedict[ir]={
                    #     'start' : ipt,
                    #     'polyline' : polylist,
                    #     'end' : endindex,
                    #     'hill' : ihill
                    # }
                    self.ridgedict[ir]['polyline'] = polylist
                    self.ridgedict[ir]['end'] = endindex
                    if thalwegside:
                        self.ridgedict[ir]['endside'] = thalwegside


    def orderRidgesAroundNodes(self):
        """
        Orders the list of ridges around each node in an anticlockwise order 
        
        To be done: holes are not considered.
    
        Returns
        -------
        None.
    
        """
        for ipt, pt in self.nodedict.items():
            # -1 is the id of the virtual pit closing outside the domain
            if ipt == -1:
                continue
            ridgelist = pt['ridge']
            nr = len(ridgelist)
            if nr < 2:
                continue
            (i,j) = pt['ij']
            neighbours = [(i+di, j+dj) for (di,dj) in df.ldr8]
            # arrange the lines in the same order as the neighbours
            try:
                ridgelist.sort(key = lambda it: neighbours.index(self.ridgedict[abs(it)]['polyline'][df.getSecondPointPosition(it)]))
            except (IndexError, ValueError) as err:
                print('Exception raised in orderRidgesAroundNodes', err, ipt, pt, neighbours)
            # check if two ridges overlap at the node
            # this only occurs for two consecutive ridges ending at a saddle
            for i in range(nr):
                ir1 = ridgelist[i-1]
                ir2 = ridgelist[i]
                hill1 = self.ridgedict[abs(ir1)]['hill']
                hill2 = self.ridgedict[abs(ir2)]['hill']
                if hill1 == hill2:
                    print("in orderRidgeAroundNodes, check", ipt, ir1, ir2)
                    continue
                if pt['type'] == df.SADDLE and ir1 < 0 and ir2 < 0:
                    secondpoint = self.ridgedict[-ir1]['polyline'][-2]
                    if secondpoint == self.ridgedict[-ir2]['polyline'][-2]:
                        thalweglist = pt['thalweg']
                        reverse = False
                        for it in thalweglist:
                            #if it < 0:
                                #print(ir1, ir2, it)
                            if it < 0:
                                t = self.thalwegdict[-it]
                                lefthill = t['righthill']
                                righthill = t['lefthill']
                            else:
                                t = self.thalwegdict[it]
                                lefthill = t['lefthill']
                                righthill = t['righthill']
                            if hill2 == lefthill and hill1 == righthill and t['polyline'][1] == secondpoint:
                                reverse = True
                                break
                        if reverse:
                            ridgelist[i-1] = ir2
                            ridgelist[i] = ir1

    def getRidgeRing(self, idnode, idridge):
        """
        For a given node and a given ridge starting at this node, the function
        returns the list of ridges that delineate a dale when walking from this
        ridge in a clockwise order. The function takes into account overlapping
        segments although it has not been checked thoroughly (if there are overlaps
        inside a ring).

        Parameters
        ----------
        idnode : integer
            Node at which the hill starts.
        idridge : integer
            First ridge in the list, index can be negative depending on the
            ridge direction.

        Returns
        -------
        ridgering: a list of lists of oriented ridges forming rings, the first list
        is the outer boundary, other lists are holes
        ridgeline: a list of ridges representing dangling ridges inside the 
        dale. They do not form a ring.

        """
        # case of a line and no polygon has not been tested yet
        idridge = abs(idridge)
        currentnode = idnode
        currentridge = idridge
        keepgoing = True
        ridgering = []
        while keepgoing:
            r = self.ridgedict[currentridge]
            if currentnode == r['start']:
                nextnode = r['end']
                ridgering.append(currentridge)
            elif currentnode == r['end']:
                nextnode = r['start']
                ridgering.append(-currentridge)
            else:
                print('error in getRidgeRing: ridge', currentridge, 'and node', currentnode, 'do not match')

            currentnode = nextnode
            ridgelist = [abs(r) for r in self.nodedict[currentnode]['ridge']]
            ridgepos = ridgelist.index(currentridge)
            currentridge = ridgelist[ridgepos - 1]
            keepgoing = not ((currentnode == idnode) and (currentridge == idridge))
        # if we have overlapping segments
        # we remove these overlapping segments to split the ring in
        # the outer and inner rings and inner lines
        
        nt = len(ridgering)
        tbr = []
        innerline = []
        innerring = []
        
        it1=0
        while it1 < nt:
            t1 = ridgering[it1]
            it2 = nt-1
    #           print 't1 = ', t1
            while it2 > it1:
                t2 = ridgering[it2]
    #               print '   t2 = ', t2
                if (t2 == -t1):
    #                   print 't1=-t2=', t1
                    while (t2 == -t1) and (it2 > it1):
                        innerline.append(abs(t2))
    #                       print 't2 =', t2, 'innerline = ', innerline
                        tbr.append(t1)
                        tbr.append(t2)
                        it1 += 1
                        it2 -= 1
                        t1 = ridgering[it1]
                        t2 = ridgering[it2]
                    if (t2 != -t1):
                        innerring.append(ridgering[it1:it2+1])
                        tbr.extend(innerring[-1])
                    it1 -= 1
                    break
    # bug: we should go iteratively in the innerring to check if any other overlapping thalweg
                it2 -= 1
            it1 += 1

        ridgeboundary = [t for t in ridgering if t not in tbr]
        if ridgeboundary == []:
            ridgering = innerring
        else:
            ridgering = [ridgeboundary] + innerring
        for ir, ring in enumerate(ridgering): # this corrects a bug but to be checked if we can do better
            for ir2, ring2 in enumerate(ridgering[ir+1:]):
                for t in ring2:
                    if t in ring:
                        ring.remove(t)
        for t in innerline: # this may correct the bug but is a fix, to be checked
            for ring in ridgering:
                if t in ring:
                    ring.remove(t)
                if -t in ring:
                    ring.remove(-t)
        while [] in ridgering:
            ridgering.remove([])
        return ridgering, innerline

    def computeDales(self):
        """
        Computes the dictionary of dales of the surface network. A dale is an
        ascending manifold. It is formed by a series of ridges. The dale can
        contain some holes.
        Each dale is defined by a series of ridges forming the boundary and
        a series of holes formed each by a list of ridges. Ridges are oriented
        in the lists. The dale also contains a list of dangling ridges, ridges
        that do not form a loop.
        The function also adds a pointer to the left and right dale to ridges
        in ridgedict. For a dangling ridge, the left and right dales are the
        same.
    
        Returns
        -------
        None.
    
        """
        # initialise the link from thalweg to hill
        for ir in self.ridgedict:
            r = self.ridgedict[ir]
            r['leftdale'] = None
            r['rightdale'] = None
    
        self.daledict = {}
        dalekey = 0
        saddlelist = [i for i in self.nodedict if self.nodedict[i]['type'] == df.SADDLE and self.nodedict[i]['z']>self.terrain.nodata]
    
        for ipt in saddlelist:
            pt = self.nodedict[ipt]
                    
            # we get the list of thalwegs
            ridgelist = pt['ridge']
            # we need at least two thalwegs to create a hill
            # if we don't have this, there's a problem in the data
            # a saddle always has two thalwegs at least
            if len(ridgelist)<2:
                print('less than two ridges', ipt, ridgelist)
                continue
            # for each thalweg, we build a new hill on its right
            for r in ridgelist:
                # computes the ring and dangling thalwegs starting from t at ipt
                absr = abs(r)
                if self.ridgedict[absr]['rightdale']:
                    continue
                ridgering, ridgeline = self.getRidgeRing(ipt, r)
                # if t is a danglilng thalweg, no new hill is created because there's none
                if absr in ridgeline:
                    continue
                # ridgering can include several rings, one of them is the boundary
                # the others are holes
                # we compare the bounding boxes to check if the first ring is the boundary
                arealist = []
                for ring in ridgering:
                    pointlist = set(reduce(lambda a,b:a+b, [self.ridgedict[abs(rr)]['polyline'] for rr in ring]))
                    xmin = min(pointlist, key = lambda a: a[0])[0]
                    xmax = max(pointlist, key = lambda a: a[0])[0]
                    ymin = min(pointlist, key = lambda a: a[1])[1]
                    ymax = max(pointlist, key = lambda a: a[1])[1]
                    area = (xmax - xmin) * (ymax - ymin)
                    arealist.append(area)
                border = ridgering.pop(0)
                holes = ridgering
                # if the first ring is not the biggest, the dale is a hole in a larger one
                # we discard everything outside the hole
                if max(arealist) != arealist[0]:
                    holes = []
                    ridgeline = []
                
                self.daledict[dalekey] = {'boundary': border, 
                                          'hole': holes,
                                          'dangling': ridgeline}
                for rr in border:
                    if rr > 0:
                        self.ridgedict[rr]['rightdale'] = dalekey
                    else:
                        self.ridgedict[-rr]['leftdale'] = dalekey
                for r in holes:
                    for rr in r:
                        if rr > 0:
                            self.ridgedict[rr]['rightdale'] = dalekey
                        else:
                            self.ridgedict[-rr]['leftdale'] = dalekey
                for rr in ridgeline:
                    self.ridgedict[rr]['leftdale'] = dalekey
                    self.ridgedict[rr]['rightdale'] = dalekey
                dalekey += 1
    
        for i in self.thalwegdict:
            self.thalwegdict[i]['dale'] = -1
            
        # correction of a bug from dale computation where some dales were computed twice
        nopit = set()
        for idale in self.daledict:
            tlist = self.daledict[idale]['boundary']
            for it in tlist:
                t = self.ridgedict[abs(it)]
                if t['leftdale'] != idale and t['rightdale'] != idale:
                    nopit.add(idale)
                    break
        for i in nopit:
            del self.daledict[i]

        
        # now, we assign a dale to each thalweg
        # that's a bit messy here and may be simplified
        thalwegset = set() 
        danglinglist = []
        noridgeneighbour = []
        for ipt in saddlelist:
            pt = self.nodedict[ipt]
            tlist = pt['thalweg']
            for i, it in enumerate(tlist):
                if it < 0:
                    continue
                idale = -1
                tlwg = self.thalwegdict[it]
                ilefthill = tlwg['lefthill']
                leftlr = [ir for ir in pt['ridge'] if self.ridgedict[abs(ir)]['hill'] == ilefthill]
                irighthill = tlwg['righthill']
                rightlr = [ir for ir in pt['ridge'] if self.ridgedict[abs(ir)]['hill'] == irighthill]
                # if the right and left ridges are the same, we cannot choose a dale
                if leftlr == rightlr:
                    # this case occurs if there is no ridge around
                    if leftlr == []:
                        noridgeneighbour.append(it)
                    else:
                    # or if we have a dangling thalweg
                        danglinglist.append(it)
                    continue
                # if there's only one ridge on the left, easy
                if len(leftlr) == 1:
                    ir = leftlr[0]
                    if ir > 0:
                        idale = self.ridgedict[ir]['rightdale']
                    else:
                        idale = self.ridgedict[-ir]['leftdale']
                    tlwg['dale'] = idale
                    iend = tlwg['end']
                    # we note down if there is a confluence to set next thalwegs later
                    if self.nodedict[iend]['type'] == df.CONFLUENCE:
                        tlwg = self.nodedict[iend]['thalweg']
                        nextthalweg = [j for j in tlwg if j > 0]
                        thalwegset.add((nextthalweg[0], idale))
                    continue
                
                # if there's only one ridge on the right, same easy
                if len(rightlr) == 1:
                    ir = rightlr[0]
                    if ir > 0:
                        idale = self.ridgedict[ir]['leftdale']
                    else:
                        idale = self.ridgedict[-ir]['rightdale']
                    tlwg['dale'] = idale
                    iend = tlwg['end']
                    if self.nodedict[iend]['type'] == df.CONFLUENCE:
                        tlwg = self.nodedict[iend]['thalweg']
                        nextthalweg = [j for j in tlwg if j > 0]
                        thalwegset.add((nextthalweg[0], idale))
                    continue
                
                # if there's no side with just one ridge,
                anglemin = -1
                anglet = df.polylineangle(tlwg['polyline'], tlwg['polyline'][0])
                angler = []
                # if there are several ridges on the left
                if len(leftlr) > 0:
                    # we sort them according to the angle
                    lr = leftlr
                    for ir in lr:
                        if ir > 0:
                            poly = self.ridgedict[ir]['polyline']
                            first = poly[0]
                        else:
                            poly = self.ridgedict[-ir]['polyline']
                            first = poly[-1]
                        angle = df.polylineangle(poly, first) - anglet
                        if angle < 0:
                            angle += 6.283185308
                        angler.append(angle)
                    anglemin = min(angler)
                    i = angler.index(anglemin)
                    ir = lr[i]
    #                print(angler, ir)
    
                    if ir > 0:
                        idale = self.ridgedict[ir]['rightdale']
                    else:
                        idale = self.ridgedict[-ir]['leftdale']
                # if there are several ridges on the right
                elif len(rightlr) > 0:
                    lr = rightlr
                    for ir in lr:
                        if ir > 0:
                            poly = self.ridgedict[ir]['polyline']
                            first = poly[0]
                        else:
                            poly = self.ridgedict[-ir]['polyline']
                            first = poly[-1]
                        angle = anglet - df.polylineangle(poly, first)
                        if angle < 0:
                            angle += 6.283185308
                        angler.append(angle)
                    if not angler:
                        print(it, rightlr, angle)
                    anglemin = min(angler)
                    i = angler.index(anglemin)
                    ir = lr[i]
    #                print(angler, ir)
    
                    if ir > 0:
                        idale = self.ridgedict[ir]['leftdale']
                    else:
                        idale = self.ridgedict[-ir]['rightdale']
                else:
                    print("no dale", it)
                    
                # but if the thalweg overlaps a ridge and we cannot choose a side
                # we'll see this later
                if irighthill == ilefthill and anglemin == 0:
                    idale = -1
                if idale > -1:
                    tlwg['dale'] = idale
                    iend = tlwg['end']
                    if self.nodedict[iend]['type'] == df.CONFLUENCE:
                        thalweg = self.nodedict[iend]['thalweg']
                        nextthalweg = [j for j in thalweg if j > 0]
                        thalwegset.add((nextthalweg[0], idale))
    
        # we complete all thalwegs located after a confluence for those we know    
        while thalwegset:
            itd = thalwegset.pop()
            tlwg = self.thalwegdict[itd[0]]
            tlwg['dale'] = itd[1]
            ipt =  tlwg['end']
            pt = self.nodedict[ipt]
            if pt['type'] == df.CONFLUENCE:
                tlist = pt['thalweg']
                nextthalweg = [j for j in tlist if j > 0]
                thalwegset.add((nextthalweg[0], itd[1]))
        
        # we look at dangling lines that do not overlap a ridge
        # we also choose a dale based on the closest ridge
        danglinglistbis = []
        for it in danglinglist:
            ipt = self.thalwegdict[it]['start']
            pt = self.nodedict[ipt]
            
            poly = self.thalwegdict[it]['polyline']
            first = poly[0]
            tangle = df.polylineangle(poly, first)
            
            rlist = pt['ridge']
            rangle = []
            for ir in rlist:
                if ir > 0:
                    poly = self.ridgedict[ir]['polyline']
                    first = poly[0]
                else:
                    poly = self.ridgedict[-ir]['polyline']
                    first = poly[-1]
                angle = df.polylineangle(poly, first) - tangle
                if angle < 0:
                    angle += 6.283185308
                rangle.append(angle)
            minangle = min(rangle)
            if minangle > 0:
                ir = rlist[rangle.index(minangle)]
                if ir > 0:
                    idale = self.ridgedict[ir]['rightdale']
                else:
                    idale = self.ridgedict[-ir]['leftdale']
                self.thalwegdict[it]['dale'] = idale
            else:
                danglinglistbis.append(it)
        danglinglist = danglinglistbis        
        
        # we do the same for remaining dangling thalwegs but starting from the end
        # if it is a saddle
        danglinglistbis = []
        for it in danglinglist:
            ipt = self.thalwegdict[it]['end']
            pt = self.nodedict[ipt]
            if pt['type'] != df.SADDLE:
                continue
            poly = self.thalwegdict[it]['polyline']
            first = poly[-1]
            tangle = df.polylineangle(poly, first)
            
            rlist = pt['ridge']
            rangle = []
            for ir in rlist:
                if ir > 0:
                    poly = self.ridgedict[ir]['polyline']
                    first = poly[0]
                else:
                    poly = self.ridgedict[-ir]['polyline']
                    first = poly[-1]
                angle = df.polylineangle(poly, first) - tangle
                if angle < 0:
                    angle += 6.283185308
                rangle.append(angle)
            if rangle:
                minangle = min(rangle)
            else:
                minangle = 0
            if minangle > 0:
                ir = rlist[rangle.index(minangle)]
                if ir > 0:
                    idale = self.ridgedict[ir]['rightdale']
                else:
                    idale = self.ridgedict[-ir]['leftdale']
                self.thalwegdict[it]['dale'] = idale
            else:
                danglinglistbis.append(it)
                
        # we look at pits and confluences if any thalweg without a dale
        # this is more of a check and may be done later
        pitconflist = [i for i in self.nodedict if self.nodedict[i]['type'] == df.PIT or self.nodedict[i]['type'] == df.CONFLUENCE]
        for i in pitconflist:
            tlist = self.nodedict[i]['thalweg']
            dale = [self.thalwegdict[abs(j)]['dale'] for j in tlist]
            unknown = [k for k in dale if k == -1]
            known = {k for k in dale if k != -1}
            if not known:
                #print('no dale in pit', i)
                continue
            if len(known) > 1:
                print('error pit', i)
            else:
                if unknown:
                    idale = known.pop()
                    for j in tlist:
                        if self.thalwegdict[abs(j)]['dale'] == -1:
                            self.thalwegdict[abs(j)]['dale'] = idale
                            
        # for thalwegs which had no ridge around,
        # we set the same dale as the thalweg next to them
        for it in noridgeneighbour:
            ipt = self.thalwegdict[it]['start']
            tlist = self.nodedict[ipt]['thalweg']
            i = tlist.index(it)
            jt = tlist[i - 1]
            self.thalwegdict[it]['dale'] = self.thalwegdict[abs(jt)]['dale']
            
        # last case are dangling thalwegs overlapping a ridge
        # only those ending at a saddle may remain, others have been processed
        # when looking at junctions and pits
        for it in danglinglistbis:
            t = self.thalwegdict[it]
            ipt = t['start']
            pt = self.nodedict[ipt]
            rlist = [-k for k in pt['ridge'] if k < 0]
            thalwegpolyline = t['polyline']
            rsublist = [k for k in rlist if self.ridgedict[k]['polyline'][-2] == thalwegpolyline[1]]
            ir = rsublist[0]
            r = self.ridgedict[ir]
            ridgepolyline = r['polyline']
            # first case, if the two lines are not equal, we can check on which side of the ridge
            # is the thalweg point outside the ridge
            nt = len(thalwegpolyline)
            nr = -len(ridgepolyline) -1
            ri = -2
            ti = 1
            b = True
            c = True
            while b:
                c = thalwegpolyline[ti] == ridgepolyline[ri]
                ti += 1
                ri -= 1
                b = (ti < nt) and (ri > nr) and c
            if not c:
                r1 = ridgepolyline[ri+2]
                r2 = ridgepolyline[ri+1]
                t2 = thalwegpolyline[ti-1]
                r12 = (r2[0] - r1[0], r2[1] - r1[1])
                t12 = (t2[0] - r1[0], t2[1] - r1[1])
                pv = r12[0]*t12[1] - r12[1]*t12[0]
                # r12 goes backward so sign of vector product is inversed for ridge
                if pv > 0:
                    t['dale'] = r['rightdale']
                else:
                    t['dale'] = r['leftdale']
            else:
            # and if both lines are equal, we check the gradient
            # the thalweg is on the steeper side downward
                p = thalwegpolyline[0]
                q = thalwegpolyline[1]
                dp = (q[0] - p[0], q[1] - p[1])
                d = df.ldr8.index(dp)
                dz = [0, 0]
                neighbours =[(0,0)] * 6
                for k in [1,2,3]:
                    ldr1 = df.ldr8[d - k]
                    ldr2 = df.ldr8[d - 8 + k]
                    neighbours[k-1] = (ldr1, ldr2)
                    ldr1 = (2*ldr1[0], 2*ldr1[1])
                    ldr2 = (2*ldr2[0], 2*ldr2[1])
                    neighbours[k+2] = (ldr1, ldr2)
                k = 0
                while dz[0] == dz[1] and k < 6:
                    (dp1, dp2) = neighbours[k]
                    ldr = [dp1, dp2]
                    z1 = self.terrain.dtm[p[0]+dp1[0], p[1]+dp1[1]]
                    z2 = self.terrain.dtm[p[0]+dp2[0], p[1]+dp2[1]]
                    dz = df.gradLength(pt['z'], [z1, z2], ldr)
                    if k > 2:
                        print(k, ipt, ldr, [z1, z2], dz)
                    k += 1
                if dz[0] < dz[1]:
                    t['dale'] = r['rightdale']
                else:
                    t['dale'] = r['leftdale']
                    
    
        # this is a check to see if any mistake or missing dale    
        danglingunknown = []
        danglingknown = []
        for it in self.thalwegdict:
            t = self.thalwegdict[it]
            ihill = t['lefthill']
            if it in self.hilldict[ihill]['dangling']:
                if t['dale'] == -1:
                    danglingunknown.append(it)
                else:
                    danglingknown.append(it)
        print(danglingunknown)

    def removeSpuriousRidges(self):
        toberemoved = set()
        for inode, node in self.nodedict.items():
            ridgelist = node['ridge']
            for i in range(len(ridgelist)):
                if ridgelist[i]*ridgelist[i-1] < 0:
                    ir1 = abs(ridgelist[i])
                    ir2 = abs(ridgelist[i-1])
                    r1 = self.ridgedict[ir1]
                    r2 = self.ridgedict[ir2]
                    if r1['start'] == r2['end'] and r1['end'] == r2['start'] and r1['hill'] == r2['hill']:
                        #print('node', inode, 'ridges', ir1, ir2)
                        toberemoved.add(max(ir1, ir2))
        print("Spurious ridges", toberemoved, "are removed")
        for ir in toberemoved:
            r = self.ridgedict[ir]
            ihill = r['hill']
            istart = r['start']
            iend = r['end']
            self.nodedict[istart]['ridge'].remove(ir)
            self.nodedict[iend]['ridge'].remove(-ir)
            self.hilldict[ihill]['ridge'].remove(ir)
            del self.ridgedict[ir]
            
def checkCells(self):
    for idale in self.daledict:
        self.daledict[idale]['pit'] = None
    
    # for each pit, assign the pit to the dale
    for ipit in self.nodedict:
        pit = self.nodedict[ipit]
        if pit['type'] != df.PIT:
            continue
        tlist = pit['thalweg']
        if not tlist:
            print("pit ", ipit, "is isolated")
        it = -tlist[0]
        idale = self.thalwegdict[it]['dale']
        for i in tlist:
            if self.thalwegdict[-i]['dale'] != idale:
                print("pit ", ipit, "problem with dale")
        if self.daledict[idale]['pit']:
            print("dale", idale, "already pit ", ipit, self.daledict[idale]['pit'])
        self.daledict[idale]['pit'] = ipit
    
    nopit = []
    for idale in self.daledict:
        if not self.daledict[idale]['pit']:
            nopit.append(idale)

    for ihill in self.hilldict:
        self.hilldict[ihill]['peak'] = None
    
    # for each pit, assign the pit to the dale
    for ipeak in self.nodedict:
        peak = self.nodedict[ipeak]
        if peak['type'] != df.PEAK:
            continue
        tlist = peak['ridge']
        if not tlist:
            print("peak ", ipeak, "is isolated")
        it = -tlist[0]
        ihill = self.ridgedict[it]['hill']
        for i in tlist:
            if self.ridgedict[-i]['hill'] != ihill:
                print("peak ", ipeak, "problem with hill")
        if self.hilldict[ihill]['peak']:
            print("hill", ihill, "already peak ", ipeak)
        self.hilldict[ihill]['peak'] = ipeak
    
    nopeak = []
    for ihill in self.hilldict:
        if not self.hilldict[ihill]['peak']:
            nopeak.append(ihill)
    
    return nopit, nopeak
