# coding=utf-8
# coding=utf-8
"""
thalwegnetwork.py:
    Class defining a thalweg network.
    
    author: Eric Guilbert
"""

from time import perf_counter_ns
from math import sqrt
from functools import reduce

import pickle
import numpy as np
import shapely.geometry

import deffunction as df

class ThalwegNetwork:
    # a class describing a thalweg network composed of
    # nodedict: a dictionary of critical points in the network (pit, saddle, confluence)
    # thalwegdict: a dictionary of thalwegs connecting at critical points
    # hilldict: a dictionary of hills (ascending cells)
    
    # saddles are detected first
    # terrain triangulation is computed
    # thalwegs are computed and hills are defined from thalwegs

    # constructor building an empty network
    def __init__(self):
        """
        Constructor creating dictionaries for nodes, ridges and thalwegs.
        All kinds of nodes are stored in the same dictionary.
        nodeidx is a hash table indexing nodes on their coordinates
        nodekey is a counter for node dictionary keys

        Returns
        -------
        None.

        """
        self.nodedict = {}
        self.nodeidx = {}
        self.thalwegdict = {}
        self.ridgedict = {}
        self.nodekey = 0
    
    # constructor building the network from a terrain
    def buildThalwegsFromTerrain(self, terrain, correction = True):
        """
        Builds a surface network on a given terrain. Start by detecting saddles.
        Follows with thalwegs and ridges. Ridges are constrained to avoid 
        conflicts with thalwegs.

        Parameters
        ----------
        terrain : Terrain
            Instance of the class Terrain on which the network is built.
        correction : boolean
            True if a correction is done to stay closer to the gradient
            False if the descent direction is not corrected as in D8

        Returns
        -------
        None.

        """
        self.nodedict = {}
        self.nodeidx = {}
        self.thalwegdict = {}
        self.ridgedict = {}
        self.terrain = terrain
        
        if correction:
            df.distancefactor = 2/3
        else:
            df.distancefactor = sqrt(2)/2
        #print('Computing all potential saddles')
        start = perf_counter_ns()
        saddledict, saddleidx = self.terrain.computeSaddles()
        #print('Removing conflicting saddles')
        self.mergeSaddles(saddledict, saddleidx)
        end = perf_counter_ns()
        print("Saddle computation time:", end - start)
        start = perf_counter_ns()
        terrain.saddleAndFlowAwareDiagonalisation(self.nodedict)
        end = perf_counter_ns()
        print("Diagonalisation time:", end - start)
        start = perf_counter_ns()
        self.computePits()
        end = perf_counter_ns()
        print("Pit computation time:", end - start)
        start = perf_counter_ns()
        self.createThalwegs(correction)
        self.orderThalwegsAroundNodes()
        end = perf_counter_ns()
        print("Thalweg computation time (including ordering):", end - start)
        start = perf_counter_ns()
        self.computeHills()
        end = perf_counter_ns()
        print("Hill computation time:", end - start)
        start = perf_counter_ns()
        terrain.saddleAndThalwegAwareDiagonalisation()
        end = perf_counter_ns()
        print("Diagonalisation time:", end - start)
        start = perf_counter_ns()
        self.computePeaks()
        end = perf_counter_ns()
        print("Peak computation time:", end - start)
        
    def buildFlowThalwegsFromTerrain(self, terrain, correction = False):
        """
        Builds a surface network on a given terrain. Start by detecting saddles.
        Follows with thalwegs and ridges. Ridges are constrained to avoid 
        conflicts with thalwegs.

        Parameters
        ----------
        terrain : Terrain
            Instance of the class Terrain on which the network is built.
        correction : boolean
            True if a correction is done to stay closer to the gradient
            False if the descent direction is not corrected as in D8

        Returns
        -------
        None.

        """
        self.nodedict = {}
        self.nodeidx = {}
        self.thalwegdict = {}
        self.ridgedict = {}
        self.terrain = terrain
        
        print("Diagonalising the raster")
        terrain.flowDiagonalisation()
        print('Computing saddles and pits')
        start = perf_counter_ns()
        self.computeSaddles()
        self.computePits()
        end = perf_counter_ns()
        print("Saddle and pit computation time:", end - start)
        start = perf_counter_ns()
        self.createThalwegs(correction)
        self.orderThalwegsAroundNodes()
        end = perf_counter_ns()
        print("Thalweg computation time (including ordering):", end - start)
        start = perf_counter_ns()
        self.computeHills()
        end = perf_counter_ns()
        print("Hill computation time:", end - start)
        start = perf_counter_ns()
        terrain.saddleAndThalwegAwareDiagonalisation()
        end = perf_counter_ns()
        print("Diagonalisation time:", end - start)
        start = perf_counter_ns()
        self.computePeaks()
        end = perf_counter_ns()
        print("Peak computation time:", end - start)
        
    def saveNetwork(self, directory, name):
        tmp = self.terrain
        self.terrain = None # we don't save the terrain, that's too big
        with open( directory + name + '.snw', 'wb') as f:
            pickle.dump(self, f, pickle.HIGHEST_PROTOCOL)
        self.terrain = tmp
    
    def getThalwegNeighbourNodes(self, ipt):
        """
        Returns nodes connected to this node through the thalwegs. Two sets are
        returned: one containing higher neighbours and one containing lower neighbours.

        ipt: key of a node defined in nodedict
        
        Returns
        -------
        (lowerneighbour, higherneighbour): two sets containing all neighbouring node ids.

        """
        pt = self.nodedict[ipt]
        higherneighbour = set()
        lowerneighbour = set()
        for it in pt['thalweg']:
            if it > 0:
                lowerneighbour.add(self.thalwegdict[it]['end'])
            else:
                higherneighbour.add(self.thalwegdict[-it]['start'])
        return lowerneighbour, higherneighbour
    
    def computeThalwegLength(self, it):
        """
        Compute the length of a thalweg

        Parameters
        ----------
        it : integer, thalweg index

        Returns
        -------
        length: length of the thalweg.

        """
        t = self.thalwegdict[it]
        return df.length(t['polyline'])
    
    def computeThalwegSlope(self, it):
        t = self.thalwegdict[it]
        if t['end'] == -1:
            return df.NOTIN
        length = df.length(t['polyline'])
        ptend = self.nodedict[t['end']]
        ptstart = self.nodedict[t['start']]
        dz = ptstart['z'] - ptend['z']
        return dz / length
        
    def computeSlope(self, ipt1, ipt2):
        """
        Computes the slope in a straight line between two nodes.
        It is done by dividing the difference of elevation by the distance

        Parameters
        ----------
        ipt1 : integer
            index of the first node.
        ipt2 : integer
            index of the second node.

        Returns
        -------
        The signed slope. The value is negative if the second point is lower 
        than the first point.

        """
        pt1 = self.nodedict[ipt1]
        pt2 = self.nodedict[ipt2]
        
        if ipt1 < 0:
            return -pt1['z']
        dz = pt2['z'] - pt1['z']
        if dz == 0:
            return 0
        delta = (pt2['ij'][0] - pt1['ij'][0], pt2['ij'][1] - pt1['ij'][1])
        distance = np.sqrt(delta[0]*delta[0] + delta[1]*delta[1])
        
        return dz / distance
        
    def mergeSaddles(self, saddledict, saddleidx):
        """
        Eliminates saddles in conflict with other saddles. The method proceeds
        first by building clusters of saddles.

        Parameters
        ----------
        saddledict : dictionary
            List of saddles previously detected on the terrain.
        saddleidx : dictionary
            Hash table on saddle coordinates.

        Returns
        -------
        None.

        """
        m,n = self.terrain.m, self.terrain.n
        pixelclass = self.terrain.pixelclass
        passdict = {}
        passid = 0
        
        # a pass is a group of (one or more) D4-connected saddles
        # we build up all passes
        for idx, pxl in saddledict.items():
            (i,j) = pxl['ij']
            if pixelclass[i,j] == df.TOBEPROCESSED:
                saddlestack = [idx]
                pixelclass[i, j] = df.SADDLE
                passdict[passid] = {
                    'node' : [],
                    'line' : [],
                }
                while saddlestack:
                    idx = saddlestack.pop()
                    (i,j) = saddledict[idx]['ij']
                    passdict[passid]['node'].append(idx)
                    neighbours4 = [(i+di, j+dj) for (di,dj) in df.ldr4]
                    for k in neighbours4:
                        if (k[0] >= 0 and k[0] < m and k[1] >= 0 and k[1] < n):
                            if (pixelclass[k] == df.TOBEPROCESSED):
                                saddlestack.append(saddleidx[k])
                                pixelclass[k] = df.SADDLE
                passid += 1

        # get the list of all passes composed of several saddles
        multisaddle = [p for p in passdict if len(passdict[p]['node'])>1]
        while multisaddle:
            p = multisaddle.pop() # p is a cluster of saddles
            saddlelist = passdict[p]['node']
            saddleneighbour = {}
            for s in saddlelist:
                saddleneighbour[s] = list(df.ldr8)
            #saddleneighbour
            while len(saddlelist)>1:
                conflictsaddle = [] # list of saddles in conflict with minsaddle
                # we get the lowest saddle
 #               iminsaddle = max(saddlelist, key=lambda x: saddledict[x]['z'])
                iminsaddle = min(saddlelist, key=lambda x: saddledict[x]['z'])
                minsaddle = saddledict[iminsaddle]
                # we get its diagonal lines
                linelist = minsaddle['line']
                # we keep the lines that are in diagonal
                dline = [l for l in linelist if abs(l[1][0]-l[0][0])+abs(l[1][1]-l[0][1]) == 2]
                if not dline: # if lines of a saddle are all edge-aligned, it is not in conflict
                    saddlelist.remove(iminsaddle)
                    continue
                # we compute other saddles intersecting its lines
                removable = True
                for si in saddlelist:
                    if si == iminsaddle:
                        continue
                    saddle = saddledict[si]
                    conflictline = []
                    for l1 in saddle['line']:
                        for l2 in dline:
                            if df.crossingSegment(l1[0:2], l2[0:2]):
                                removable = False
                                conflictline.append(l1)
                    if conflictline:
                        conflictsaddle.append((si, conflictline))
                if removable: # the saddle is not in conflict, it is removed from the cluster
                    # lines from dline must be preserved if other saddles are modified
                    for l in dline:
                        # get l endpoints
                        p1 = l[0]
                        p2 = l[1]
                        # take the points on the other diagonal
                        q1 = (p1[0], p2[1])
                        q2 = (p2[0], p1[1])
                        # check if there is a saddle at each point
                        for si in saddlelist:
                            saddle = saddledict[si]
                            if saddle['ij'] == q1:
                                dq = (q2[0] - q1[0], q2[1] - q1[1])
                                saddleneighbour[si].remove(dq)
                            elif saddle['ij'] == q2:
                                dq = (q1[0] - q2[0], q1[1] - q2[1])
                                saddleneighbour[si].remove(dq)
                            
                    saddlelist.remove(iminsaddle)
                    continue
                # we have the list of saddles in conflict with minsaddle
                # and we know which line is in conflict
#                print(iminsaddle, conflictsaddle)
                for conflict in conflictsaddle:
                    si = conflict[0]
                    conflictlist = conflict[1]
                    for conflictline in conflictlist:
                        p = conflictline[0]
                        q = conflictline[1]
                        pq = (q[0]-p[0],q[1]-p[1])
                        saddleneighbour[si].remove(pq) # remove direction from possible line directions
                    (i,j) = saddledict[si]['ij']
                    # check if the point is still a saddle with the remaining directions
                    criticalline = self.terrain.checkSaddle(i,j,saddleneighbour[si])
                    if criticalline:
#                        print('cannot remove', si)
                        saddledict[si]['line'] = criticalline
                    else:
                        saddlelist.remove(si)
                        del saddledict[si]
                        self.terrain.pixelclass[i,j] = df.SLOPE
                saddlelist.remove(iminsaddle)
#                break
        
        self.nodekey = 0
        for i in saddledict:
            self.nodedict[self.nodekey] = saddledict[i]
            self.nodekey += 1
        del saddledict
#        self.saddlelist = [*passdict] # define list of saddle keys

    def computeSaddles(self):
        """
        Detection of saddles from the terrain. The terrain has already been diagonalised
        Pits are detected at the same time. Peaks are not because the terrain may be
        rediagonalised to optimise peak and ridge computation

        Returns
        -------
        saddledict : dictionary
            Contains all the saddles with their critical lines.
        saddleidx : dictionary
            Hash table storing saddle coordinates for indexing.

        """
        # matrix dimensions
        A = self.terrain.dtm
        pixelclass = self.terrain.pixelclass
        m = self.terrain.m
        n = self.terrain.n

        # these are counters for counting saddle points
        count_p = 0  # number of saddles (counting multiplicity)
        count_n = 0  # number of saddles (without multiplicity)

        saddledict = {}  # dictionary of saddles
        saddleidx = {}  # dictionary indexing saddle coordinates


        # here starts the computation of saddle points
        # we go through each point and check its elevation against its neighbours
        for i in range(m):
            for j in range(n):
                if pixelclass[i, j] > df.OUT:
                    z = A[i, j]
                    ldr = self.terrain.neighbour[i, j]
                    nv = len(ldr)
                    v = self.terrain.getNeighbourHeights(i, j)  # neighbour elevations
                    signe = [0] * nv  # sign of elevation difference with centre
                    Nc = 0

                    # compute the points above and below
                    for k in range(nv):
                        dr = ldr[k]
                        if (v[k] - z > 0):
                            signe[k] = 1

                        elif (v[k] - z < 0):
                            signe[k] = -1

                        else:
                            # to be processed separately: flat areas?
                            if (dr[0] > 0):  # di > 0
                                signe[k] = 1
                            elif (dr[0] < 0):
                                signe[k] = -1
                            elif (dr[1] > 0):
                                signe[k] = 1
                            else:
                                signe[k] = -1
                    # compute how many times we go above and below
                    for k in range(nv):
                        if (signe[k] != signe[k - 1]):
                            Nc += 1

                    if Nc >= 4:  # classified as a saddle
                        dz = df.gradLength(z, v, ldr)
                        # compute the points above and below
                        indx = list(range(nv))
                        # if two consecutive neighbours have the same sign, we keep the one with the steepest slope
                        k = nv - 1
                        while k >= 0:
                            if (signe[k] == signe[k - 1]):
                                tmp = signe[k] * (dz[k] - dz[k - 1])
                                if tmp == 0:
                                    # we fix the sign according to the lexicographical order
                                    ik = indx[k]
                                    ik1 = indx[k - 1]
                                    drk = ldr[ik]
                                    drk1 = ldr[ik1]
                                    tmp = -signe[k] * df.lexicoSign(drk1[0] - drk[0], drk1[1] - drk[1])
                                if (tmp > 0):
                                    del signe[k - 1]
                                    del dz[k - 1]
                                    del indx[k - 1]
                                elif (tmp < 0):
                                    del signe[k]
                                    del dz[k]
                                    del indx[k]
                                if k == len(signe):
                                    k -= 1
                            else:
                                k -= 1
                        criticalline = [[(i, j), (i + di, j + dj), 0, 0] for (di, dj) in [ldr[k] for k in indx]]
                        k = 0
                        for l in criticalline:
                            l[3] = v[indx[k]]
                            k += 1
                        firstline = criticalline[0]
                        if firstline[3] > z:  # first line is a ridge
                            s = 1
                        elif firstline[3] < z:
                            s = -1
                        else:
                            s = df.lexicoSign(firstline[1][0] - firstline[0][0], firstline[1][1] - firstline[0][1])
                        for l in criticalline:
                            l[2] = s
                            s = -1 * s
                        saddledict[count_n] = {
                            'ij': (i, j),
                            'z': A[i, j],
                            'line': criticalline,
                            'thalweg': [],
                            'ridge': [],
                            'type': df.SADDLE
                        }
                        pixelclass[i, j] = df.SADDLE
                        saddleidx[i, j] = count_n
                        count_n += 1
                        count_p += (Nc / 2 - 1)
        self.nodekey = 0
        for i in saddledict:
            self.nodedict[self.nodekey] = saddledict[i]
            self.nodekey += 1
        del saddledict

    # computation of pits considering saddles
    def computePits(self):
        """
        Compute all pits on the raster. To be done after the triangulation.
        This is not necessary to compute the thalwegs but it helps to check the
        consistency of the result. It also adds pits to the node dictionary and
        add the virtual pit. Only one virtual pit is created. If holes are handled,
        there should be one virtual pit per hole (or one virtual node, not 
        necessarily a pit depending how it is defined).

        Returns
        -------
        None.

        """
        # matrix dimensions
        m = self.terrain.m
        n = self.terrain.n
        pixelclass = self.terrain.pixelclass

        # here starts the computation of pit points
        # we go through each point and check its elevation against its neighbours
        for i in range(m):
            for j in range(n):
                if (pixelclass[i,j] == df.SLOPE):
                    z = self.terrain.dtm[i,j]
                    ldr = self.terrain.neighbour[i,j]
                    nv = len(ldr)
                    v=self.terrain.getNeighbourHeights(i,j) # elevation of 8 neighbours
                    delta_moins = 0
                 
                    # compute the points above and below
                    for k in range(nv) :
                        dr=ldr[k]
                        if (z-v[k])>0:
                            delta_moins = 1
                            break 
                        if z == v[k]:
                            # to be processed before: flat areas
                            if (dr[0]<0):
                                delta_moins = 1
                                break
                            if (dr[0] == 0) and (dr[1]<0):
                                delta_moins = 1
                                break
                                
                    if delta_moins == 0: #classified as a pit
                        pixelclass[i,j] = df.PIT
                        self.nodedict[self.nodekey] = {
                            'ij': (i,j),
                            'z': self.terrain.dtm[i,j],
                            'ridge': [],
                            'thalweg': [],
                            'type': df.PIT
                        }
                        self.nodekey += 1

        # create an index of points to nodedict
        self.nodeidx = {}
        for i, node in self.nodedict.items():
            self.nodeidx[node['ij']] = i

        # add the virtual pit
        self.nodedict[-1] = {
            'ij' : [],
            'thalweg' : [], 
            'ridge' : [],
            'type' : df.PIT,
            'z' : self.terrain.nodata
        }

    # computation of peaks considering saddles and thalwegs
    def computePeaks(self):
        """
        Compute all peaks on the raster. To be done after the triangulation and
        normally, after thalweg computation and diagonal flipping.
        This is not necessary to compute the ridges but it helps to check the
        consistency of the result. It also adds peaks to the node dictionary . 
        If holes are handled, some virtual peaks may be created also.

        Returns
        -------
        None.

        """
        # matrix dimensions
        m = self.terrain.m
        n = self.terrain.n
        pixelclass = self.terrain.pixelclass

        # here starts the computation of peak points
        # we go through each point and check its elevation against its neighbours
        for i in range(m):
            for j in range(n):
                if (pixelclass[i,j] == df.SLOPE):
                    z = self.terrain.dtm[i,j]
                    ldr = self.terrain.neighbour[i,j]
                    nv = len(ldr)
                    v=self.terrain.getNeighbourHeights(i,j) # elevation of 8 neighbours
                    delta_plus = 0
                 
                    # compute the points above and below
                    for k in range(nv) :
                        dr=ldr[k]
                        if (z-v[k])<0 : # higher neighbour, cannot have a peak
                            delta_plus = 1
                            break
                        if z == v[k]: # same height but on the right or above is higher
                            if (dr[0]>0):
                                delta_plus = 1
                                break
                            if (dr[0] == 0) and (dr[1]>0):
                                delta_plus = 1
                                break
                                
                    if delta_plus == 0: #classified as a peak
                        pixelclass[i,j] = df.PEAK
                        self.nodedict[self.nodekey] = {
                            'ij': (i,j),
                            'z': self.terrain.dtm[i,j],
                            'ridge': [],
                            'thalweg': [],
                            'type': df.PEAK
                        }
                        self.nodeidx[(i,j)] = self.nodekey
                        self.nodekey += 1

    #creation of the thalweg dictionary
    def createThalwegs(self, correction = True):
        """
        Computation of thalwegs starting from saddles. Thalwegs do not overlap or
        go through saddles. If two thalwegs join, a confluence is inserted and a
        new thalweg runs from the confluence to the next node. Thalwegs go
        normally along the line of steepest descent. Since we stick to the vertices,
        we go to the closest vertice but record the offset and adjust it when
        too big by shifting the thalweg.

        Parameters
        ----------
        correction : boolean
            True if a correction is done to stay closer to the gradient
            False if the descent direction is not corrected as in D8

        Returns
        -------
        None.

        """
        # no thalwegkey 0 because key sign indicates direction
        thalwegkey=1
        nextthalweg = thalwegkey + 1
        self.floorkey={}
        
        # thalwegs all start at saddles
        saddlelist = [i for i in self.nodedict if self.nodedict[i]['type'] == df.SADDLE]
        for ipt in saddlelist:
            pt = self.nodedict[ipt]
#            if ipt>14:
#                continue
#            print('pass ', ipt)
            linelist = pt['line']
            for l in linelist:
                if l[2] != -1: # if the line is not a thalweg
                    continue
                    
                polylist = l[0:2] # the first two points are already computed
                (i,j) = polylist[1]
#                print('creating thalweg', thalwegkey, i,j)
                keepgoing = True # no node yet
                outbound = False # not out of the terrain
                endindex = -1
                # check if we are out of the domain or on a node
                if self.terrain.isInside(i,j):
                    if self.terrain.pixelclass[i,j] == df.SLOPE:
                        keepgoing = True
                        self.terrain.pixelclass[i,j] = df.THALWEG
                        self.floorkey[i,j] = thalwegkey
                    else:
                        keepgoing = False
                else:
                    keepgoing = False
                    outbound = True # we’ll connect the thalweg to the virtual pit
                # if not a node and inbound, we look for the next point
                oldlag = 0
                oldldr = (0,0)
                while keepgoing:
                    ldr = self.terrain.neighbour[i,j]
                    v = self.terrain.getNeighbourHeights(i, j)
                    z = self.terrain.dtm[i,j]

                    # gradient computation
                    # gradLength does not compute the gradient for each triangle
                    # but the dz for each node along an edge and 3/2dz for each
                    # node along a diagonal
                    dz = df.gradLength(z, v, ldr)
                    minimum = 1
                    kmin = -1
                    for k,value in enumerate(dz):
                        if value<minimum:
                            minimum=value
                            kmin=k
                        elif value==minimum:
                            if df.lexicoSign(ldr[k][0]-ldr[kmin][0], ldr[k][1]-ldr[kmin][1])<0:
                                minimum=value
                                kmin=k

                    # the steepest downward gradient is for v[kmin]
                    # we check if we need to modify it according to the lag
                    
                    ########### beginning of the correction block #############
                    if correction:            
                        # there are two possible triangles, one on each side of the edge
                        # we first check in which triangle the gradient is
                        lag = 0
                        vleft = 0
                        vright = 0
                        kleft = kmin - 1
                        kright = kmin + 1
                        if kright == len(ldr):
                            kright = 0
                        # cases where the diagonal does not connect to the node z
                        if abs(ldr[kmin][0]-ldr[kleft][0])+abs(ldr[kmin][1]-ldr[kleft][1]) == 1:
                            vleft = v[kleft] - z
                        if abs(ldr[kmin][0]-ldr[kright][0])+abs(ldr[kmin][1]-ldr[kright][1]) == 1:
                            vright = v[kright] - z
    
                        if vleft == vright or min(vleft, vright) >= 0: # no lag, it goes along ldr[kmin]
                            lag = 0
                        else: # either to the left or the right triangle
                            k1 = kmin
                            if vleft < vright: # if the left point is lower, gradient is on the left
                                k1 = kleft
                                v1 = vleft
                            if vright < vleft:
                                k1 = kright
                                v1 = vright
                            # here we compute the gradient direction by solving a 2x2 system
                            # may be simplified noting that one coefficient is zero, others are ones
                            # but would need to separate all cases
                            x1 = ldr[kmin][0]
                            y1 = ldr[kmin][1]
                            x2 = ldr[k1][0]
                            y2 = ldr[k1][1]
                            z1 = v[kmin] - z
                            z2 = v1
                            if z1 != self.terrain.nodata and z2 != self.terrain.nodata:
                                det = 1./(x1*y2 - x2*y1)
                                a = det*(y2*z1 - y1*z2)
                                b = det*(x1*z2 - x2*z1)
                                # the gradient direction is given by (a,b)
                                # now we compute the lag
                                # the lag is the difference between ldr[kmin] and (a,b)
                                if (ldr[kmin][0] == ldr[k1][0]): # compute b for a = 1
                                    b = -b/a
                                    a = ldr[kmin][0]
                                    if abs(b) < abs(a):
                                        lag = b - ldr[kmin][1]
                                    else:
                                        lag = 0
                                else:
                                    a = -a/b
                                    b = ldr[kmin][1]
                                    if abs(a) < abs(b):
                                        lag = a - ldr[kmin][0]
                                    else: 
                                        lag = 0
                            else:
                                lag = 0
                            # if we keep going the same way, correct with the lag
                            if (ldr[kmin] == oldldr):
                                lag += oldlag
                                if lag > 0.5:
                                    kmin = k1
                                    lag -= 1
                                if lag < -0.5:
                                    kmin = k1
                                    lag += 1
                            oldldr = ldr[kmin]
                            oldlag = lag

                    ############# end of the correction block #################

                    i += ldr[kmin][0]
                    j += ldr[kmin][1]
#                    print(i,j)
               
                    # check again if the point is a node or out
                    if self.terrain.isInside(i,j):
                        if self.terrain.pixelclass[i,j] == df.SLOPE:
                            keepgoing = True
                            self.terrain.pixelclass[i,j] = df.THALWEG
                            self.floorkey[i,j] = thalwegkey
                        else:
                            keepgoing = False
                            outbound = False
                    else:
                        keepgoing = False
                        outbound = True

                    # this point is on the thalweg and is added to the polyline
                    polylist.append((i,j))

                # if we reach a node
                # the thalweg id is marked <0 at its end point
                if outbound:
                    z=self.terrain.nodata
                    self.nodeidx[(i,j)] = -1 # virtual pit
                    self.nodedict[-1]['ij'].append((i,j))
                    endindex = -1
                else:
                    if self.terrain.pixelclass[i,j] == df.THALWEG:
                        # if we reach an existing thalweg, we have to create a confluence
                        idthalweg = self.floorkey[i,j]
                        # we create the new node
                        endindex = self.nodekey
                        self.nodeidx[(i,j)] = self.nodekey
                        self.nodedict[self.nodekey] = {
                            'ij' : (i,j),
                            'z': self.terrain.dtm[i,j],
                            'thalweg' : [-idthalweg, nextthalweg],
                            'ridge' : [],
                            'type' : df.CONFLUENCE
                        }
                        self.terrain.pixelclass[i,j] = df.CONFLUENCE
                        self.nodekey += 1
                        # split the existing thalweg in two
                        t = self.thalwegdict[idthalweg]
                        ptlist = t['polyline']
                        iend = t['end']
                        icfl = ptlist.index((i,j))
                        ptlist1 = ptlist[:icfl+1]
                        ptlist2 = ptlist[icfl:]
                        # modify the first thalweg to stop at the confluence
                        t['polyline'] = ptlist1
                        t['end'] = endindex
                        # create a new thalweg
                        self.thalwegdict[nextthalweg] = {
                            'start' : endindex,
                            'polyline' : ptlist2,
                            'end' : iend
                        }
                        # we update the node at the end of the new thalweg
                        ndt = self.nodedict[iend]['thalweg']
                        ndt.remove(-idthalweg)
                        ndt.append(-nextthalweg)
                        # we update floorkey with nextthalweg
                        for k in ptlist2[1:-1]:
                            self.floorkey[k] = nextthalweg
                        nextthalweg += 1
# this can be simplified, all non thalweg cases are treated the same way
                    elif self.terrain.pixelclass[i,j] == df.CONFLUENCE: 
                        # reach a confluence: add the thalweg to the confluence
                        endindex = self.nodeidx[(i,j)]
                    elif self.terrain.pixelclass[i,j] == df.SADDLE:
                        endindex = self.nodeidx[(i,j)]
                    elif self.terrain.pixelclass[i,j] == df.PIT:
                        endindex = self.nodeidx[(i,j)]

                # the last point is a critical point
                # we add its index to the thalweglist
                pt['thalweg'].append(thalwegkey)
                self.nodedict[endindex]['thalweg'].append(-thalwegkey)
                self.thalwegdict[thalwegkey]={
                    'start' : ipt,
                    'polyline' : polylist,
                    'end' : endindex,
                }
#                    print(ipt, endindex, pt)
#                    l[3] = thalwegkey # we don't use it
                thalwegkey = nextthalweg
                nextthalweg += 1

    # if it is not, we should use something like Freeman chaincode to order points along
    # the domain boundary
    def orderThalwegsAroundNodes(self):
        """
        Orders the list of thalwegs around each node in an anticlockwise order except
        for the virtual pit where the order is reversed (it is in anticlockwise order
        if we see it from under). This is required to compute the hill around each peak:
        by moving from one thalweg to its next thalweg in a loop, we can build the ring
        that delineates a hill and in which a ridge must be contained.
        
        To be done: holes are not considered.

        Returns
        -------
        None.

        """
        for ipt, pt in self.nodedict.items():
            # -1 is the id of the virtual pit outside the domain
            if ipt == -1: # to be handled separately
                thalweglist = pt['thalweg']
                boundary = self.terrain.boundary
                inbound = self.terrain.inbound
                try:
                    thalweglist.sort(key = lambda it: df.getBoundaryIndex(self.thalwegdict[abs(it)]['polyline'], boundary, inbound))
                except:
                    print('Exception orderThalwegs around virtual pit')
                    break
            elif ipt < -1:
                # holes must be virtual pits with id<-1
                # get all the neighbours around the saddles
                continue
            else:
                (i,j) = pt['ij']
                neighbours = [(i+di, j+dj) for (di,dj) in df.ldr8]
                # arrange the lines in the same order as the neighbours
                thalweglist = pt['thalweg']
                try:
                    thalweglist.sort(key = lambda it: neighbours.index(self.thalwegdict[abs(it)]['polyline'][df.getSecondPointPosition(it)]))
                except ValueError:
                    print('Exception raised in orderThalwegsAroundNodes', ipt, pt, neighbours)


    # returns lists of thalwegs forming a polygon (including holes) starting at idnode along idthalweg and a list of dangling thalwegs inside the polygon
    def getThalwegRing(self, idnode, idthalweg):
        """
        For a given node and a given thalweg starting at this node, the function
        returns the list of thalwegs that delineate a hill when walking from this
        thalweg in a clockwise order. The function takes into account overlapping
        segments although it has not been checked thoroughly (if there are overlaps
        inside a ring).

        Parameters
        ----------
        idnode : integer
            Node at which the hill starts.
        idthalweg : integer
            First thalweg in the list, index can be negative depending on the
            thalweg direction.

        Returns
        -------
        thalwegring a list of lists of oriented thalwegs forming, the first list
        is the boundary ring, other lists are ring holes
        thalwegline a list of thalwegs representing dangling thalwegs inside the 
        hill. They do not form a ring.

        """
        # case of a line and no polygon has not been tested yet
        idthalweg = abs(idthalweg)
        currentnode = idnode
        currentthalweg = idthalweg
        keepgoing = True
        thalwegring = []
        while keepgoing:
            t = self.thalwegdict[currentthalweg]
            if currentnode == t['start']:
                nextnode = t['end']
                thalwegring.append(currentthalweg)
            elif currentnode == t['end']:
                nextnode = t['start']
                thalwegring.append(-currentthalweg)
            else:
                print('error in getPolygon: thalweg', currentthalweg, 'and node', currentnode, 'do not match')

    #        print(currentthalweg)
            currentnode = nextnode
            thalweglist = [abs(t) for t in self.nodedict[currentnode]['thalweg']]
            thalwegpos = thalweglist.index(currentthalweg)
            currentthalweg = thalweglist[thalwegpos - 1]
            keepgoing = not ((currentnode == idnode) and (currentthalweg == idthalweg))
        # if we have overlapping segments
        # we remove these overlapping segments to split the ring in
        # the outer and inner rings and inner lines
        
        nt = len(thalwegring)
        tbr = []
        innerline = []
        innerring = []
        
        it1=0
    #       print thalwegring
        while it1 < nt:
            t1 = thalwegring[it1]
            it2 = nt-1
    #           print 't1 = ', t1
            while it2 > it1:
                t2 = thalwegring[it2]
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
                        t1 = thalwegring[it1]
                        t2 = thalwegring[it2]
                    if (t2 != -t1):
                        innerring.append(thalwegring[it1:it2+1])
                        tbr.extend(innerring[-1])
                    it1 -= 1
                    break
    # bug: we should go iteratively in the innerring to check if any other overlapping thalweg
                it2 -= 1
            it1 += 1

        thalwegboundary = [t for t in thalwegring if t not in tbr]
        if thalwegboundary == []:
            thalwegring = innerring
        else:
            thalwegring = [thalwegboundary] + innerring
        for ir, ring in enumerate(thalwegring): # this corrects a bug but to be checked if we can do better
            for ir2, ring2 in enumerate(thalwegring[ir+1:]):
                for t in ring2:
                    if t in ring:
                        ring.remove(t)
        for t in innerline: # this may correct the bug but is a fix, to be checked
            for ring in thalwegring:
                if t in ring:
                    ring.remove(t)
                if -t in ring:
                    ring.remove(-t)
        while [] in thalwegring:
            thalwegring.remove([])
        return thalwegring, innerline

    def createHillPolygons(self):
        """
        Creates a shapely polygon from a set of thalwegs forming a hill. The 
        method is used to maintain a ridge within a hill.

        Parameters
        ----------
        thalwegring : list of integers
            List of thalwegs indexes delineating a hill.
        thalwegline : list of integers
            List of thalwegs forming dangling lines inside the hill.

        Returns
        -------
        hillpolygon : list of polygons
            Polygons of hills defined as a shapely geometry
        danglingmap : np.array
            map where each pixel has the id of the dangling thalweg passing there

        """
        danglingmap = np.zeros((self.terrain.m, self.terrain.n))
        hillpolygon = {}
        for ihill in self.hilldict:
            hill = self.hilldict[ihill]
            thalwegring = [hill['boundary']] + hill['hole']
            ringlist = []
            plgnlist = []
            arealist = []

            for tring in thalwegring:
                if not tring:
                    continue
                pring = []
                # get the list of points forming the ring
                for it in tring:
                    if it > 0:
                        t = self.thalwegdict[it]
                        if pring:
                            pring.extend(t['polyline'][1:])
                        else:
                            pring.extend(t['polyline'][:])
                    else:
                        t = self.thalwegdict[-it]
                        if t['end'] == -1:
        #                        print('it', it, pring, t['polyline'][-1])
                            if not pring:
                                lastpoint = self.thalwegdict[abs(tring[-1])]['polyline'][-1]
                                istart = self.terrain.boundary.index(lastpoint)
                                iend = self.terrain.boundary.index(t['polyline'][-1])
                            else:
                                istart = self.terrain.boundary.index(pring[-1])
                                iend = self.terrain.boundary.index(t['polyline'][-1])
                            if istart > iend:
                                pring.extend(self.terrain.boundary[istart-1:iend-1:-1])
                            elif istart < iend:
                                pring.extend(self.terrain.boundary[istart-1::-1])
                                pring.extend(self.terrain.boundary[-1:iend-1:-1])
        #                    print(istart, iend)
                        tmp = []
                        if pring:
                            tmp = t['polyline'][:-1]
                        else:
                            tmp = t['polyline'][:]
                        tmp.reverse()
                        pring.extend(tmp)
                ringlist.append(pring)
                plgn = shapely.geometry.Polygon(pring)
                arealist.append(plgn.area)
                plgnlist.append(plgn)
            boundaryid = max(enumerate(arealist), key = lambda a: a[1])[0]
            boundary = ringlist.pop(boundaryid)
        #    print(boundary)
        #    print(ringlist)
            plgn = shapely.geometry.Polygon(boundary, ringlist)
            
            hillpolygon[ihill] = plgn
            
            danglingthalweg = hill['dangling']
            for it in danglingthalweg:
                t = self.thalwegdict[it]
                polyline = t['polyline'][:-1]
                for p in polyline:
                    danglingmap[p] = it
        return hillpolygon, danglingmap

    def computeHills(self):
        """
        Computes the dictionary of hills of the surface network. A hill is a
        descending manifold. It is formed by a series of thalwegs. The hill can
        contain some holes.
        Each hill is defined by a series of thalwegs forming the boundary and
        a series of holes formed each by a list of thalwegs. Thalwegs are oriented
        in the lists. The hill also contains a list of dangling thalwegs, thalwegs
        that do not form a loop.
        The function also adds a pointer to the left and right hill to thalwegs
        in thalwegdict. For a dangling thalweg, the left and right hills are the
        same.

        Returns
        -------
        None.

        """
        # initialise the link from thalweg to hill
        for it in self.thalwegdict:
            t = self.thalwegdict[it]
            t['lefthill'] = None
            t['righthill'] = None
    
        self.hilldict = {}
        hillkey = 0
        saddlelist = [i for i in self.nodedict if self.nodedict[i]['type'] == df.SADDLE and self.nodedict[i]['z']>self.terrain.nodata]
    
        for ipt in saddlelist:
            pt = self.nodedict[ipt]
            # we get the list of thalwegs
            thalweglist = pt['thalweg']
            # we need at least two thalwegs to create a hill
            # if we don't have this, there's a problem in the data
            # a saddle always has two thalwegs at least
            if len(thalweglist)<2:
                print('less than two thalwegs', ipt, thalweglist)
                continue
            # for each thalweg, we build a new hill on its right
            for t in thalweglist:
                # computes the ring and dangling thalwegs starting from t at ipt
                abst = abs(t)
                if self.thalwegdict[abst]['righthill']:
                    continue
                thalwegring, thalwegline = self.getThalwegRing(ipt, t)
                # if t is a dangling thalweg, no new hill is created because there's none
                if abst in thalwegline:
                    continue
                # thalwegring can include several rings, one of them is the boundary
                # the others are holes
                # we compare the bounding boxes to check if the first ring is the boundary
                arealist = []
                for ring in thalwegring:
                    pointlist = set(reduce(lambda a,b:a+b, [self.thalwegdict[abs(tt)]['polyline'] for tt in ring]))
                    xmin = min(pointlist, key = lambda a: a[0])[0]
                    xmax = max(pointlist, key = lambda a: a[0])[0]
                    ymin = min(pointlist, key = lambda a: a[1])[1]
                    ymax = max(pointlist, key = lambda a: a[1])[1]
                    area = (xmax - xmin) * (ymax - ymin)
                    arealist.append(area)
                # if the first ring is not the biggest, the hill is a hole in a larger one
                # we discard everything outside the hole
                maxarea = max(arealist)
                maxindex = arealist.index(maxarea)
                border = thalwegring.pop(maxindex)
                holes = thalwegring
                
                self.hilldict[hillkey] = {'boundary': border, 
                                          'hole': holes,
                                          'dangling': thalwegline}
                for tt in border:
                    if tt > 0:
                        self.thalwegdict[tt]['righthill'] = hillkey
                    else:
                        self.thalwegdict[-tt]['lefthill'] = hillkey
                for r in holes:
                    for tt in r:
                        if tt > 0:
                            self.thalwegdict[tt]['righthill'] = hillkey
                        else:
                            self.thalwegdict[-tt]['lefthill'] = hillkey
                for tt in thalwegline:
                    self.thalwegdict[tt]['lefthill'] = hillkey
                    self.thalwegdict[tt]['righthill'] = hillkey
                hillkey += 1

        # correction of a bug from hill computation where some hills were computed twice
        nopeak = set()
        for ihill in self.hilldict:
            tlist = self.hilldict[ihill]['boundary']
            for it in tlist:
                t = self.thalwegdict[abs(it)]
                if t['lefthill'] != ihill and t['righthill'] != ihill:
                    nopeak.add(ihill)
                    break
        for i in nopeak:
            del self.hilldict[i]

