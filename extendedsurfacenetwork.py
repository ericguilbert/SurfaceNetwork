# -*- coding: utf-8 -*-
"""
Created on Tue Aug 5 21:30:05 2020

@author: ERGUI19
"""
from time import perf_counter_ns
from math import sqrt, atan, pi
from collections import deque
import numpy as np

import shapely
import deffunction as df
from celltype import CellType
from streamnetwork import StreamNetwork

class ExtendedSurfaceNetwork(StreamNetwork):
    """
    ExtendedSurfaceNetwork stores a ESN defined by a set of thalwegs and a
    set of ridges that are the drainage divides. It inherits from StreamNetwork 
    which inherits from SurfaceNetwork and adds functions for flow direction and
    accumulation computation.
    ExtendedSurfaceNetwork adds ridges at pits and confluences, allowing
    for the definition of a dale around each thalweg and the computation of 
    flow accumulation and drainage bassins.
    """
    
    def buildESNFromTerrain(self, terrain, sd = True, smooth = False):
        """
        Builds a ESN network on a given terrain. Starts by building the
        thalweg network and the hills. Follows with ridges and dales.

        Parameters
        ----------
        terrain : Terrain
            Instance of the class Terrain on which the network is built.
        sd : boolean
            True if the network does not have diffluences (single direction)
            (not validated for sd = False)
        smooth : boolean
            True if the surface is considered smooth. Critical lines are computed
            by using the method from [Josis 21] to fit with the gradient
            False if the surface is noisy and the direction is given by the steepest
            slope computed from neighbours

        Returns
        -------
        None.

        """
        self.singleflowdirection = sd
        
        if smooth:
            df.distancefactor = 2/3
            correction = True
        else:
            df.distancefactor = sqrt(2)/2
            correction = False
        self.buildFlowThalwegsFromTerrain(terrain, correction)
        start = perf_counter_ns()
        divide = True
        print("Build ridges")
        self.instantiateRidgesAtSaddles(divide)
        # erase line property in saddles
        for ipt in self.nodedict:
            if 'criticalline' in self.nodedict[ipt]:
                del self.nodedict[ipt]['criticalline']
        self.instantiateRidgesAtPitsAndTJunctions()
        self.traceRidges(divide, correction)
        self.correctDivideToRidge()
        self.orderRidgesAroundNodes()
        end = perf_counter_ns()
        print("Ridge computation time:", end - start)
        start = perf_counter_ns()
        self.computeThalwegDales()
        end = perf_counter_ns()
        print("Dale computation time:", end - start)
        
    def getSurfaceWatercoursesAtNode(self, ipt):
        pt = self.nodedict[ipt]
        line = self.linedict
        return [i for i in pt['line'] if line[abs(i)]['type'] & CellType.SURFACE_FLAG]

    def previousSurfaceWatercourseIndexAtNode(self, ipt, it):
        pt = self.nodedict[ipt]
        line = pt['line']
        linedict = self.linedict
        n = len(line)
        jt = it - 1
        while jt != it:
            if jt == -1:
                jt = n - 1
            if linedict[abs(line[jt])]['type'] & CellType.SURFACE_FLAG:
                return jt
            jt -= 1
            
    def nextSurfaceWatercourseIndexAtNode(self, ipt, it):
        pt = self.nodedict[ipt]
        line = pt['line']
        linedict = self.linedict
        n = len(line)
        jt = it + 1
        while jt != it:
            if jt == n:
                jt = 0
            if linedict[abs(line[jt])]['type'] & CellType.SURFACE_FLAG:
                return jt
            jt += 1
            
    def getDividesAtNode(self, ipt):
        pt = self.nodedict[ipt]
        line = self.linedict
        return [i for i in pt['line'] if line[abs(i)]['type'] & CellType.DIVIDE_FLAG]

    def previousDivideIndexAtNode(self, ipt, it):
        pt = self.nodedict[ipt]
        line = pt['line']
        linedict = self.linedict
        n = len(line)
        jt = it - 1
        while jt != it:
            if jt == -1:
                jt = n - 1
            if linedict[abs(line[jt])]['type'] & CellType.DIVIDE_FLAG:
                return jt
            jt -= 1
            
    def nextDivideIndexAtNode(self, ipt, it):
        pt = self.nodedict[ipt]
        line = pt['line']
        linedict = self.linedict
        n = len(line)
        jt = it + 1
        while jt != it:
            if jt == n:
                jt = 0
            if linedict[abs(line[jt])]['type'] & CellType.DIVIDE_FLAG:
                return jt
            jt += 1
            
    def instantiateRidgesAtPitsAndTJunctions(self):
        """
        Instantiates ridge objects between two thalwegs at all pits and tjunctions. 
        
        Returns
        -------
        None.
        
        """
        # key for new ridges        
        ridgestart = self.linekey
        ridgekey = ridgestart
        divideside = {}

        # we check z>nodata to avoid starting a ridge in a hole 
        # we start ridges between two consecutive thalwegs from pits and confluences
        # pits with one thalweg are excluded
#        pitlist = [i for i in self.nodedict \
#                   if ((self.nodedict[i]['type'] == df.PIT and len(self.nodedict[i]['thalweg'])>1)\
#                   or self.nodedict[i]['type'] == df.CONFLUENCE) and self.nodedict[i]['z']>self.terrain.nodata]
        pitlist = [i for i in self.nodedict if (self.nodedict[i]['type'] & CellType.PIT_FLAG and i!=-1) or \
                   self.nodedict[i]['type'] & CellType.TJUNCTION_FLAG]
        
        for ipt in pitlist:
            #if ipt!=testipt:
            #    continue
            pt = self.nodedict[ipt]
                    
            # we get the list of thalwegs
#            thalweglist = self.getThalwegsAtNode(ipt)
            thalweglist = self.getDownslopeLinesAtNode(ipt)
            #print('Node ', ipt, thalweglist)

            (i,j) = pt['ij']
            # orderedneighbours is the neighbours of (i,j) in clockwise order
            orderedneighbours = [(i+di, j+dj) for (di,dj) in self.terrain.neighbour[i,j]]
            # for each thalweg, we start a ridge on its right
            for it, t in enumerate(thalweglist):
                #print('thalweg', it, t)
                if t > 0:
                    ihill = self.linedict[t]['righthill']
                else:
                    ihill = self.linedict[-t]['lefthill']
                hill = self.hilldict[ihill]
                
                thalwegline = hill['dangling']

                nt = len(thalweglist)-1
                # tright is the next thalweg in the list
                # because the list is in clockwise order, it is the thalweg on the right
                if it == nt:
                    tright = thalweglist[0]
                else:
                    tright = thalweglist[it+1]
                # pleft and pright are coordinates of thalwegs' second points
                if t>0:
                    pleft = self.linedict[t]['polyline'][1]
                else:
                    pleft = self.linedict[-t]['polyline'][-2]
                if tright>0:
                    pright = self.linedict[tright]['polyline'][1]
                else:
                    pright = self.linedict[-tright]['polyline'][-2]
                # ileft and iright are the neighbour pixels where the thalwegs pass through
                ileft = orderedneighbours.index(pleft)
                iright = orderedneighbours.index(pright)
                #print('tleft, tright', t, tright, ihill)
                #print('ileft, iright', ileft, iright)
                # slopeneighbours is the list of pixels between the thalwegs
                slopeneighbours = []
                if ileft > iright:
                    slopeneighbours = orderedneighbours[ileft:]+orderedneighbours[:iright+1]
                elif nt > 0:
                    slopeneighbours = orderedneighbours[ileft:iright+1]
                else:
                    slopeneighbours = orderedneighbours[ileft:]+orderedneighbours[:iright]
                slopeneighbours = [v for v in slopeneighbours if self.terrain.isInside(v[0],v[1])]
                # the ridge in this hill must go through one of these pixels
                
                # take the steepest slope in pixels of slopeneighbours
                z = self.terrain.dtm[i, j]
                slope = [2*(self.terrain.dtm[v] - z)/(abs(v[0] - i) + abs(v[1] - j) + 1) for v in slopeneighbours]
                maxslope = max(slope)
                maxcount = slope.count(maxslope)
                point = None
                if maxcount > 1:
                    maxcandidates = [slopeneighbours[i] for i in range(len(slope)) if slope[i] == maxslope]
                    maxcandidates.sort(key = lambda x: df.order8.index((x[0]-i, x[1]-j)))
                    point = maxcandidates[-1]
                else:
                    maxpos = slope.index(maxslope)
                    point = slopeneighbours[maxpos]
                polylist = [(i,j), point]
                #print('create new ridge', ipt, t, tright, slopeneighbours[maxpos])

                # if the ridge overlaps a dangling thalweg, the side of the thalweg must be recorded
                if pleft == polylist[1] and -t in thalwegline:
                    divideside[self.linekey] = -1
                if pright == polylist[1] and -tright in thalwegline:
                    divideside[self.linekey] = 1
                
                self.addNewRidge(self.linekey, ihill, ipt, t, polylist, CellType.DIVIDE)
                self.linekey += 1

        # deal with ridges starting at the virtual pit
        ipt = -1
        pt = self.nodedict[ipt]
        thalweglist = self.getDownslopeLinesAtNode(ipt)
        # orderedneighbours is the neighbours of (i,j) in clockwise order
        innerneighbours = self.terrain.inbound
        outerneighbours = self.terrain.boundary
        # for each thalweg, we start a ridge on its right
        for it, t in enumerate(thalweglist):
            # at a pit, thalweg ids are all negative
            ihill = self.linedict[-t]['lefthill']
            hill = self.hilldict[ihill]
            
            thalwegline = hill['dangling']

            nt = len(thalweglist)-1
            # tright is the next thalweg in the list
            # around the virtual pit, the order is reversed
            if it == nt:
                tleft = thalweglist[0]
            else:
                tleft = thalweglist[it+1]
            # pleft and pright are coordinates of thalwegs' second points
            pright = self.linedict[-t]['polyline'][-2]
            pleft = self.linedict[-tleft]['polyline'][-2]
            # ileft and iright are the neighbour pixels where the thalwegs pass through
            iright = innerneighbours.index(pright)
            if pright == pleft:
                ileft = len(innerneighbours) - 1 - innerneighbours[::-1].index(pleft)
            else:
                ileft = innerneighbours.index(pleft)
            #print('tleft, tright', t, tright, ihill)
            #print('ileft, iright', ileft, iright)
            # slopeneighbours is the list of pixels between the thalwegs
            slopeneighbours = []
            if ileft >= iright:
                if iright == 0:
                    slopeneighbours = innerneighbours[ileft::-1]
                else:
                    slopeneighbours = innerneighbours[ileft:iright-1:-1]
            else:
                slopeneighbours = innerneighbours[iright:]+innerneighbours[:ileft+1]
            # the ridge in this hill must go through one of these pixels
            if (slopeneighbours[0] in slopeneighbours[1:]):
                start = slopeneighbours[1:].index(slopeneighbours[0])
                slopeneighbours = slopeneighbours[start:]
            if (slopeneighbours[-1] in slopeneighbours[:-1]):
                end = slopeneighbours[:-1].index(slopeneighbours[-1])+1
                slopeneighbours = slopeneighbours[:end]
            # take the max z in pixels of slopeneighbours
            slope = [self.terrain.dtm[v[0], v[1]] for v in slopeneighbours]
            maxslope = max(slope)
            maxpos = slope.index(maxslope)
            point = slopeneighbours[maxpos]
            (i,j) = point
            
            # get the first pixel outside the domain
            qleft = self.linedict[-tleft]['polyline'][-1]
            qright = self.linedict[-t]['polyline'][-1]
            ileft = outerneighbours.index(qleft)
            iright = outerneighbours.index(qright)
            slopeneighbours = []
            if ileft > iright:
                slopeneighbours = outerneighbours[ileft:iright-1:-1]
            else:
                slopeneighbours = outerneighbours[iright:]+outerneighbours[:ileft+1]
            distance = [(pt[0]-i)**2+(pt[1]-j)**2 for pt in slopeneighbours]
            mindistance = min(distance)
            argdist = distance.index(mindistance)
            if pright == pleft:
                argdist = len(distance) - 1 - distance[::-1].index(mindistance)
            firstpoint = slopeneighbours[argdist]
            polylist = [firstpoint, point]
            #print('create new ridge', ipt, t, tright, slopeneighbours[maxpos])

            # if the ridge overlaps a dangling thalweg, the side of the thalweg must be recorded
            if pright == polylist[1] and -t in thalwegline:
                divideside[self.linekey] = -1
            if pleft == polylist[1] and -tleft in thalwegline:
                divideside[self.linekey] = 1

            self.addNewRidge(self.linekey, ihill, ipt, t, polylist, CellType.DIVIDE)
            self.linekey += 1
        
        for ir in divideside:
            self.linedict[ir]['startside'] = divideside[ir]

    def correctDivideToRidge(self):
        ridge = {ir for ir in self.linedict if self.linedict[ir]['type'] == CellType.RIDGE}
        for ir in ridge:
            r = self.linedict[ir]
            i = r['end']
            ihill = r['hill']
            lr = self.getDividesAtNode(i)
            lr = [j for j in lr if j > 0] # get all the divides that leave i
            lr = [j for j in lr if self.linedict[j]['hill'] == ihill] # keep those that are in the same hill
            lr = [j for j in lr if self.linedict[j]['type'] == CellType.DIVIDE]
            while lr:
                r = self.linedict[lr[0]]
                r['type'] = CellType.RIDGE
                i = r['end']
                lr = self.getDividesAtNode(i)
                lr = [j for j in lr if j > 0] # get all the divides that leave i
                lr = [j for j in lr if self.linedict[j]['hill'] == ihill] # keep those that are in the same hill
                lr = [j for j in lr if self.linedict[j]['type'] == CellType.DIVIDE]
                
    # def orderRidgesAroundVirtualPit(self):
    #     pt = self.nodedict[-1]
    #     ridgelist = pt['ridge']
    #     if len(ridgelist) < 2:
    #         return
    #     boundary = self.terrain.boundary
    #     inbound = self.terrain.inbound
    #     try:
    #         boundaryposition = {ir: boundary.index(self.ridgedict[ir]['polyline'][0]) for ir in ridgelist}
    #         inboundposition = {ir: inbound.index(self.ridgedict[ir]['polyline'][1]) for ir in ridgelist}
    #         outlist = sorted(ridgelist, key = lambda ir: boundaryposition[ir])
    #         inlist = sorted(ridgelist, key = lambda ir: inboundposition[ir])
    #         mismatched = [(outlist[i],inlist[i]) for i in range(len(outlist)) if outlist[i]!=inlist[i]]
    #         mismatched = {(min(i),max(i)) for i in mismatched}
    #         for position in mismatched:
    #             ir0 = position[0]
    #             ir1 = position[1]
    #             indice0in = [i for i, p in enumerate(inbound) if p == self.ridgedict[ir0]['polyline'][1]]
    #             indice1in = [i for i, p in enumerate(inbound) if p == self.ridgedict[ir1]['polyline'][1]]
    #             okindice0in = len(indice0in) == 1
    #             okindice1in = len(indice1in) == 1
    #             if not okindice0in or not okindice1in:
    #                 indice0out = [i for i, p in enumerate(boundary) if p == self.ridgedict[ir0]['polyline'][0]]
    #                 indice1out = [i for i, p in enumerate(boundary) if p == self.ridgedict[ir1]['polyline'][0]]
    #                 okindice0out = len(indice0out) == 1
    #                 okindice1out = len(indice1out) == 1
    #                 if okindice0out and okindice1out:
    #                     tmp = inboundposition[ir0]
    #                     inboundposition[ir0] = inboundposition[ir1]
    #                     inboundposition[ir1] = tmp
    #                 else:
    #                     if okindice0in and okindice1out:
    #                         i0 = outlist.index(ir0)
    #                         i1 = inlist.index(ir1)
    #                         if i0 == i1:
    #                             tmp = inboundposition[ir0]
    #                             inboundposition[ir0] = inboundposition[ir1]
    #                             inboundposition[ir1] = tmp
    #                         else:
    #                             print('x', position)
    #                     elif okindice1in and okindice0out:
    #                         i0 = inlist.index(ir0)
    #                         i1 = outlist.index(ir1)
    #                         if i0 == i1:
    #                             tmp = inboundposition[ir0]
    #                             inboundposition[ir0] = inboundposition[ir1]
    #                             inboundposition[ir1] = tmp
    #                         else:
    #                             print('+', position)
    #                     else:    
    #                         print(ir0, ir1, indice0in, indice1in, indice0out, indice1out)
    #         inlist = sorted(ridgelist, key = lambda ir: inboundposition[ir])
    #         mismatched = [(outlist[i],inlist[i]) for i in range(len(outlist)) if outlist[i]!=inlist[i]]
    #         mismatched = {(min(i),max(i)) for i in mismatched}
    #         ridgelist = inlist
    #     except:
    #         print('Exception orderRidges around virtual pit')
    #     self.nodedict[-1]['ridge'] = ridgelist
        
    # def orderRidgesAroundNodes(self):
    #     """
    #     Orders the list of ridges around each node in a clockwise order 
        
    #     To be done: holes are not considered.
    
    #     Returns
    #     -------
    #     None.
    
    #     """
    #     for ipt, pt in self.nodedict.items():
    #         # -1 is the id of the virtual pit closing outside the domain
    #         ridgelist = pt['ridge']
    #         nr = len(ridgelist)
    #         if nr < 2:
    #             continue
    #         if ipt == -1:
    #             self.orderRidgesAroundVirtualPit()
    #         elif ipt < -1:
    #             # holes must be virtual pits with id<-1
    #             # get all the neighbours around the saddles
    #             continue
    #         else:    
    #             (i,j) = pt['ij']
    #             neighbours = [(i+di, j+dj) for (di,dj) in df.ldr8]
    #             # arrange the lines in the same order as the neighbours
    #             try:
    #                 ridgelist.sort(key = lambda it: neighbours.index(self.ridgedict[abs(it)]['polyline'][df.getSecondPointPosition(it)]))
    #             except (IndexError, ValueError) as err:
    #                 print('Exception raised in orderRidgesAroundNodes', err, ipt, pt, neighbours)
    #             # check if two ridges overlap at the node
    #             for i in range(nr):
    #                 ir1 = ridgelist[i-1]
    #                 ir2 = ridgelist[i]
    #                 absir1 = abs(ir1)
    #                 absir2 = abs(ir2)
    #                 hill1 = self.ridgedict[absir1]['hill']
    #                 hill2 = self.ridgedict[absir2]['hill']
    #                 secondpoint = self.ridgedict[absir1]['polyline'][df.getSecondPointPosition(ir1)]
    #                 if secondpoint == self.ridgedict[absir2]['polyline'][df.getSecondPointPosition(ir2)]:
    #                     thalweglist = pt['thalweg']
    #                     if hill1 == hill2:
    #                         for itt in thalweglist:
    #                             t = self.thalwegdict[abs(itt)]
    #                             if t['lefthill'] == hill1 and t['righthill'] == hill1 and secondpoint == t['polyline'][df.getSecondPointPosition(itt)]:
    #                                 it = itt
    #                                 break
    #                         try:
    #                             if ir1 > 0:
    #                                 ts1 = self.ridgedict[absir1]['startside']
    #                             else:
    #                                 ts1 = self.ridgedict[absir1]['endside']
    #                             if ir2 > 0:
    #                                 ts2 = self.ridgedict[absir2]['startside']
    #                             else:
    #                                 ts2 = self.ridgedict[absir2]['endside']
    #                         except (IndexError, ValueError, KeyError) as err:
    #                             print('Exception raised in orderRidgesAroundNodes', err)
    #                             print('if hill1 == hill2', ipt, it, ir1, ir2)
    #                         if ts1 == ts2:
    #                             print("Error at node", ipt, "ridges", ir1, "and", ir2, "on the same side")
    #                         # the one on the left should have a thalwegside of -1
    #                         # if not, we swap them
    #                         if (ts1 == 1 and it > 0) or (ts1 == -1 and it < 0):
    #                             ridgelist[i] = ir1
    #                             ridgelist[i-1] = ir2
    #                     else:
    #                         for it in thalweglist:
    #                             if it < 0:
    #                                 t = self.thalwegdict[-it]
    #                                 if t['polyline'][-2] != secondpoint:
    #                                     continue
    #                                 lefthill = t['righthill']
    #                                 righthill = t['lefthill']
    #                             else:
    #                                 t = self.thalwegdict[it]
    #                                 if t['polyline'][1] != secondpoint:
    #                                     continue
    #                                 lefthill = t['lefthill']
    #                                 righthill = t['righthill']
    #                             if hill2 == lefthill and hill1 == righthill:
    #                                 ridgelist[i-1] = ir2
    #                                 ridgelist[i] = ir1
    
    def computeThalwegDales(self):
        """
        Computes the dictionary of dales of the surface network. A dale is 
        formed by a series of ridges. The dale can contain some holes.
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
        # initialise left and right dales in ridges
        for ir in self.linedict:
            r = self.linedict[ir]
            if r['type'] & CellType.DIVIDE_FLAG:
                r['leftdale'] = None
                r['rightdale'] = None
    
        self.daledict = {}
    
        # compute a dale around each downslope line
#        thalweglist = [it for it in self.linedict if self.linedict[it]['type'] & CellType.THALWEG_FLAG]
        thalweglist = [it for it in self.linedict if self.isDownslopeLine(self.linedict[it]['type'])]
        for it in thalweglist:
            self.computeDale(it)
            
    def computeDale(self, it):
            t = self.linedict[it]
            #if t['end'] == -1: # no dale for thalwegs ending outside
            #    t['dale'] = -1
            #    continue
            ipt = t['start'] # we start from the start node of the thalweg
            pt = self.nodedict[ipt]
    
            # we will start from the ridge on the left
            tpos = pt['line'].index(it)
            rpos = self.previousDivideIndexAtNode(ipt, tpos)
#            print(ipt, tpos, rpos)
            ir = pt['line'][rpos]
            
            # we have the first ridge of the dale, we can turn around
            ridgering, ridgeline = self.getDivideRing(ipt, ir)
            # if there are several rings, the boundary is the biggest one        
            # the others are holes
            border = ridgering[0]
            holes = []
            arealist = []
            #print(it, ir, "There is a hole", ridgering, ridgeline)
            # we compute the area of each ring
            for ring in ridgering:
                pointlist = []
                #for rr in ring:
                rr = ring[-1]
                for i in range(len(ring)):
                    lastrr = rr
                    rr = ring[i]
                    if rr > 0:
                        r = self.linedict[rr]
                        if r['start'] == -1:
                            firstpoint = self.linedict[-lastrr]['polyline'][0]
                            lastpoint = r['polyline'][0]
                            istart = self.terrain.boundary.index(firstpoint)
                            iend = self.terrain.boundary.index(lastpoint)
                            if istart > iend:
                                pointlist += self.terrain.boundary[istart:iend:-1]
                            elif iend > istart:
                                pointlist += self.terrain.boundary[istart::-1]
                                pointlist += self.terrain.boundary[-1:iend:-1]
                        pointlist += r['polyline'][:-1]
                    else:
                        pointlist += reversed(self.linedict[-rr]['polyline'][1:])
                crosslist = []
                q = pointlist[-1]
                for p in pointlist:
                    crosslist.append(p[0]*q[1]-p[1]*q[0])
                    q = p
                area = sum(crosslist)
                arealist.append(area)
            iborder = arealist.index(max(arealist))
            border = ridgering.pop(iborder)
            holes = ridgering
            
            area = sum(arealist)
            if area < 0:
                print('dale', it, "is negative")
                area = -area
            dalekey = it
            self.daledict[dalekey] = {'boundary': border, 
                                      'hole': holes,
                                      'dangling': ridgeline,
                                      'area': 0.5*area}
            t['dale'] = dalekey
            for rr in border:
                if rr > 0:
                    self.linedict[rr]['rightdale'] = dalekey
                else:
                    self.linedict[-rr]['leftdale'] = dalekey
            for r in holes:
                for rr in r:
                    if rr > 0:
                        self.linedict[rr]['rightdale'] = dalekey
                    else:
                        self.linedict[-rr]['leftdale'] = dalekey
            for rr in ridgeline:
                self.linedict[rr]['leftdale'] = dalekey
                self.linedict[rr]['rightdale'] = dalekey

    def computeDale2(self, it):
            t = self.linedict[it]
            #if t['end'] == -1: # no dale for thalwegs ending outside
            #    t['dale'] = -1
            #    continue
            ipt = t['start'] # we start from the start node of the thalweg
            pt = self.nodedict[ipt]
            if it == 21:
                print("in computeDale", it, ipt)
    
            # we will start from the ridge on the left
            lefthill = t['lefthill']
            righthill = t['righthill']
            noderidge = self.getDividesAtNode(ipt)
            ridgelist = [ir for ir in noderidge if self.linedict[abs(ir)]['hill'] == lefthill]
            if it == 21:
                print("in computeDale", noderidge, ridgelist)
            # if there is more than one ridge, we have to choose the closest one on the left
            ir = ridgelist[0] # take the one ridge in the list if there's only one
            lenridgelist = len(ridgelist)
            # if more than one ridge in the hill
            if lenridgelist > 1:
                (i,j) = pt['ij']
                # take the second point of each ridge
                secondpoints = [self.linedict[abs(jr)]['polyline'][df.getSecondPointPosition(jr)] for jr in ridgelist]
                thalwegpoint = t['polyline'][1]
                # if the thalweg overlaps a ridge
                if it == 21:
                    print("in computeDale test", thalwegpoint in secondpoints)
                if thalwegpoint in secondpoints:
                    # we take the overlapping ridge
                    k = secondpoints.index(thalwegpoint)
                    jr = ridgelist[k]
                    # if it is not a dangling thalweg, it is the closest we can get
                    if lefthill != righthill:
                        ir = jr
                    else:
                        # if it is a dangling thalweg, it is the same hill on both sides
                        # if the ridge starts at node ipt
                        if jr > 0:
                            r = self.linedict[jr]
                            if 'startside' in r: 
                                if r['startside'] == 1: # on the right of the thalweg
                                    ir = ridgelist[k-1] # take the previous one
                                else: # on the left
                                    ir = jr # take this one
                        # if the ridge ends at node ipt
                        if jr < 0:
                            r = self.linedict[-jr]
                            if 'endside' in r:
                                if r['endside'] == 1: # on the right of the thalweg
                                    ir = ridgelist[k-1] # take the previous one
                                else: # on the left
                                    ir = jr # take this one
                # if the thalweg is between two ridges
                else:
                    # get the coordinates of neighbouring pixels
                    neighbours = [(i+di, j+dj) for (di,dj) in df.ldr8]
                    jt = neighbours.index(thalwegpoint)
                    found = False
                    for k in range(lenridgelist): # lenridgelist = len(secondpoints)
                        ptr1 = secondpoints[k-1]
                        ptr2 = secondpoints[k]
                        if neighbours.index(ptr1) < jt and neighbours.index(ptr2) > jt:
                            # the talweg is between the two ridges
                            ir = ridgelist[k-1] # we take the first one
                            found = True
                            break
                    if it == 21:
                        print("in computeDale test", neighbours, jt, found)
                    if not found:
                        ir = ridgelist[-1]
            
            # we have the first ridge of the dale, we can turn around
            ridgering, ridgeline = self.getDivideRing(ipt, ir)
            if it == 21:
                print("in computeDale", ir, ridgering, ridgeline)
            # if there are several rings, the boundary is the biggest one        
            # the others are holes
            border = ridgering[0]
            holes = []
            arealist = []
            #print(it, ir, "There is a hole", ridgering, ridgeline)
            # we compute the area of each ring
            for ring in ridgering:
                pointlist = []
                #for rr in ring:
                rr = ring[-1]
                for i in range(len(ring)):
                    lastrr = rr
                    rr = ring[i]
                    if rr > 0:
                        r = self.linedict[rr]
                        if r['start'] == -1:
                            firstpoint = self.linedict[-lastrr]['polyline'][0]
                            lastpoint = r['polyline'][0]
                            istart = self.terrain.boundary.index(firstpoint)
                            iend = self.terrain.boundary.index(lastpoint)
                            if istart > iend:
                                pointlist += self.terrain.boundary[istart:iend:-1]
                            elif iend > istart:
                                pointlist += self.terrain.boundary[istart::-1]
                                pointlist += self.terrain.boundary[-1:iend:-1]
                        pointlist += r['polyline'][:-1]
                    else:
                        pointlist += reversed(self.linedict[-rr]['polyline'][1:])
                crosslist = []
                q = pointlist[-1]
                for p in pointlist:
                    crosslist.append(p[0]*q[1]-p[1]*q[0])
                    q = p
                area = sum(crosslist)
                arealist.append(area)
            iborder = arealist.index(max(arealist))
            border = ridgering.pop(iborder)
            holes = ridgering
            
            area = sum(arealist)
            if area < 0:
                print('dale', it, "is negative")
                area = -area
            dalekey = it
            self.daledict[dalekey] = {'boundary': border, 
                                      'hole': holes,
                                      'dangling': ridgeline,
                                      'area': 0.5*area}
            t['dale'] = dalekey
            for rr in border:
                if rr > 0:
                    self.linedict[rr]['rightdale'] = dalekey
                else:
                    self.linedict[-rr]['leftdale'] = dalekey
            for r in holes:
                for rr in r:
                    if rr > 0:
                        self.linedict[rr]['rightdale'] = dalekey
                    else:
                        self.linedict[-rr]['leftdale'] = dalekey
            for rr in ridgeline:
                self.linedict[rr]['leftdale'] = dalekey
                self.linedict[rr]['rightdale'] = dalekey

    def getDivideRing(self, idnode, idridge):
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
            r = self.linedict[currentridge]
            if currentnode == r['start']:
                nextnode = r['end']
                ridgering.append(currentridge)
            elif currentnode == r['end']:
                nextnode = r['start']
                ridgering.append(-currentridge)
            else:
                print('error in getDivideRing: ridge', currentridge, 'and node', currentnode, 'do not match')

            currentnode = nextnode
            ridgelist = [abs(r) for r in self.getDividesAtNode(currentnode)]
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
    
    def addThalwegSlopeInDegree(self):
        for it, t in self.linedict.items():
            if not (t['type'] & CellType.THALWEG_FLAG):
                continue
            ptstart = self.nodedict[t['start']]
            if t['end'] == -1:
                if len(t['polyline']) <= 2:
                    t['slope'] = 0
                    continue
                length = df.length(t['polyline'][:-1]) * self.terrain.x_size
                pt = t['polyline'][-2]
                z = self.terrain.dtm[pt[0],pt[1]]
                dz = ptstart['z'] - z
                slope = dz/length
                t['slope'] = 180*atan(slope)/pi
                continue
            # we assume that x_size = y_size, otherwise, rewrite the length function
            length = df.length(t['polyline']) * self.terrain.x_size
            ptend = self.nodedict[t['end']]
            dz = ptstart['z'] - ptend['z']
            slope = dz/length
            slope = 180*atan(slope)/pi
            t['slope'] = slope
    
    def addDivideSlopeInDegree(self):
        for ir, r in self.linedict.items():
            if not (r['type'] & CellType.DIVIDE_FLAG):
                continue
            if r['start'] == -1:
                r['slope'] = 90
                continue
            # we assume that x_size = y_size, otherwise, rewrite the length function
            length = df.length(r['polyline']) * self.terrain.x_size
            ptend = self.nodedict[r['end']]
            ptstart = self.nodedict[r['start']]
            dz = ptend['z'] - ptstart['z']
            slope = dz/length
            slope = 180*atan(slope)/pi
            r['slope'] = slope

    def getRidgeMinZ(self, iridge):
        """
        Return the z of the start point of the ridge, which is the lowest point

        Parameters
        ----------
        iridge : integer
            Ridge id.

        Returns
        -------
        Float
            The z value of the ridge start node.

        """
        return self.nodedict[self.linedict[iridge]['start']]['z']
    
    def getRidgeMaxZ(self, iridge):
        """
        Return the z of the end point of the ridge, which is the highest point

        Parameters
        ----------
        iridge : integer
            Ridge id.

        Returns
        -------
        Float
            The z value of the ridge end node.

        """
        return self.nodedict[self.linedict[iridge]['end']]['z']
    
    def getHillHeight(self, idhill):
        """
        returns the height of a hill measured between its peak and the lowest
        point on its boundary

        Parameters
        ----------
        idhill : integer
            the id of the hill.

        Returns
        -------
        float
            the height of the hill.

        """
        hill = self.hilldict[idhill]
        peak = self.getPeak(idhill)
        zmax = self.nodedict[peak]['z']
        thalwegs = hill['boundary']
        nodeset = set()
        for it in thalwegs:
            inode = self.linedict[abs(it)]['end']
            nodeset.add(inode)
        # get the lowest node
        zmin = self.terrain.maxz
        for inode in nodeset:
            z = self.nodedict[inode]['z']
            if z < zmin:
                zmin = z
        return zmax - zmin

    def getNeighbourHills(self, idhill):
        """
        Returns a set of hills neighbouring the hill. The set is not ordered
        and each hill appears only once. The set does not include hills within
        that form a hole.

        Parameters
        ----------
        idhill : integer
            id of the hill.

        Returns
        -------
        neighbours : set
            the ids of the neighbouring hills.

        """
        hill = self.hilldict[idhill]
        neighbours = set()
        boundary = hill['boundary']
        for it in boundary:
            if it > 0:
                neighbours.add(self.thalwegdict[it]['lefthill'])
            else:
                neighbours.add(self.thalwegdict[-it]['righthill'])
        return neighbours
    
    def getInnerHills(self, idhill):
        """
        Returns a set of hills located within the hill. The set is not ordered
        and each hill appears only once. These hills are holes in the hill.

        Parameters
        ----------
        idhill : integer
            id of the hill.

        Returns
        -------
        neighbours : set
            the ids of the hills within.

        """
        hill = self.hilldict[idhill]
        neighbours = set()
        holes = hill['hole']
        for it in holes:
            if it > 0:
                neighbours.add(self.thalwegdict[it]['lefthill'])
            else:
                neighbours.add(self.thalwegdict[it]['righthill'])
        return neighbours
    
    def getPathToPeak(self, idridge):
        """
        Computes a path along the ridges up to a peak starting at a given ridge.
        The path is given by a succession of ridges and goes upward to the peak
        located in the same hill as the initial ridge.

        Parameters
        ----------
        idridge : integer
            ridge from which the path starts.

        Returns
        -------
        path : list of integers
            list of ridge ids forming a path to a peak.

        """
        r = self.ridgedict[idridge]
        jnode = r['end']
        jhill = r['hill']
        node = self.nodedict[jnode]
        path = [idridge]
        while node['type'] != df.PEAK:
            rlist = [ir for ir in node['ridge'] if ir > 0 and self.ridgedict[ir]['hill'] == jhill]
            # if there are several ridges, take the steepest one
            kmax = rlist[0]
            maxslope = -1
            for kr in rlist:
                slope = self.ridgedict[kr]['slope']
                if slope > maxslope:
                    maxslope = slope
                    kmax = kr
            jr = kmax
            jnode = self.ridgedict[jr]['end']
            node = self.nodedict[jnode]
            path.append(jr)
        return path
    
    def computePathToPeak(self, inode, ihill):
        """
        Compute the ascending path from a node to a peak inside a given hill.

        Parameters
        ----------
        inode : integer
            index of the starting node.
        ihill : integer
            index of the hill.

        Returns
        -------
        path : list of integers
            list of ridge ids forming a path to a peak.

        """
        node = self.nodedict[inode]
        ridges = [ir for ir in node['ridge'] if ir>0]
        maxslope = -1
        kmax = ridges[0]
        for kr in ridges:
            r = self.ridgedict[kr]
            if r['hill'] != ihill:
                continue
            slope = r['slope']
            if slope > maxslope:
                maxslope = slope
                kmax = kr
        path = self.getPathToPeak(kmax)
        return path

    def getPassBetweenHills(self, ihill, jhill):
        """
        For two adjacent hills, return the highest node located on the boundary
        of the two hills. This node is located on the path connecting the two peaks

        Parameters
        ----------
        ihill : integer
            index of the first hill.
        jhill : integer
            index of the second hill.

        Returns
        -------
        ipass : integer
            index of the pass between the two hills.

        """
        # get the thalwegs between the two hills
        ibound = self.hilldict[ihill]['boundary']
        jbound = self.hilldict[jhill]['boundary']
        boundary = [i for i in ibound if -i in jbound]
        # get the nodes along this boundary
        nodes = set()
        for i in boundary:
            t = self.thalwegdict[abs(i)]
            nodes.add(t['start'])
            nodes.add(t['end'])
        # get the highest node of them
        zmax = self.terrain.nodata
        ipass = -1
        for ip in nodes:
            z = self.nodedict[ip]['z']
            if z > zmax:
                zmax = z
                ipass = ip
        return ipass

    def orderRidgesInHill(self, ihill):
        orderedridge = []
        hill = self.hilldict[ihill]
        nodeset = set()
        allthalwegs = hill['boundary'] + hill['dangling']
        for hole in hill['hole']:
            allthalwegs.append(hole)
        for it in allthalwegs:
            t = self.thalwegdict[abs(it)]
            nodeset.add(t['start'])
            nodeset.add(t['end'])
        ipeak = self.getPeak(ihill)
        peak = self.nodedict[ipeak]
        stack = peak['ridge'][:]
        while stack:
            ir = abs(stack.pop())
            orderedridge.append(ir)
            inode = self.ridgedict[ir]['start']
            if inode in nodeset:
                continue
            ridgelist = self.nodedict[inode]['ridge'][:]
            i = ridgelist.index(ir)
            if i != len(ridgelist)-1:
                ridgelist = ridgelist[i+1:]+ridgelist[:i]
            else:
                ridgelist = ridgelist[:-1]
            stack.extend(ridgelist)
        return orderedridge
    
    def getHillChain(self, ihill):
        # get all the nodes on the border of the hill
        hill = self.hilldict[ihill]
        border = hill['boundary']
        thalweg = [abs(it) for it in border]
        node = {self.thalwegdict[it]['start'] for it in thalweg}
        node.update({self.thalwegdict[it]['end'] for it in thalweg})
#        if -1 in node:
#            node.remove(-1)
        z = [self.nodedict[s]['z'] for s in node]
        zmin = min(z)
        zmax = max(z)
        ipeak = self.getPeak(ihill)
        zpeak = self.nodedict[ipeak]['z']
#        if zpeak - zmin < df.precision:
#            return None
        saddle = {p for p in node if self.nodedict[p]['type'] == df.SADDLE}
        saddle = {p for p in saddle if self.isSpring(p)}
        path = {}
        for isaddle in saddle:
            ridges = [ir for ir in self.nodedict[isaddle]['ridge'] if ir>0]
            ridges = [ir for ir in ridges if self.ridgedict[ir]['hill'] == ihill]
            if len(ridges) != 1:
                print("getHillChain, at saddle", isaddle, "ridges are", ridges)
            path[isaddle] = self.getPathToPeak(ridges[0])
#            path[isaddle] = (zpeak - self.nodedict[isaddle]['z']) / (zpeak - zmin)
        #saddle = sorted(saddle, key = lambda s: self.nodedict[s]['z'], reverse = True)
        
        return path
        # one path to keep is the path to the keysaddle
        
        # if len(path) > 2:
        #     ridgeset = set()
        #     for p in path:
        #         ridgeset.update(p)
        #     ipeak = self.getPeak(ihill)
        #     peak = self.nodedict[ipeak]
        #     ridge = [-ir for ir in peak['ridge']]
        #     ridge = [ir for ir in ridge if ir in ridgeset]
        #     ridge = sorted(ridge, key = lambda ir: self.ridgedict[ir]['slope'])[:2]
        #     chains = []
        #     while ridge:
        #         ir = -ridge.pop()
        #         chain = [ir]
        #         r = self.ridgedict[ir]
        #         p = self.nodedict[r['start']]
        #         while p['type'] == df.JUNCTION:
        #             rlist = [-i for i in p['ridge'] if i < 0]
        #             rlist = [i for i in rlist if i in ridgeset]
        #             ir = min(rlist, key = lambda i: self.ridgedict[i]['slope'])
        #             chain.append(ir)
        #             r = self.ridgedict[ir]
        #             p = self.nodedict[r['start']]
        #         chains.append(chain)
        # ridgechain = chains[0][::-1]
        # if len(chains) > 1:
        #     ridgechain.extend(chains[1])
        # return ridgechain
            
    def getDalesInBasin(self, ithalweg):
        t = self.thalwegdict[ithalweg]
        if t['dale'] < 1:
            return None
        # collect all the dales above the thalweg
        visited = {} # mark if a thalweg has been visited or not
        # this is necessary if distributaries are allowed
        for it in self.thalwegdict:
            visited[it] = False # mark all thalwegs as unvisited
        dalelist = []
        stack = [ithalweg]
        # the loop visits all thalwegs upstream of ithalweg to collect all the
        # dales that form the basin
        while stack:
            it = stack.pop()
            if visited[it]:
                continue
            visited[it] = True
            dalelist.append(self.thalwegdict[it]['dale'])
            if 'headwater' in self.thalwegdict[it]:
                continue
            inode = self.fromNode(it)
            thalweglist = self.upstreamThalwegs(inode)
            stack += thalweglist
        # dalelist contains all the dales composing the drainage basin
        return dalelist

    def getRidgesInBasin(self, ithalweg):
        """
        Returns all the ridges found in a basin. Ridge id are all positive.

        Parameters
        ----------
        ithalweg : integer
            Thalweg at the bottom of the basin.

        Returns
        -------
        Set of integer
            Ridges found in the basin. No duplicates and all id are positive.

        """
        t = self.thalwegdict[ithalweg]
        if t['dale'] < 1:
            return None
        # collect all the dales above the thalweg
        visited = {} # mark if a thalweg has been visited or not
        # this is necessary if distributaries are allowed
        for it in self.thalwegdict:
            visited[it] = False # mark all thalwegs as unvisited
        dalelist = []
        stack = [ithalweg]
        # the loop visits all thalwegs upstream of ithalweg to collect all the
        # dales that form the basin
        while stack:
            it = stack.pop()
            if visited[it]:
                continue
            visited[it] = True
            dalelist.append(self.thalwegdict[it]['dale'])
            if 'headwater' in self.thalwegdict[it]:
                continue
            inode = self.fromNode(it)
            thalweglist = self.upstreamThalwegs(inode)
            stack += thalweglist
        # dalelist contains all the dales composing the drainage basin
        
        # we need to find the ridges on the border of the basin
        # we take the boundary of each dale
        boundarylist = [self.daledict[idale]['boundary'] for idale in dalelist]
        # each boundary is a list of ridges
        # we take all the holes in the dale, each hole is a list of rings
        holelist = [self.daledict[idale]['hole'] for idale in dalelist]
        # we turn the holes into a list of rings
        holelist = [iring for sublist in holelist for iring in sublist]
        # now we can turn everything into sets of ridges
        outerridgeset = {abs(iridge) for sublist in boundarylist for iridge in sublist}
        innerridgeset = {abs(iridge) for sublist in holelist for iridge in sublist}
        return outerridgeset.union(innerridgeset)

    def getDivideSetBoundingThalwegs(self, thalwegset):
        daleset = {self.linedict[it]['dale'] for it in thalwegset}
        # we need to find the ridges on the border of the basin
        # we take the boundary of each dale
        boundarylist = [self.daledict[idale]['boundary'] for idale in daleset]
        # each boundary is a list of ridges
        # we take all the holes in the dale, each hole is a list of rings
        holelist = [self.daledict[idale]['hole'] for idale in daleset]
        # we turn the holes into a list of rings
        holelist = [iring for sublist in holelist for iring in sublist]
        # now we can turn everything into lists of ridges
        outerridgeset = {iridge for sublist in boundarylist for iridge in sublist}
        innerridgeset = {iridge for sublist in holelist for iridge in sublist}
        
        # remove all ridges of innerridgeset from outerridgeset
        # ridges should have opposite signs
        # inner ridges are also considered
        ridgeset = outerridgeset.union(innerridgeset)
        # remove from ridgelist all ridges that appear twice with different sign
        plusset = {iridge for iridge in ridgeset if iridge > 0}
        minusset = {iridge for iridge in ridgeset if iridge < 0}
        plusridgeset = {iridge for iridge in plusset if -iridge not in minusset}
        minusridgeset = {iridge for iridge in minusset if -iridge not in plusset}
        return plusridgeset, minusridgeset
        #ridgeset = plusridgeset | minusridgeset
    
    def getRidgeSetBoundingThalwegs2(self, thalwegset):
        daleset = {self.thalwegdict[it]['dale'] for it in thalwegset}
        # we need to find the ridges on the border of the basin
        # we take the boundary of each dale
        boundarylist = [self.daledict[idale]['boundary'] for idale in daleset]
        # each boundary is a list of ridges
        # we take all the holes in the dale, each hole is a list of rings
        holelist = [self.daledict[idale]['hole'] for idale in daleset]
        # we turn the holes into a list of rings
        holelist = [iring for sublist in holelist for iring in sublist]
        # now we can turn everything into lists of ridges
        outerridgeset = {iridge for sublist in boundarylist for iridge in sublist}
        innerridgeset = {iridge for sublist in holelist for iridge in sublist}
        
        # remove all ridges of innerridgeset from outerridgeset
        # ridges should have opposite signs
        # all inner ridges should be removed unless there are some endorheic bassins
        ridgeset = {iridge for iridge in outerridgeset if -iridge not in innerridgeset}
        # remove from ridgelist all ridges that appear twice with different sign
        plusset = {iridge for iridge in ridgeset if iridge > 0}
        minusset = {iridge for iridge in ridgeset if iridge < 0}
        plusridgeset = {iridge for iridge in plusset if -iridge not in minusset}
        minusridgeset = {iridge for iridge in minusset if -iridge not in plusset}
        return plusridgeset, minusridgeset
        #ridgeset = plusridgeset | minusridgeset

    def turnRidgeSetToRidgeRing(self, ridgeset):
        # ridgeset contains all the ridges on the basin boundary
        
        # take any ridge from ridgeset
        for ir in ridgeset:
            break
        depression = [ir]
        iridge = ir
        startnode = None
        endnode = None
        if iridge > 0:
            ridge = self.ridgedict[iridge]
            startnode = ridge['start']
            endnode = ridge['end']
        else:
            ridge = self.ridgedict[-iridge]
            startnode = ridge['end']
            endnode = ridge['start']
        currentnode = endnode
        while currentnode != startnode:
            ridgeset.remove(iridge)
            node = self.nodedict[currentnode]
            ridges = set(node['ridge']) & ridgeset
            if len(ridges) > 1:
                print("problem, two ridges at node", currentnode, iridge, ridges)
                print(self.nodedict[currentnode]['ridge'])
                ihill = self.ridgedict[abs(iridge)]['hill']
                for lr in ridges:
                    if self.ridgedict[abs(lr)]['hill'] == ihill:
                        iridge = lr
                        break
                #index = node['ridge'].index(-iridge) + 1
                #if index >= len(node['ridge']):
                #    index = 0
                #iridge = node['ridge'][index]
                if not iridge in ridges:
                    print("problem finding the next ridge", currentnode, iridge, startnode)
            else:
                iridge = ridges.pop()
            depression.append(iridge)
            if iridge > 0:
                ridge = self.ridgedict[iridge]
                currentnode = ridge['end']
            else:
                ridge = self.ridgedict[-iridge]
                currentnode = ridge['start']
        return depression

    def ridgeRingToPolygon(self, ridgering):
        pring = []
        # get the list of points forming the ring
        for ir in ridgering:
            if ir > 0:
                r = self.ridgedict[ir]
                if pring:
                    pring.extend(r['polyline'][1:])
                else:
                    pring.extend(r['polyline'][:])
            else:
                r = self.ridgedict[-ir]
                if r['start'] == -1:
                    print('ir', ir, pring, "we don't deal with this case yet")
                tmp = []
                if pring:
                    tmp = r['polyline'][:-1]
                else:
                    tmp = r['polyline'][:]
                tmp.reverse()
                pring.extend(tmp)
        plgn = shapely.geometry.Polygon(pring)
        return plgn

    def fillPuddle(self, idpuddle):
        terrain = self.terrain
        # get the ridge list bounding the depression
        thalwegset = self.getThalwegsInPuddle(idpuddle)
        plusridgeset, minusridgeset = self.getRidgeSetBoundingThalwegs(thalwegset)
        depression = self.turnRidgeSetToRidgeRing(plusridgeset|minusridgeset)
        plgn = self.ridgeRingToPolygon(depression)
        shapely.prepare(plgn)
        
        outlet = self.puddledict[idpuddle]['outlet']
        zmax = self.nodedict[outlet]['z']
        pit = self.nodedict[idpuddle]
        (i,j) = pit['ij']
        stack = {(i,j)}
        visited = set()
        water = []
        while stack:
            (i,j) = stack.pop()
            if (i,j) not in visited and terrain.isInside(i,j) and terrain.dtm[i,j] <= zmax and shapely.intersects_xy(plgn, i, j):
                water.append((i,j))
                neighbours = [(i+di, j+dj) for (di, dj) in df.ldr8]
                stack.update(neighbours)
            visited.add((i,j))
        return water

        # create a collection of points
#        pond = shapely.geometry.MultiPoint(water)
#        pond = shapely.buffer(pond, 0.5)
#        return pond
        
    def createClassMatrices(self):
        terrain = self.terrain
        # rebuild the nodeidx without the virtual node
        self.nodeidx = {node['ij']: key for key, node in self.nodedict.items() if key >= 0}
        # build the thalweg index
        # first and last nodes of thalwegs are not recorded since there are several
        # thalwegs and there is a node connecting the thalwegs
        self.thalwegidx = {}
        for key, thalweg in self.thalwegdict.items():
            line = thalweg['polyline']
            for ij in line[1:-1]:
                self.thalwegidx[ij] = key
        # build the ridge index
        # works like the thalweg index except that two ridges can overlap so we have to store pairs
        self.ridgeidx = {}
        for key, ridge in self.ridgedict.items():
            line = ridge['polyline']
            for ij in line[1:-1]:
                if ij in self.ridgeidx:
                    self.ridgeidx[ij].append(key)
                else:
                    self.ridgeidx[ij] = [key]
        # build the hill matrix
        # done by taking each hill and finding all the pixels strictly inside the hill
        boundary = set(terrain.boundary)
        print("compute hill matrix")
        self.hillmatrix = np.zeros((terrain.m, terrain.n), dtype=int)
        for key, hill in self.hilldict.items():
            # get all the thalwegs bounding the hill, including the holes
            # to be checked: add the dangling thalwegs or not
            thalwegs = hill['boundary'][:]
            for h in hill['hole']:
                thalwegs.extend(h)
            thalwegs = [abs(t) for t in thalwegs]
            tpointset = set()
            # get all the thalweg points
            for t in thalwegs:
                tpointset.update(self.thalwegdict[t]['polyline'])
            # get all the ridge points
            ridges = hill['ridge']
            rpointset = set()
            for r in ridges:
                line = self.ridgedict[r]['polyline'][1:]
                rpointset.update(line)
            # remove all points of ridges that are also on thalwegs
            rpointset.difference_update(tpointset)
            stack = set(rpointset)
            while stack:
                i,j = stack.pop()
                neighbours = [(i+di, j+dj) for di,dj in df.ldr4]
                neighbours = [p for p in neighbours if p not in rpointset and p not in tpointset and p not in boundary]
                stack.update(neighbours)
                rpointset.update(neighbours)
            for (i,j) in rpointset:
                self.hillmatrix[i,j] = key
        print("compute dale matrix")
        self.dalematrix = np.zeros((terrain.m, terrain.n), dtype=int)
        for key, dale in self.daledict.items():
            # get all the ridges bounding the dale, including the holes
            # to be checked: add the dangling ridges or not
            ridges = dale['boundary'][:]
            for h in dale['hole']:
                ridges.extend(h)
            ridges = [abs(r) for r in ridges]
            rpointset = set()
            # get all the ridge points
            for r in ridges:
                rpointset.update(self.ridgedict[r]['polyline'])
            # get all the thalweg points
            line = self.thalwegdict[key]['polyline'][1:-1]
            tpointset = set(line)
            # remove all points of thalwegs that are also on ridges
            tpointset.difference_update(rpointset)
            stack = set(tpointset)
            counter = 0
            while stack:
                counter += 1
                i,j = stack.pop()
                neighbours = [(i+di, j+dj) for di,dj in df.ldr4]
                neighbours = [p for p in neighbours if p not in tpointset and p not in rpointset and p not in boundary]
                stack.update(neighbours)
                tpointset.update(neighbours)
            for (i,j) in tpointset:
                self.dalematrix[i,j] = key
