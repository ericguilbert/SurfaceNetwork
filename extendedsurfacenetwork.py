# -*- coding: utf-8 -*-
"""
Created on Tue Aug 5 21:30:05 2020

@author: ERGUI19
"""
from time import perf_counter_ns
from math import sqrt, atan, pi

import deffunction as df
from surfacenetwork import SurfaceNetwork
from streamnetwork import StreamNetwork

class ExtendedSurfaceNetwork(SurfaceNetwork, StreamNetwork):
    """
    DrainageNetwork stores a drainage network defined by a set of streams and a
    set of ridges that are the drainage divides. It inherits from SurfaceNetwork
    and StreamNetwork which produces the stream network as a directed graph and 
    computes the Strahler order of each stream.
    DrainageNetwork adds the ridges defined in a similar way as in SurfaceNetwork
    but adding those drainage divides and it defines catchment areas, allowing
    for the computation of flow accumulation at each confluence or saddle point.
    """
    
    def buildESNFromTerrain(self, terrain, sd = True, smooth = False):
        """
        Builds a drainage network on a given terrain. Starts by building the
        stream network. Follows with ridges. Then computes the catchment area of
        each stream based on the ridges.

        Parameters
        ----------
        terrain : Terrain
            Instance of the class Terrain on which the network is built.
        sd : boolean
            True if the network does not have diffluences (single direction)
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
            if 'line' in self.nodedict[ipt]:
                del self.nodedict[ipt]['line']
        self.instantiateRidgesAtPitsAndConfluences()
        self.traceRidges(divide, correction)
        self.orderRidgesAroundNodes()
        #self.removeSpuriousRidges()
        end = perf_counter_ns()
        print("Ridge computation time:", end - start)
        start = perf_counter_ns()
        self.computeThalwegDales()
        end = perf_counter_ns()
        print("Dale computation time:", end - start)
        start = perf_counter_ns()
        self.addPuddles()
        self.mergePuddlesAndAssignFlow()
        end = perf_counter_ns()
        print("Flow direction computation time:", end - start)
        
        start = perf_counter_ns()
        self.accumulateFlow()
        end = perf_counter_ns()
        print("Flow accumulation computation time:", end - start)
        
    def instantiateRidgesAtPitsAndConfluences(self):
        """
        Instantiates ridge objects between two thalwegs at pits of degree three
        or more. 
        
        Returns
        -------
        None.
        
        """
        # key for new ridges        
        ridgestart = max(self.ridgedict) + 1
        ridgekey = ridgestart
        divideside = {}

        # we check z>nodata to avoid starting a ridge in a hole 
        # we start ridges between two consecutive thalwegs from pits and confluences
        # pits with one thalweg are excluded
#        pitlist = [i for i in self.nodedict \
#                   if ((self.nodedict[i]['type'] == df.PIT and len(self.nodedict[i]['thalweg'])>1)\
#                   or self.nodedict[i]['type'] == df.CONFLUENCE) and self.nodedict[i]['z']>self.terrain.nodata]
        pitlist = [i for i in self.nodedict if (self.nodedict[i]['type'] == df.PIT and i!=-1) or \
                   self.nodedict[i]['type'] == df.CONFLUENCE]
        
        for ipt in pitlist:
            #if ipt!=testipt:
            #    continue
            pt = self.nodedict[ipt]
                    
            # we get the list of thalwegs
            thalweglist = pt['thalweg']
            #print('Node ', ipt, thalweglist)

            (i,j) = pt['ij']
            # orderedneighbours is the neighbours of (i,j) in clockwise order
            orderedneighbours = [(i+di, j+dj) for (di,dj) in self.terrain.neighbour[i,j]]
            # for each thalweg, we start a ridge on its right
            for it, t in enumerate(thalweglist):
                #print('thalweg', it, t)
                # at a pit, thalweg ids are all negative
                if t > 0:
                    ihill = self.thalwegdict[t]['righthill']
                else:
                    ihill = self.thalwegdict[-t]['lefthill']
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
                    divideside[ridgekey] = -1
                if pright == polylist[1] and -tright in thalwegline:
                    divideside[ridgekey] = 1
                    
                self.addNewRidge(ridgekey, ihill, ipt, polylist)
                ridgekey += 1

        # deal with ridges starting at the virtual pit
        ipt = -1
        pt = self.nodedict[ipt]
        thalweglist = pt['thalweg']
        # orderedneighbours is the neighbours of (i,j) in clockwise order
        innerneighbours = self.terrain.inbound
        outerneighbours = self.terrain.boundary
        # for each thalweg, we start a ridge on its right
        for it, t in enumerate(thalweglist):
            #print('thalweg', it, t)
            # at a pit, thalweg ids are all negative
            ihill = self.thalwegdict[-t]['lefthill']
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
            pright = self.thalwegdict[-t]['polyline'][-2]
            pleft = self.thalwegdict[-tleft]['polyline'][-2]
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
            qleft = self.thalwegdict[-tleft]['polyline'][-1]
            qright = self.thalwegdict[-t]['polyline'][-1]
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
                divideside[ridgekey] = -1
            if pleft == polylist[1] and -tleft in thalwegline:
                divideside[ridgekey] = 1

            self.addNewRidge(ridgekey, ihill, ipt, polylist)
            ridgekey += 1
        
        for ir in divideside:
            self.ridgedict[ir]['startside'] = divideside[ir]

    def orderRidgesAroundVirtualPit(self):
        pt = self.nodedict[-1]
        ridgelist = pt['ridge']
        if len(ridgelist) < 2:
            return
        boundary = self.terrain.boundary
        inbound = self.terrain.inbound
        try:
            boundaryposition = {ir: boundary.index(self.ridgedict[ir]['polyline'][0]) for ir in ridgelist}
            inboundposition = {ir: inbound.index(self.ridgedict[ir]['polyline'][1]) for ir in ridgelist}
            outlist = sorted(ridgelist, key = lambda ir: boundaryposition[ir])
            inlist = sorted(ridgelist, key = lambda ir: inboundposition[ir])
            mismatched = [(outlist[i],inlist[i]) for i in range(len(outlist)) if outlist[i]!=inlist[i]]
            mismatched = {(min(i),max(i)) for i in mismatched}
            for position in mismatched:
                ir0 = position[0]
                ir1 = position[1]
                indice0in = [i for i, p in enumerate(inbound) if p == self.ridgedict[ir0]['polyline'][1]]
                indice1in = [i for i, p in enumerate(inbound) if p == self.ridgedict[ir1]['polyline'][1]]
                okindice0in = len(indice0in) == 1
                okindice1in = len(indice1in) == 1
                if not okindice0in or not okindice1in:
                    indice0out = [i for i, p in enumerate(boundary) if p == self.ridgedict[ir0]['polyline'][0]]
                    indice1out = [i for i, p in enumerate(boundary) if p == self.ridgedict[ir1]['polyline'][0]]
                    okindice0out = len(indice0out) == 1
                    okindice1out = len(indice1out) == 1
                    if okindice0out and okindice1out:
                        tmp = inboundposition[ir0]
                        inboundposition[ir0] = inboundposition[ir1]
                        inboundposition[ir1] = tmp
                    else:
                        if okindice0in and okindice1out:
                            i0 = outlist.index(ir0)
                            i1 = inlist.index(ir1)
                            if i0 == i1:
                                tmp = inboundposition[ir0]
                                inboundposition[ir0] = inboundposition[ir1]
                                inboundposition[ir1] = tmp
                            else:
                                print('x', position)
                        elif okindice1in and okindice0out:
                            i0 = inlist.index(ir0)
                            i1 = outlist.index(ir1)
                            if i0 == i1:
                                tmp = inboundposition[ir0]
                                inboundposition[ir0] = inboundposition[ir1]
                                inboundposition[ir1] = tmp
                            else:
                                print('+', position)
                        else:    
                            print(ir0, ir1, indice0in, indice1in, indice0out, indice1out)
            inlist = sorted(ridgelist, key = lambda ir: inboundposition[ir])
            mismatched = [(outlist[i],inlist[i]) for i in range(len(outlist)) if outlist[i]!=inlist[i]]
            mismatched = {(min(i),max(i)) for i in mismatched}
            ridgelist = inlist
        except:
            print('Exception orderRidges around virtual pit')
        self.nodedict[-1]['ridge'] = ridgelist
        print(self.nodedict[-1]['ridge'][21:25])
        
    def orderRidgesAroundNodes(self):
        """
        Orders the list of ridges around each node in a clockwise order 
        
        To be done: holes are not considered.
    
        Returns
        -------
        None.
    
        """
        for ipt, pt in self.nodedict.items():
            # -1 is the id of the virtual pit closing outside the domain
            ridgelist = pt['ridge']
            nr = len(ridgelist)
            if nr < 2:
                continue
            if ipt == -1:
                self.orderRidgesAroundVirtualPit()
            elif ipt < -1:
                # holes must be virtual pits with id<-1
                # get all the neighbours around the saddles
                continue
            else:    
                (i,j) = pt['ij']
                neighbours = [(i+di, j+dj) for (di,dj) in df.ldr8]
                # arrange the lines in the same order as the neighbours
                try:
                    ridgelist.sort(key = lambda it: neighbours.index(self.ridgedict[abs(it)]['polyline'][df.getSecondPointPosition(it)]))
                except (IndexError, ValueError) as err:
                    print('Exception raised in orderRidgesAroundNodes', err, ipt, pt, neighbours)
                # check if two ridges overlap at the node
                for i in range(nr):
                    ir1 = ridgelist[i-1]
                    ir2 = ridgelist[i]
                    absir1 = abs(ir1)
                    absir2 = abs(ir2)
                    hill1 = self.ridgedict[absir1]['hill']
                    hill2 = self.ridgedict[absir2]['hill']
                    secondpoint = self.ridgedict[absir1]['polyline'][df.getSecondPointPosition(ir1)]
                    if secondpoint == self.ridgedict[absir2]['polyline'][df.getSecondPointPosition(ir2)]:
                        thalweglist = pt['thalweg']
                        if hill1 == hill2:
                            for itt in thalweglist:
                                t = self.thalwegdict[abs(itt)]
                                if t['lefthill'] == hill1 and t['righthill'] == hill1 and secondpoint == t['polyline'][df.getSecondPointPosition(itt)]:
                                    it = itt
                                    break
                            try:
                                if ir1 > 0:
                                    ts1 = self.ridgedict[absir1]['startside']
                                else:
                                    ts1 = self.ridgedict[absir1]['endside']
                                if ir2 > 0:
                                    ts2 = self.ridgedict[absir2]['startside']
                                else:
                                    ts2 = self.ridgedict[absir2]['endside']
                            except (IndexError, ValueError, KeyError) as err:
                                print('Exception raised in orderRidgesAroundNodes', err)
                                print('if hill1 == hill2', ipt, it, ir1, ir2)
                            if ts1 == ts2:
                                print("Error at node", ipt, "ridges", ir1, "and", ir2, "on the same side")
                            # the one on the left should have a thalwegside of -1
                            # if not, we swap them
                            if (ts1 == 1 and it > 0) or (ts1 == -1 and it < 0):
                                ridgelist[i] = ir1
                                ridgelist[i-1] = ir2
                        else:
                            for it in thalweglist:
                                if it < 0:
                                    t = self.thalwegdict[-it]
                                    if t['polyline'][-2] != secondpoint:
                                        continue
                                    lefthill = t['righthill']
                                    righthill = t['lefthill']
                                else:
                                    t = self.thalwegdict[it]
                                    if t['polyline'][1] != secondpoint:
                                        continue
                                    lefthill = t['lefthill']
                                    righthill = t['righthill']
                                if hill2 == lefthill and hill1 == righthill:
                                    ridgelist[i-1] = ir2
                                    ridgelist[i] = ir1
    
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
        for ir in self.ridgedict:
            r = self.ridgedict[ir]
            r['leftdale'] = None
            r['rightdale'] = None
    
        self.daledict = {}
    
        # compute a dale around each thalweg
        for it in self.thalwegdict:
            t = self.thalwegdict[it]
            #if t['end'] == -1: # no dale for thalwegs ending outside
            #    t['dale'] = -1
            #    continue
            ipt = t['start'] # we start from the start node of the thalweg
            pt = self.nodedict[ipt]
    
            # we will start from the ridge on the left
            lefthill = t['lefthill']
            righthill = t['righthill']
            ridgelist = [ir for ir in pt['ridge'] if self.ridgedict[abs(ir)]['hill'] == lefthill]
            # if there is more than one ridge, we have to choose the closest one on the left
            ir = ridgelist[0] # take the one ridge in the list if there's only one
            lenridgelist = len(ridgelist)
            # if more than one ridge in the hill
            if lenridgelist > 1:
                (i,j) = pt['ij']
                # take the second point of each ridge
                secondpoints = [self.ridgedict[abs(jr)]['polyline'][df.getSecondPointPosition(jr)] for jr in ridgelist]
                thalwegpoint = t['polyline'][1]
                # if the thalweg overlaps a ridge
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
                            r = self.ridgedict[jr]
                            if 'startside' in r: 
                                if r['startside'] == 1: # on the right of the thalweg
                                    ir = ridgelist[k-1] # take the previous one
                                else: # on the left
                                    ir = jr # take this one
                        # if the ridge ends at node ipt
                        if jr < 0:
                            r = self.ridgedict[-jr]
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
                    if not found:
                        ir = ridgelist[-1]
            
            # we have the first ridge of the dale, we can turn around
            ridgering, ridgeline = self.getRidgeRing(ipt, ir)
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
                        r = self.ridgedict[rr]
                        if r['start'] == -1:
                            firstpoint = self.ridgedict[-lastrr]['polyline'][0]
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
                        pointlist += reversed(self.ridgedict[-rr]['polyline'][1:])
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

    def accumulateFlow(self):
        """
        computes the flow accumulation at the mouth of each thalweg. The flow
        accumulation is defined by the area of all dales upstream of the mouth.
        If a node has several downstream thalwegs, the accumulation can be sent
        to the thalwegs in proportion of the slope or to only one thalweg.

        Parameters
        ----------
        sd : boolean
            True if the network does not have diffluences (single direction).
            In that case, all the flow follows the steepest slope.
            Otherwise, the flow to each thalweg depends on the slope

        Returns
        -------
        None.

        """
        sd = self.singleflowdirection
        visited = {}
        for it in self.thalwegdict:
            t = self.thalwegdict[it]
            idale = t['dale']
            visited[it] = len(self.upstreamThalwegs(self.fromNode(it)))
            if idale != -1:
                t['accumulation'] = self.daledict[idale]['area']
            else:
                t['accumulation'] = 0
    
        streamlist = [self.downstreamThalwegs(ipt) for ipt in self.nodedict if self.isSpring(ipt)]
        thalweglist = [it for sublist in streamlist for it in sublist]
        for ist in thalweglist:
            stack = [ist]
            while stack:
                it = stack.pop()
                inode = self.toNode(it)
                downstream = self.downstreamThalwegs(inode)
                if downstream:
                    n = len(downstream)
                    slopelist = [min(0, self.computeStreamSlope(idt)) for idt in downstream]
                    if df.NOTIN in slopelist:
                        for i in range(n):
                            if slopelist[i] != df.NOTIN:
                                slopelist[i] = 0
                    accu = [0]*n
                    accutot = self.thalwegdict[it]['accumulation']
                    if sd:
                        maxslope = min(slopelist)
                        imax = slopelist.index(maxslope)
                        for i in range(n):
                            if i == imax:
                                accu[imax] = accutot
                            else:
                                accu[i] = 0
                                jt = downstream[i]
                                self.thalwegdict[jt]['headwater'] = True
                    else:
                        totalslope = sum(slopelist)
                        if totalslope == 0:
                            accu = [accutot/n]*n
                        else:
                            accu = [accutot*slopelist[i]/totalslope for i in range(n)]
                    for i in range(n):
                        idt = downstream[i]
                        self.thalwegdict[idt]['accumulation'] += accu[i]
                        visited[idt] -= 1
                        if visited[idt] == 0:
                            stack.append(idt)

    def computeDrainageBasins(self, thalweglist):
    	basins = []
    	for it in thalweglist:
            basin = self.computeDrainageBasin(it)
            basins.append(basin)
    	return basins

    def computeDrainageBasin(self, ithalweg):
        """
        Computes the drainage basin of a thalweg from the ridges. The method
        first finds all the ridges delineating the basin and then orders then
        to form a polygon.

        Parameters
        ----------
        ithalweg : integer
            The index of the thalweg whose drainage basin is computed.

        Returns
        -------
        basin : list of integers
            The list of directed ridges that form the drainage basin of ithalweg.

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
        ridgeset = plusridgeset | minusridgeset
        # ridgeset contains all the ridges on the basin boundary
        
        # take the ridge on the left of ithalweg
        inode = self.toNode(ithalweg) # the ridge is connected to inode
        node = self.nodedict[inode]
        ridges = node['ridge']
        #idale = t['dale']
        #boundary = self.daledict[idale]['boundary']
        # the first ridge is both in boundary and in ridges with the same sign
        #basin = [ir for ir in ridges if ir in boundary]
        basin = [ir for ir in ridges if ir in ridgeset]
        if len(basin) > 1:
            print("problem finding ridge left of thalweg", ithalweg)
        iridge = basin[0]
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
                    print("problem finding the next ridge", currentnode, iridge)
            else:
                iridge = ridges.pop()
            basin.append(iridge)
            if iridge > 0:
                ridge = self.ridgedict[iridge]
                currentnode = ridge['end']
            else:
                ridge = self.ridgedict[-iridge]
                currentnode = ridge['start']
        return basin

    def computeStrahlerOrderDownstream(self, idspring, accumulation):
        """
        Computation of the Strahler order starting from a spring and down to the outlet
        The Strahler order is assigned to both nodes and thalwegs.
        A thalweg has the same order as its from node
        For each node, we record the list of springs from which the streams flow
        We can then note if streams are coming from the same source, in which 
        case the order is not increased

        Parameters
        ----------
        idspring : TYPE
            DESCRIPTION.

        Returns
        -------
        None.

        """
        # start from spring idspring
        stack = []
        order = 0
        self.nodedict[idspring]['order'] = order
        
        # take all streams flowing from idspring
        downstream = self.downstreamThalwegs(idspring)
        for it in downstream:
            # a triplet formed by the tonode of the thalweg, the thalweg id and the from node order
            if self.thalwegdict[it]['accumulation'] < accumulation:
                order = 0
            else:
                order = 1
            stack.append((self.toNode(it), it, order))
        
        # depth first search of all the nodes below the spring
        while stack:
            # we have a node, an upper stream and the order of the from node
            (inode, iupthalweg, order) = stack.pop()
            node = self.nodedict[inode]
            # if the spring is already in node source, that means we already 
            # went through this node from the same spring in another way
            if idspring in node['source']:
                self.thalwegdict[iupthalweg]['order'] = order
                # the node takes the highest order of the streams above
                # if the order has been modified, we go down to update other 
                # streams and nodes below
                if node['order'] < order:
                    node['order'] = order
                    downstream = self.downstreamThalwegs(inode)
                    for it in downstream:
                        if 'headwater' in self.thalwegdict[it]:
                            torder = 0
                        else:
                            torder = order
                        if torder == 0 and self.thalwegdict[it]['accumulation'] >= accumulation:
                            stack.append((self.toNode(it), it, 1))
                        else:
                            stack.append((self.toNode(it), it, torder))
            else:
                # first time we go through this node from this spring
                node['source'].add(idspring)
                upstream = self.upstreamThalwegs(inode)
                upstream.remove(iupthalweg)
                self.thalwegdict[iupthalweg]['order'] = order
                # if the node order is higher, we stop there
                # if the order is equal, we increase the order and go down
                # if it's lower, assign the new order and keep going
                if upstream:
                    orders = [self.thalwegdict[k]['order'] for k in upstream]
                    oldorder = max(orders)
                    if order < oldorder:
                        continue # we stop going down, the stream reached a bigger stream
                    elif order == oldorder and order > 0:
                        order += 1
                node['order'] = order
                downstream = self.downstreamThalwegs(inode)
                for it in downstream:
                    if 'headwater' in self.thalwegdict[it]:
                        torder = 0
                    else:
                        torder = order
                    if torder == 0 and self.thalwegdict[it]['accumulation'] >= accumulation:
                        stack.append((self.toNode(it), it, 1))
                    else:
                        stack.append((self.toNode(it), it, torder))
                
    def assignStrahlerOrder(self, threshold = 0):
        """
        Compute the Strahler order for all streams in the network
        The order is computed starting from each spring
        All streams that have an accumulation lower than the threshold have an
        order of 0. Threshold is given in square units.

        Returns
        -------
        None.

        """
        # convert square units in square pixels
        accumulation = threshold/(self.terrain.x_size*self.terrain.x_size)
        # initialise all nodes with an order 0 and no spring
        for inode, node in self.nodedict.items():
            node['order'] = 0
            node['source'] = set()
    
        # set all thalweg orders to 0 by default 
        for it, t in self.thalwegdict.items():
            t['order'] = 0
        
        # take all nodes that are springs
        springlist = [key for key in self.nodedict if self.isSpring(key)]
        # we need to deal with the case where saddles are at the same height
        
        # go down the streams from the springs to assign the orders
        for inode in springlist:
            self.computeStrahlerOrderDownstream(inode, accumulation)

    def assignHortonOrder(self):
        # set all thalweg Horton orders to 0 by default 
        for it, t in self.thalwegdict.items():
            if 'order' not in t:
                print("Compute the Strahler order first")
                return
            t['horton'] = 0
        
        # get all sinks of the network
        outletlist = [ipt for ipt in self.nodedict if not self.downstreamThalwegs(ipt)]
        
        # set a stack of thalwegs to start from
        thalweglistlist = [self.upstreamThalwegs(ipt) for ipt in outletlist]
        thalwegstack = []
        for ilist in thalweglistlist:
            for it in ilist:
                order = self.thalwegdict[it]['order']
                if order > 0:
                    thalwegstack.append((it, order))
                
        while thalwegstack:
            (it, order) = thalwegstack.pop()
            if order <= self.thalwegdict[it]['horton']:
                continue
            self.thalwegdict[it]['horton'] = order
            upnode = self.fromNode(it)
            thalweglist = self.upstreamThalwegs(upnode)
            thalweglist = [i for i in thalweglist if self.thalwegdict[i]['order']>0]
            if not thalweglist:
                continue
            n = len(thalweglist)
            # get the thalwegs with the highest order
            orders = [self.thalwegdict[i]['order'] for i in thalweglist]
            maxorder = max(orders)
            sublist = [thalweglist[i] for i in range(n) if orders[i] == maxorder]
            # get the thalweg with the highest accumulation
            acculist = [self.thalwegdict[i]['accumulation'] for i in sublist]
            maxaccu = max(acculist)
            accuindex = acculist.index(maxaccu)
            mainthalweg = sublist[accuindex]
            # the thalweg is put in the stack with the horton order
            thalwegstack.append((mainthalweg, order))
            # other thalwegs are added with their own order
            for it in thalweglist:
                if it != mainthalweg and self.thalwegdict[it]['order'] > 0:
                    thalwegstack.append((it, self.thalwegdict[it]['order']))

    def addThalwegSlopeInDegree(self):
        for it, t in self.thalwegdict.items():
            if t['end'] == -1:
                t['slope'] = -90
                continue
            # we assume that x_size = y_size, otherwise, rewrite the length function
            length = df.length(t['polyline']) * self.terrain.x_size
            ptend = self.nodedict[t['end']]
            ptstart = self.nodedict[t['start']]
            dz = ptstart['z'] - ptend['z']
            slope = dz/length
            slope = 180*atan(slope)/pi
            t['slope'] = slope
    
    def addRidgeSlopeInDegree(self):
        for ir, r in self.ridgedict.items():
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

    def getHill(self, idpeak):
        """
        Returns the hill id that contains the peak idpeak

        Parameters
        ----------
        idpeak : integer
            id of the peak.

        Returns
        -------
        integer
            id of the hill defined in the hilldict.

        """
        ridges = self.nodedict[idpeak]['ridge']
        ridge = -ridges[0]
        return self.ridgedict[ridge]['hill']
    
    def getPeak(self, idhill):
        """
        Returns the peak contained in a hill

        Parameters
        ----------
        idhill : integer
            id of the hill.

        Returns
        -------
        ipt : integer
            id of the peak contained by the hill.

        """
        ridges = self.hilldict[idhill]['ridge']
        for ir in ridges:
            r = self.ridgedict[ir]
            ipt = r['end']
            if self.nodedict[ipt]['type'] == df.PEAK:
                return ipt
        print("Error no peak in hill")
        
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
        jnode : integer
            the id of the peak

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
        return path, jnode
    
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
        peaknode : integer
            the id of the peak

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
        path, peaknode = self.getPathToPeak(kmax)
        return path, peaknode

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
