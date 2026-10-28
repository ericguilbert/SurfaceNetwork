# -*- coding: utf-8 -*-
"""
Created on Tue Feb  4 15:19:05 2020

@author: ERGUI19
"""
from collections import deque
import array
import deffunction as df
from celltype import CellType
from surfacenetwork import SurfaceNetwork

class StreamNetwork(SurfaceNetwork):
    """
    StreamNetwork inherits from SurfaceNetwork and adds some functions 
    to compute flow direction and accumulation. This class is abstract and is
    there to define some functions later used by ExtendedSurfaceNetwork.
    
    In a thalweg, the starting node is the highest point (saddle or confluence)
    but the flow direction can be reversed to come out of puddles.
    """
    
    def getWaterCoursesAtNode(self, ipt):
        pt = self.nodedict[ipt]
        line = self.linedict
        return [i for i in pt['line'] if line[abs(i)]['type'] & CellType.WATERCOURSE_FLAG]

    def fromNode(self, it):
        """
        Return the upstream node of a stream. This depends on the flowdirection 
        attribute. If the flow direction is the same as the thalweg direction,
        from node is the start node otherwise it's the end node.

        Parameters
        ----------
        it : integer
            thalweg id in thalwegdict.

        Returns
        -------
        integer
            node id in nodedict.

        """
        if it < 0:
            it = -it
        t = self.linedict[it]
        if t['flowdirection']:
            return t['start']
        return t['end']
    
    def toNode(self, it):
        """
        Return the downstream node of a stream. This depends on the flowdirection 
        attribute. If the flow direction is the same as the thalweg direction,
        to node is the end node otherwise it's the start node.

        Parameters
        ----------
        it : integer
            thalweg id in thalwegdict.

        Returns
        -------
        integer
            node id in nodedict.

        """
        if it < 0:
            it = -it
        t = self.linedict[it]
        if t['flowdirection']:
            return t['end']
        return t['start']
    
    def isSink(self, ipit):
        """
        A node is a sink if it is a pit and is connected to no more than one
        thalweg, in opposition to other pits which are an outlet or a spurous
        pit through which the stream can flow. The only outlet is the virtual pit

        Parameters
        ----------
        ipit : integer
            node id in nodedict.

        Returns
        -------
        boolean
            true if the critical node is a sink.

        """
        node = self.nodedict[ipit]
        return node['type'] == CellType.PIT and len(self.getThalwegsAtNode(ipit)) == 1
    
    def isSpring(self, idsaddle):
        """
        A node is a spring if it is the from node of all the thalwegs it is connected to

        Parameters
        ----------
        idsaddle : integer
            node id in nodedict (not necessarily a saddle).

        Returns
        -------
        boolean
            true if the node is a spring.

        """
        tlist = self.getWaterCoursesAtNode(idsaddle)
        if not tlist:
            return False
        b = True
        for it in tlist:
            b = b and ((it > 0 and self.linedict[it]['flowdirection']) or (it < 0 and not self.linedict[-it]['flowdirection']))
        return b
    
    def computeStartSlope(self, it):
        """
        Compute the slope at the fromnode towards the next pixel.
        Downward slope is negative but can be positive if the stream is out
        of a puddle

        Parameters
        ----------
        it : integer
            index of the thalweg to compute.

        Returns
        -------
        TYPE: float
            Slope measured along the thalweg. If the thalweg is connected to
            the virtual pit, return a constant.

        """
        istart = self.fromNode(it)
        thalweg = self.linedict[it]
        if istart == thalweg['start']:
            ptend = thalweg['polyline'][1]
        else:
            ptend = thalweg['polyline'][-2]
        if not self.terrain.isInside(ptend[0], ptend[1]):
            return CellType.INNaN
        ptstart = self.nodedict[istart]['ij']
        length = df.length([ptstart, ptend])
        dz = self.terrain.dtm[ptend[0], ptend[1]] - self.nodedict[istart]['z']
        return dz / length

    def computeStreamSlope(self, it):
        """
        Compute the slope along a stream from the fromnode to the tonode.
        Divide the elevation difference by the thalweg length
        Downward slope is negative but can be positive if the stream is out
        of a puddle

        Parameters
        ----------
        it : integer
            index of the thalweg to compute.

        Returns
        -------
        TYPE: float
            Slope measured along the thalweg. If the thalweg is connected to
            the virtual pit, return a constant.

        """
        iend = self.toNode(it)
        if iend == -1:
            return CellType.INNaN
        length = df.length(self.linedict[it]['polyline'])
        ptend = self.nodedict[iend]
        ptstart = self.nodedict[self.fromNode(it)]
        dz = ptend['z'] - ptstart['z']
        return dz / length

    def downstreamThalwegs(self, idnode):
        """
        returns the list of streams that flow from this node, that is all the 
        thalwegs for which the node is a from node

        Parameters
        ----------
        idnode : integer
            node id in nodedict.

        Returns
        -------
        downstream : list
            the list of thalweg id in thalwegdict.

        """
        tlist = self.getWaterCoursesAtNode(idnode)
        downstream = []
        for it in tlist:
            if self.fromNode(it) == idnode:
                downstream.append(abs(it))
        return downstream
    
    def upstreamThalwegs(self, idnode):
        """
        returns the list of streams that flow to this node, that is all the 
        thalwegs for which the node is a to node

        Parameters
        ----------
        idnode : integer
            node id in the nodedict.

        Returns
        -------
        upstream : list
            the list of thalweg id in thalwegdict.

        """
        tlist = self.getWaterCoursesAtNode(idnode)
        upstream = []
        for it in tlist:
            if self.toNode(it) == idnode:
                upstream.append(abs(it))
        return upstream
        
    def preserveSinks(self):
        """
        Set the flow direction towards the sink for all thalwegs connecting to
        a sink to avoid having sinks appearing as springs

        Returns
        -------
        None.

        """
        sinklist = [key for key in self.nodedict if self.isSink(key)]
        for i in sinklist:
            it = self.getThalwegsAtNode(i)[0]
            self.linedict[-it]['flowdirection'] = True
                   
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
            if self.linedict[it]['accumulation'] < accumulation:
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
                self.linedict[iupthalweg]['order'] = order
                # the node takes the highest order of the streams above
                # if the order has been modified, we go down to update other 
                # streams and nodes below
                if node['order'] < order:
                    node['order'] = order
                    downstream = self.downstreamThalwegs(inode)
                    for it in downstream:
                        if 'headwater' in self.linedict[it]:
                            torder = 0
                        else:
                            torder = order
                        if torder == 0 and self.linedict[it]['accumulation'] >= accumulation:
                            stack.append((self.toNode(it), it, 1))
                        else:
                            stack.append((self.toNode(it), it, torder))
            else:
                # first time we go through this node from this spring
                node['source'].add(idspring)
                upstream = self.upstreamThalwegs(inode)
                upstream.remove(iupthalweg)
                self.linedict[iupthalweg]['order'] = order
                # if the node order is higher, we stop there
                # if the order is equal, we increase the order and go down
                # if it's lower, assign the new order and keep going
                if upstream:
                    orders = [self.linedict[k]['order'] for k in upstream]
                    oldorder = max(orders)
                    if order < oldorder:
                        continue # we stop going down, the stream reached a bigger stream
                    elif order == oldorder and order > 0:
                        order += 1
                node['order'] = order
                downstream = self.downstreamThalwegs(inode)
                for it in downstream:
                    if 'headwater' in self.linedict[it]:
                        torder = 0
                    else:
                        torder = order
                    if torder == 0 and self.linedict[it]['accumulation'] >= accumulation:
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
        # convert square units in pixels
        accumulation = threshold/(self.terrain.x_size*self.terrain.x_size)
        # initialise all nodes with an order 0 and no spring
        for inode, node in self.nodedict.items():
            node['order'] = 0
            node['source'] = set()
    
        # set all thalweg orders to 0 by default 
        for it, t in self.linedict.items():
            if t['type'] & CellType.WATERCOURSE_FLAG:
                t['order'] = 0
        
        # take all nodes that are springs
        springlist = [key for key in self.nodedict if self.isSpring(key)]
        # we need to deal with the case where saddles are at the same height
        
        # go down the streams from the springs to assign the orders
        for inode in springlist:
            self.computeStrahlerOrderDownstream(inode, accumulation)

        # remove order and spring node
        for inode, node in self.nodedict.items():
            del node['order']
            del node['source']

    def assignHortonOrder(self):
        """
        Computation of the Horton order for each thalweg. The Strahler order
        must have been computed first. If two upstream thalwegs have the same
        Strahler order, the main stream is the one with highest accumulation.
        Thalwegs with a Strahler order of 0, meaning no water, have a Horton 
        order of 0 also.

        Returns
        -------
        None.

        """
        # set all thalweg Horton orders to 0 by default 
        for it, t in self.linedict.items():
            if not (t['type'] & CellType.WATERCOURSE_FLAG):
                continue
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
                order = self.linedict[it]['order']
                if order > 0:
                    thalwegstack.append((it, order))
                
        while thalwegstack:
            (it, order) = thalwegstack.pop()
            if order <= self.linedict[it]['horton']:
                continue
            self.linedict[it]['horton'] = order
            upnode = self.fromNode(it)
            thalweglist = self.upstreamThalwegs(upnode)
            thalweglist = [i for i in thalweglist if self.linedict[i]['order']>0]
            if not thalweglist:
                continue
            n = len(thalweglist)
            # get the thalwegs with the highest order
            orders = [self.linedict[i]['order'] for i in thalweglist]
            maxorder = max(orders)
            sublist = [thalweglist[i] for i in range(n) if orders[i] == maxorder]
            # get the thalweg with the highest accumulation
            acculist = [self.linedict[i]['accumulation'] for i in sublist]
            maxaccu = max(acculist)
            accuindex = acculist.index(maxaccu)
            mainthalweg = sublist[accuindex]
            # the thalweg is put in the stack with the horton order
            thalwegstack.append((mainthalweg, order))
            # other thalwegs are added with their own order
            for it in thalweglist:
                if it != mainthalweg and self.linedict[it]['order'] > 0:
                    thalwegstack.append((it, self.linedict[it]['order']))
                
    def assignFlowDirection(self, puddledict):
        """
        Compute the flow direction of each thalweg with consideration of puddles.
        Puddles are previously computed, normally from the Depression Graph.

        Parameters
        ----------
        puddledict : dictionary
            The dictionary contains puddles in which the flow is directed towards
            the outlet. Keys are the puddle indexes. If they are computed from 
            the Depression Graph, they are the depression ids and are the same
            as the lowest pit of the depression. Values are pairs containing
            the id of the outlet and the set of thalwegs of the depression
            below the outlet.

        Returns
        -------
        None.

        """
        # by default, the flow direction is oriented along the thalweg and the
        # stream order is 0 (that means no order is assigned)
        for it, t in self.linedict.items():
            if t['type'] & CellType.WATERCOURSE_FLAG:
                t['flowdirection'] = True
                t['order'] = 0
    
        # We go through all puddles to assign the correct flow direction
        for ip, p in puddledict.items():
            # get the outlet of the depression
            inode = p[0]
            # get all the thalwegs whose highest node is strictly lower than the outlet
            thalwegs = p[1]
            # We take the thalwegs connecting to s and set their flow direction towards s
            nodestack = deque([inode])
            while nodestack:
                inode = nodestack.pop()
                s = self.nodedict[inode]
                tlist = [it for it in s['line'] if abs(it) in thalwegs]
                for it in tlist:
                    thalwegs.remove(abs(it))
                    if it > 0:
                        t = self.linedict[it]
                        tend = t['end']
                        t['flowdirection'] = False
                        nodestack.appendleft(tend)
                    if it < 0:
                        t = self.linedict[-it]
                        tstart = t['start']
                        t['flowdirection'] = True
                        nodestack.appendleft(tstart)
   
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
        for it in self.linedict:
            t = self.linedict[it]
            if not (t['type'] & CellType.WATERCOURSE_FLAG):
                continue
#            if t['type'] & CellType.UNDERGROUND_FLAG:
#                t['accumulation'] = 0
#                continue
            visited[it] = len(self.upstreamThalwegs(self.fromNode(it)))
            if 'headwater' in t:
                del t['headwater']
            if 'dale' in t and t['dale'] != -1:
                idale = t['dale']
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
                    if CellType.INNaN in slopelist:
                        for i in range(n):
                            if slopelist[i] != CellType.INNaN:
                                slopelist[i] = 0
                    accu = [0]*n
                    accutot = self.linedict[it]['accumulation']
                    if sd:
                        maxslope = min(slopelist)
                        imax = slopelist.index(maxslope)
                        for i in range(n):
                            if i == imax:
                                accu[imax] = accutot
                            else:
                                accu[i] = 0
                                jt = downstream[i]
                                self.linedict[jt]['headwater'] = True
                    else:
                        totalslope = sum(slopelist)
                        if totalslope == 0:
                            accu = [accutot/n]*n
                        else:
                            accu = [accutot*slopelist[i]/totalslope for i in range(n)]
                    for i in range(n):
                        idt = downstream[i]
                        self.linedict[idt]['accumulation'] += accu[i]
                        visited[idt] -= 1
                        if visited[idt] == 0:
                            stack.append(idt)

    def computeDrainageBasins(self, thalweglist):
    	basins = []
    	for it in thalweglist:
            basin = self.computeDrainageBasin(it)
            basins.append(basin)
    	return basins
### On est ici
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
        watercourse = [t for t in self.linedict if self.linedict[t] & CellType.WATERCOURSE_FLAG]
        t = self.linedict[ithalweg]
        if t['dale'] < 1:
            return None
        # collect all the dales above the thalweg
        visited = {} # mark if a thalweg has been visited or not
        # this is necessary if distributaries are allowed
        for it in watercourse:
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
            if 'dale' in self.linedict[it]:
                dalelist.append(self.linedict[it]['dale'])
            if 'headwater' in self.linedict[it]:
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
                    print("problem finding the next ridge", currentnode, iridge, ithalweg)
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

    def getMaxOrder(self, order = 'order'):
        return max([self.thalwegdict[i][order] for i in self.thalwegdict])
    
    def aggregateStreams(self):
        """
        Remove nodes and aggregate continuous streams
        At the end of the process, the stream only contains confluences,
        saddles corresponding to springs and pits corresponding to sinks.
        All other nodes should have been removed.
        Such a network is no longer a thalweg network as part of a surface
        network but a river network as used in cartography
        
        The method cannot be called for a surface network because connection
        with ridges is lost and hill and dale definition becomes inconsistent
    
        Returns
        -------
        None.
    
        """
        
        # select all the nodes that are connected to exactly two streams
        # directed the same way
        flownodestack = [key for key in self.nodedict if len(self.nodedict[key]['thalweg']) == 2]
        while flownodestack:
            fn = flownodestack.pop()
            lt = self.nodedict[fn]['thalweg']
            it1 = lt[0]
            it2 = lt[1]
            flow1 = self.thalwegdict[abs(it1)]['flowdirection']
            flow2 = self.thalwegdict[abs(it2)]['flowdirection']
            flowfrom1 = ( (it1 > 0) == flow1 ) # true if the flow starts at the node
            flowfrom2 = ( (it2 > 0) == flow2 )
            if flowfrom1 != flowfrom2: # one flow arrives, the other leaves at the node
                # the node is removed
                startnode = 0
                endnode = 0
                poly1 = []
                poly2 = []
                order = self.nodedict[fn]['order']
                if flowfrom1: # flow from node to thalweg 1
                    # flow goes from thalweg 2 to thalweg 1    
                    startnode = self.fromNode(it2)
                    endnode = self.toNode(it1)
                    firstit = abs(it2)
                    lastit = abs(it1)
                    if flow2:
                        poly1 = self.thalwegdict[-it2]['polyline'][:]
                    else:
                        poly1 = self.thalwegdict[it2]['polyline'][::-1]
                    if flow1:
                        poly2 = self.thalwegdict[it1]['polyline'][:]
                    else:
                        poly2 = self.thalwegdict[-it1]['polyline'][::-1]
                else: # flow from thaweg 1 to node
                    # flow goes from thalweg 1 to thalweg 2
                    startnode = self.fromNode(it1)
                    endnode = self.toNode(it2)
                    firstit = abs(it1)
                    lastit = abs(it2)
                    if flow2:
                        poly2 = self.thalwegdict[it2]['polyline'][:]
                    else:
                        poly2 = self.thalwegdict[-it2]['polyline'][::-1]
                    if flow1:
                        poly1 = self.thalwegdict[-it1]['polyline'][:]
                    else:
                        poly1 = self.thalwegdict[it1]['polyline'][::-1]
                # now we can build the new thalweg
                # it replaces the first thalweg
                self.thalwegdict[firstit] = { 'start': startnode,
                                              'polyline': poly1 + poly2,
                                              'end': endnode,
                                              'flowdirection': True,
                                              'order': order
                    }
                # the second thalweg is deleted
                del self.thalwegdict[lastit]
                # update the start node
                t = self.nodedict[startnode]['thalweg']
                if -firstit in t:
                    i = t.index(-firstit)
                    t[i] = firstit
                # update the end node
                t = self.nodedict[endnode]['thalweg']
                if -lastit in t:
                    i = t.index(-lastit)
                else:
                    i = t.index(lastit)
                t[i] = -firstit
                # delete the middle node
                del self.nodedict[fn]

    def pruneStreams(self, lengthmax):
        """
        Prune the streams by removing all the streams of order 1 and shorter 
        than lengthmax. Streams should have been previously aggregated.
    
        Parameters
        ----------
        ordermax : integer
            maximum order of streams to be removed.
        lengthmax : float
            maximum length to be removed.
    
        Returns
        -------
        None.
    
        """
        
        # get all the springs
        springlist = [key for key in self.nodedict if self.isSpring(key)]
        for key in springlist:
            tlist = self.nodedict[key]['thalweg'][:]
            for it in tlist:
                ait = abs(it)
                # calculate thalweg length
                length = df.length(self.thalwegdict[ait]['polyline'])
                if length <= lengthmax:
                   # the thalweg is removed from the spring's thalweg list
                   self.nodedict[key]['thalweg'].remove(it)
                   # the thalweg is removed from the other node's list
                   iend = self.toNode(it)
                   self.nodedict[iend]['thalweg'].remove(-it)
                   # the thalweg is removed from the thalweg list
                   del self.thalwegdict[ait]
                   # if the end node is no longer connected to a thalweg
                   if not self.nodedict[iend]['thalweg']:
                       del self.nodedict[iend]
            # if all thalwegs were removed, remove the node
            if not self.nodedict[key]['thalweg']:
                del self.nodedict[key]
    
    def pruneNetwork(self, order, length):
        o = order
        while o > 0:
            self.aggregateStreams()
            self.pruneStreams(length)
            self.assignStrahlerOrder()
            o -= 1
        self.aggregateStreams()

    def addStreamLength(self):
        """
        Add to each stream an attribute with its length

        Returns
        -------
        None.

        """
        thalwegdict = self.thalwegdict
        cellsize = self.terrain.x_size
    
        for id_key, thalweg in thalwegdict.items():
            polyline = thalweg['polyline']
            tlength = df.length(polyline)
            thalweg['length'] = cellsize * tlength
    