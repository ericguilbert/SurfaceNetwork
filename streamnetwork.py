# -*- coding: utf-8 -*-
"""
Created on Tue Feb  4 15:19:05 2020

@author: ERGUI19
"""
from collections import deque
import array
import deffunction as df
from thalwegnetwork import ThalwegNetwork

class StreamNetwork(ThalwegNetwork):
    """
    StreamNetwork stores a network of thalwegs. It inherits from
    ThalwegNetwork and adds some functions to compute a flow direction
    
    In a thalweg, the starting node is the highest point (saddle or confluence)
    In a stream, the starting and end nodes depend on the flow direction
    This flow direction may be reversed to handle spurious pits
    """
    
    def buildStreamsFromTerrain(self, terrain):
        self.buildThalwegsFromTerrain(terrain)
        print("Create puddles from network")
        self.addPuddles()
        print("Merge puddles and assign flow direction to streams")
        self.mergePuddlesAndAssignFlow()
        print("Compute the Strahler order of each stream")
        self.assignStrahlerOrder()

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
        t = self.thalwegdict[it]
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
        t = self.thalwegdict[it]
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
        return node['type'] == df.PIT and len(node['thalweg']) == 1
    
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
        s = self.nodedict[idsaddle]
        if not s['thalweg']:
            return False
        b = True
        for it in s['thalweg']:
            b = b and ((it > 0 and self.thalwegdict[it]['flowdirection']) or (it < 0 and not self.thalwegdict[-it]['flowdirection']))
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
        thalweg = self.thalwegdict[it]
        if istart == thalweg['start']:
            ptend = thalweg['polyline'][1]
        else:
            ptend = thalweg['polyline'][-2]
        if not self.terrain.isInside(ptend[0], ptend[1]):
            return df.NOTIN
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
            return df.NOTIN
        length = df.length(self.thalwegdict[it]['polyline'])
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
        node = self.nodedict[idnode]
        tlist = node['thalweg']
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
        node = self.nodedict[idnode]
        tlist = node['thalweg']
        upstream = []
        for it in tlist:
            if self.toNode(it) == idnode:
                upstream.append(abs(it))
        return upstream

    def floodPuddle(self, idseed):
        """
        Floods an area starting from a point until an outlet is found. The area
        is a puddle composed of all points around idseed.
        Only one outlet is possible, if there are two candidate nodes, we keep the
        one with the lowest neighbour. This can be discussed since we could also 
        choose the steepest slope or any other rule.
        One question: what to do if passes have neighbours at the same height?
    
        Parameters
        ----------
        idpt : list of nodes forming a puddle.
    
        Returns
        -------
        puddle defined by a list of nodes.
        outlet is the node (a pass) where the flow pours out.
    
        """
        puddle = list(idseed) # list of points forming the puddle
        puddlestack = list(idseed) # stack of points to visit
        outlet = puddle[0]
        #outlet = []
        outflow = outlet
        zoutlet = self.terrain.maxz
        zoutflow = zoutlet
        while (puddlestack):
            ipt = puddlestack.pop()
            z = self.nodedict[ipt]['z']
            
            #print("début", ipt, puddlestack)
            # séparer les voisins au dessus et en dessous
            lowerneighbours, higherneighbours = self.getThalwegNeighbourNodes(ipt)
            #print("voisins", lowerneighbours, higherneighbours)
            for iv in lowerneighbours: # water can flow outside the puddle from there
                if not iv in puddle: # if the neighbour is already in the puddle, nothing to do
                    if ipt == outlet: # we are dealing with the outlet
                        if self.computeSlope(iv, ipt) > self.computeSlope(outflow, ipt):
                        #if self.nodedict[iv]['z'] < zoutflow: 
                            outflow = iv # and it flows into iv
                            zoutflow = self.nodedict[outflow]['z']                        
                    elif z == zoutlet and outflow > -1: # problem can occur if holes are not nan
                        if self.computeSlope(iv, ipt) > self.computeSlope(outflow, outlet):
                        #if self.nodedict[iv]['z'] < zoutflow: 
                            outlet = ipt # ipt is added as an outlet
                            outflow = iv # and it flows into iv
                            zoutflow = self.nodedict[outflow]['z']
                        elif iv != outflow and iv > -1 and self.nodedict[iv]['z'] == zoutflow:
                            #print("Cas spécial zoutlet deux solutions", iv, outflow)
                            po = self.nodedict[outflow]['ij']
                            pi = self.nodedict[iv]['ij']
                            if df.lexicoSign(pi[0] - po[0], pi[1] - po[1]) < 0:
                                outlet = ipt # ipt is added as an outlet
                                outflow = iv # and it flows into iv
                                zoutflow = self.nodedict[outflow]['z']
                    elif z < zoutlet: # ipt is a new outlet
                        zoutlet = z # so we have a new outlet height
                        outlet = ipt
                        outflow = iv
                        zoutflow = self.nodedict[outflow]['z']
    
            for iv in higherneighbours:
                if not iv in puddle: # water cannot flow, we extend the puddle
                    # à vérifier: peut-on simplifier en ajoutant and v < zoutlet?
                    v = self.nodedict[iv]['z']
                    if v < zoutlet and iv != outflow:
                        puddle.append(iv)
                        puddlestack.append(iv)
            #print("fin", outlet, outflow, puddle)
    
        # we remove all nodes that are higher than the outlet
        for ipt in puddle:
            puddle={i for i in puddle if self.nodedict[i]['z']<=zoutlet}
        puddle = puddle.union(idseed)
        
        return puddle, outlet, outflow
        
    def addPuddles(self):
        """
        Add puddles to the drainage. A puddle is defined by a pit and the area 
        up to its lowest saddle.
        To be checked: if a pit is connected to only one thalweg, should we 
        define a puddle since it does not interrupt any flow?
    
        Returns
        -------
        None.
    
        """
        self.puddledict = {}
        #print("Warning: Pits with only one thalwegs are not considered as puddles")
        for ipit, pit in self.nodedict.items():
#            if pit['type'] == df.PIT and ipit != -1 and len(pit['thalweg']) > 1:
            if pit['type'] == df.PIT and ipit != -1:
                puddle, outlet, outflow = self.floodPuddle([ipit])
                self.puddledict[ipit] = {'nodes': puddle,
                                         'outlet': outlet,
                                         'outflow': outflow}

    def mergePuddlesSameOutlet(self):
        """
        Merge puddles sharing the same outlet into larger puddles.
        At the end of the process, an outlet can be located in another puddle.
        That means that both puddles overlap. It should not be a problem unless
        in some test, puddles must be disjoint.
    
        Returns
        -------
        None.
    
        """
        # sort puddles for each outlet
        outletidx = {}
        for ipuddle, puddle in self.puddledict.items():
            outlet = puddle['outlet']    
            if outlet in outletidx:
                outletidx[outlet].append(ipuddle)
            else:
                outletidx[outlet] = [ipuddle]
    
        # list of all outlets of puddles to be merged
        listtomerge = []
        for idx in outletidx:
            lpuddle = outletidx[idx]
            if len(lpuddle) > 1:
                listtomerge.append(idx)
            
        while listtomerge:
            for idx in listtomerge:
                lpuddle = outletidx[idx]
                seed = []
                # put all the nodes of the puddles in the same set
                for jpuddle in lpuddle:
                    seed.extend(self.puddledict[jpuddle]['nodes'])
                nodes, outlet, outflow = self.floodPuddle(seed)
                # the new puddle is stored in place of the first old puddle
                ipuddle = lpuddle.pop()
                self.puddledict[ipuddle] = {'nodes': nodes,
                                            'outlet': outlet,
                                            'outflow': outflow}
                # other puddles are deleted
                for jpuddle in lpuddle:
                    del self.puddledict[jpuddle]
                del outletidx[idx]
                # if the new outlet is already the outlet of another puddle
                # we add it to the puddle list
                if outlet in outletidx:
                    outletidx[outlet].append(ipuddle)
                else:
                    outletidx[outlet] = [ipuddle]
            # rebuild the list of puddles to merge
            listtomerge = []
            for idx in outletidx:
                lpuddle = outletidx[idx]
                if len(lpuddle) > 1:
                    listtomerge.append(idx)
    
    def mergePuddlesSameOutflow(self):
        """
        merge together puddles that have the same outflow and no other stream
        to flow to, meaning that the flow is stuck there. The outflow is the
        node located downstream the outlet.

        Returns
        -------
        None.

        """
        # sort puddles per outflow
        outflowidx = {}
        for ipuddle, puddle in self.puddledict.items():
            outflow = puddle['outflow']
            if outflow in outflowidx:
                outflowidx[outflow].append(ipuddle)
            else:
                outflowidx[outflow] = [ipuddle]
        
        outflowstack = []
        listtomerge = []
        for idx in outflowidx:
            lpuddle = outflowidx[idx]
            # if several puddles have the same outflow
            if len(lpuddle) > 1:
                outflowstack.append(idx)

        # outflowstack contains the list of outflows to be processed
        while outflowstack:
            outflow = outflowstack.pop()
            # get the list of nodes connected to the outflow by a thalweg
            lowerneighbours, higherneighbours = self.getThalwegNeighbourNodes(outflow)
            nodeset = set()
            # get the list of all nodes in the puddles connected to outflow
            for ipuddle in outflowidx[outflow]:
                nodeset = nodeset.union(self.puddledict[ipuddle]['nodes'])
            # remainder are the neighbouring nodes of the outflow but are not in the puddles
            remainder = lowerneighbours.union(higherneighbours).difference(nodeset)
            # if remainder is empty, no exit, the puddles have to be merged
            if not remainder:
                listtomerge.append(outflow)
        # create the new puddles
        for outf in listtomerge:
            # create a set of nodes containing the outflow and all puddle nodes
            seed = {outf}
            lpuddle = outflowidx[outf]
            print(lpuddle)
            for jpuddle in lpuddle:
                seed = seed.union(self.puddledict[jpuddle]['nodes'])
            # build the new puddle
            nodes, outlet, outflow = self.floodPuddle(seed)
            ipuddle = lpuddle.pop()
            # store the new puddle
            self.puddledict[ipuddle] = {'nodes': nodes,
                                        'outlet': outlet,
                                        'outflow': outflow}
            # delete the old puddles
            for jpuddle in lpuddle:
                del self.puddledict[jpuddle]
    
    def mergePuddleLoops(self):
        """
        Merge puddles that flow into each other, creating a loop in the flow

        Returns
        -------
        None.

        """
        outflowidx = {}
        # sort puddles by outflow
        for ipuddle, puddle in self.puddledict.items():
            outflow = puddle['outflow']
            if outflow in outflowidx:
                outflowidx[outflow].append(ipuddle)
            else:
                outflowidx[outflow] = [ipuddle]
    
        pairstobemerged = set()
        for ipuddle in self.puddledict:
            # puddlefilter contains the outflows located in ipuddle
            puddlefilter = filter(lambda i: i in outflowidx, self.puddledict[ipuddle]['nodes'])
            for puddlef in puddlefilter:
                # the outflow of ipuddle
                ioutflow = self.puddledict[ipuddle]['outflow']
                # list of puddles that flow in ipuddle
                puddlelist = outflowidx[puddlef]
                # look if ipuddle flows inside one of these puddles
                for jpuddle in puddlelist:
                    if ioutflow in self.puddledict[jpuddle]['nodes']:
                        # it means that ipuddle flows in jpuddle and jpuddle flows in ipuddle
                        tmp = (ipuddle, jpuddle)
                        pairstobemerged.add(tmp)
        # merge the puddles that flow into each other
        #print(pairstobemerged)
        tobemerged = []
        while pairstobemerged:
            loop = set(pairstobemerged.pop())
            tobedeleted = []
            for pair in pairstobemerged:
                if not loop.isdisjoint(pair):
                    loop.update(pair)
                    tobedeleted.append(pair)
            for pair in tobedeleted:    
                pairstobemerged.discard(pair)
            tobemerged.append(loop)
            
        print(tobemerged)
        for lpuddle in tobemerged:
            ipuddle = lpuddle.pop()
            seed = self.puddledict[ipuddle]['nodes']
            for jpuddle in lpuddle:
                # put all the nodes of both puddles in one set
                seed.update(self.puddledict[jpuddle]['nodes'])
            nodes, outlet, outflow = self.floodPuddle(seed)
            self.puddledict[ipuddle] = {'nodes': nodes,
                                        'outlet': outlet,
                                        'outflow': outflow}
            for jpuddle in lpuddle:
                del self.puddledict[jpuddle]

    def mergeIntersectingPuddles(self):
        #keys = [key for key in self.puddledict]
        keylist = [key for key in self.puddledict]
        keys = array.array('l', keylist)
        
        tobemerged = []
        n = len(keys)
        for i in range(n-1):
            ipuddle = keys[i]
            nodes = self.puddledict[ipuddle]['nodes']
            for j in range(i+1, n):
                jpuddle = keys[j]
                intersection = nodes.intersection(self.puddledict[jpuddle]['nodes'])
                ninter = len(intersection)
                if ninter > 1:
                        tobemerged.append((ipuddle, jpuddle))
        replace = {}
        print(tobemerged)
        for (i,j) in tobemerged:
            try:
                if i in replace:
                    i = replace[i]
                while j in replace:
                    j = replace[j]
                if i == j:
                    continue
                seed = self.puddledict[i]['nodes'].union(self.puddledict[j]['nodes'])
            except KeyError as err:
                print("mergeIntersectingPuddles KeyError", err, "keys", i, j)
            # build the new puddle
            nodes, outlet, outflow = self.floodPuddle(seed)
            # store the new puddle
            self.puddledict[i] = {'nodes': nodes,
                                  'outlet': outlet,
                                  'outflow': outflow}
            # delete the old puddles
            del self.puddledict[j]
            replace[j] = i
            print(replace)
        
    # def mergeIntersectingPuddles(self):
    #     outletidx = {self.puddledict[ip]['outlet']: ip for ip in self.puddledict}
    #     # sort the outlets based on their elevation, starting from the lowest
    #     outletlist = sorted(outletidx.keys(), key = lambda i: self.nodedict[i]['z'])
        
    #     tobedeleted = []
    #     n = len(outletlist)
    #     for i in range(n-1):
    #         iout = outletlist[i]
    #         ipuddle = outletidx[iout]
    #         for j in range(i+1, n):
    #             jout = outletlist[j]
    #             jpuddle = outletidx[jout]
    #             if iout in self.puddledict[jpuddle]['nodes']:
    #                 if self.puddledict[ipuddle]['nodes'].issubset(self.puddledict[jpuddle]['nodes']):
    #                     tobedeleted.append(ipuddle)
    #     for ip in tobedeleted:
    #         del self.puddledict[ip]

    def mergePuddles(self):
        """
        merges puddles until no more merging can be done.
        In a loop, merge puddles with the same outlet then with the same outflow
        without an exit. In outer loop, merge those looping puddles and repeat

        Returns
        -------
        None.

        """
        number = len(self.puddledict)
        number1 = number + 1
        while (number<number1):
            while (number < number1):
                self.mergePuddlesSameOutlet()
                self.mergePuddlesSameOutflow()
                number1 = number
                number = len(self.puddledict)
            self.mergePuddleLoops()
            self.mergeIntersectingPuddles()
            number1 = number
            number = len(self.puddledict)
        
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
            it = self.nodedict[i]['thalweg'][0]
            self.thalwegdict[-it]['flowdirection'] = True
        # wrongspringlist = [key for key in self.nodedict if self.isSpring(key) and self.nodedict[key]['type'] != df.SADDLE]
        # for i in wrongspringlist:
        #     print("wronglist:", i)
        #     thalweg = [it for it in self.nodedict[i]['thalweg'] if it < 0]
        #     endpoint = [self.toNode(i) for i in thalweg]
        #     imax = max(endpoint, key = lambda i: self.nodedict[i]['z'])
        #     lit = [it for it in thalweg if self.thalwegdict[-it]['start'] == imax]
        #     it = -lit[0]
        #     #lit = [abs(it) for it in thalweg if self.thalwegdict[abs(it)]['start'] == imax]
        #     #it = lit[0]
        #     self.thalwegdict[it]['flowdirection'] = True
                
    def assignFlowDirection(self):
        # by default, the flow direction is oriented along the thalweg and the
        # stream order is 0 (that means no order is assigned)
        for it, t in self.thalwegdict.items():
            t['flowdirection'] = True
            t['order'] = 0
    
        # We go through all puddles to assign the correct flow direction
        for ip, p in self.puddledict.items():
            inode = p['outlet']
            nodes = list(p['nodes'])
            nodestack = deque([inode])
            # We take the thalwegs connecting to s and set their flow direction towards s
            while nodestack:
                inode = nodestack.pop()
                nodes.remove(inode)
                s = self.nodedict[inode]
                tlist = s['thalweg']
                for it in tlist:
                    if it > 0:
                        t = self.thalwegdict[it]
                        tend = t['end']
                        if tend in nodes:
                            t['flowdirection'] = False
                            if tend not in nodestack:
                                nodestack.appendleft(tend)
                    if it < 0:
                        t = self.thalwegdict[-it]
                        tstart = t['start']
                        if tstart in nodes:
                            t['flowdirection'] = True
                            if tstart not in nodestack:
                                nodestack.appendleft(tstart)
        #self.preserveSinks()
   
    def flowDownStream(self, idspring):
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
        order = 1
        self.nodedict[idspring]['order'] = order
        
        # take all streams flowing from idspring
        downstream = self.downstreamThalwegs(idspring)
        for it in downstream:
            # a triplet formed by the tonode of the thalweg, the thalweg id and the from node order
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
                        stack.append((self.toNode(it), it, order))
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
                    elif order == oldorder:
                        order += 1
                node['order'] = order
                downstream = self.downstreamThalwegs(inode)
                for it in downstream:
                    stack.append((self.toNode(it), it, order))
                
    def assignStrahlerOrder(self):
        """
        Compute the Strahler order for all streams in the network
        The order is computed starting from each spring

        Returns
        -------
        None.

        """
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
        
        print(len(springlist))
        # go down the streams from the springs to assign the orders
        for inode in springlist:
            self.flowDownStream(inode)

    # def assignStrahlerOrder(self):
    #     """
    #     Compute the Strahler order for all streams in the network
    #     The order is computed starting from each spring

    #     Returns
    #     -------
    #     None.

    #     """
    #     nodeorder = {}
    #     spring = {}
    #     # initialise all nodes with an order 0 and no spring
    #     for inode in self.nodedict.items():
    #         spring[inode] = set()
    #         nodeorder[inode] = 0
    
    #     # set all thalweg orders to 0 by default 
    #     for it, t in self.thalwegdict.items():
    #         t['order'] = 0
        
    #     # take all nodes that are springs
    #     nodestack = [key for key in self.nodedict if self.isSpring(key)]
    #     # we need to deal with the case where saddles are at the same height
        
    #     while nodestack:
    #         inode = nodestack.pop()
    #         uplist = self.upstreamThalwegs(inode)
    #         if uplist:
    #             for it in uplist:
    #                 t = self.thalwegdict[it]
                    
    #         else:
    #             nodeorder[inode] = 1
    #         tlist = self.downstreamThalwegs(inode)
    #         for t in tlist:
    #             t.order = nodeorder[inode]

    def detectLoops2(self):
        """
        Detects loops where a stream leaves a puddle and returns to this puddle
        This test has to be done after computing the flow direction
        Loops are corrected after by a merge
    
        Returns
        -------
        puddleloop : set
            set of puddle id in which a stream returns.
    
        """
        # index of puddles from outlet, puddles flowing to virtual pits are excluded (no loop there)
        puddleidx = {self.puddledict[ipuddle]['outlet']:ipuddle for ipuddle in self.puddledict if self.puddledict[ipuddle]['outflow'] != -1}
        
        # list of springs
        springlist = [key for key in self.nodedict if self.isSpring(key)]
        # depth search of loops
        puddleloop = set()
        for ispring in springlist: 
            # take the downstream nodes and put them in a stack
            downstream = self.downstreamThalwegs(ispring)
            stack = []
            # ontheway stores outlet we go through
            ontheway = set()
            for i in downstream:
                stack.append(self.toNode(i))
            # we go downstream until reaching the end of a stream or finding a loop
            # a loop occurs if we pass twice through the same outlet
            while stack:
                inode = stack.pop()
                if inode == -1:
                    break
                # outlets are saddles
                typenode = self.nodedict[inode]['type']
                if typenode == df.SADDLE or typenode == df.CONFLUENCE:
                #if self.nodedict[inode]['type'] == df.SADDLE:
                    # check if we've already seen this node
                    if inode in ontheway:
                        puddleloop.add(puddleidx[inode])
                        break
                    # if the node is an outlet, add it to the set of visited outlets
                    if inode in puddleidx:
                        ontheway.add(inode)
                # keep going down
                downstream = self.downstreamThalwegs(inode)
                for i in downstream:
                    stack.append(self.toNode(i))
        return puddleloop

    def detectLoops(self):
        """
        Detects loops where a stream leaves a puddle and returns to this puddle
        This test has to be done after computing the flow direction
        Loops are corrected after by a merge
        The algortihm performs a DFS starting from puddle outlets and checks if
        a node is visited twice.
        Because a search starting from a puddle can loop in another puddle, the
        search starts from lower puddles so that their nodes are visited first
    
        Returns
        -------
        puddleloop : list
            list of puddle id in which a stream returns.
    
        """
        visited = {} # mark if a node has been visited, is being visited or not
        puddleloop = [] # list of puddles where loops were found
        for ip in self.nodedict:
            visited[ip] = -1 # mark all nodes as unvisited
    
        # get the dict of all outlets from where DFS will start
        # the dict is used to find back the pudlle where we are
        outletidx = {self.puddledict[ip]['outlet']: ip for ip in self.puddledict}
        # sort the outlets based on their elevation, starting from the lowest
        outletlist = sorted(outletidx.keys(), key = lambda i: self.nodedict[i]['z'])
        for ioutlet in outletlist:
            if visited[ioutlet] == -1:
                stack = [(0, ioutlet)] # mark as entered the outlet
                while stack:
                    act, inode = stack.pop()
                    if act == 1: # exit the node
                        visited[inode] = 1 # mark the node as visited
                        continue
                    visited[inode] = 0 # the node is being visited
                    stack.append((1, inode)) # add the node exit to the stack
                    # we go down to the next nodes
                    thalweglist = self.downstreamThalwegs(inode)
                    for it in thalweglist:
                        ipt = self.toNode(it)
                        # if this node is already being visited, we are looping
                        if visited[ipt] == 0: 
                            if ipt == ioutlet: # usually, we are still in the same puddle
                                puddleloop.append(outletidx[ipt])
                            else: # otherwise we have to look into other puddles
                                found = False
                                for ipuddle in self.puddledict:
                                    for i in self.puddledict[ipuddle]['nodes']:
                                        if ipt == i:
                                            found = True
                                            puddleloop.append(ipuddle)
                                            break
                                    if found:
                                        break
                                if not found:
                                    print("Found a loop at node", ipt, "but not the puddle")
                                    newpuddle = [ipt]
                                    node = self.nodedict[ipt]
                                    while node['type']!=df.PIT:
                                        downthalweg = [ii for ii in node['thalweg'] if ii>0]
                                        jpt = self.thalwegdict[downthalweg[0]]['end']
                                        newpuddle.append(jpt)
                                        node = self.nodedict[jpt]
                                    puddle, outlet, outflow = self.floodPuddle(newpuddle)
                                    self.puddledict[ipt] = {'nodes': puddle,
                                                             'outlet': outlet,
                                                             'outflow': outflow}
                                    puddleloop.append(ipt)
                        # if the node has not been visited yet, we add it in the stack as an entrance
                        elif visited[ipt] == -1:
                            stack.append((0, ipt))
        return puddleloop
    
    def extendPuddles(self, puddleloop):
        """
        Correct puddles where streams are looping.
        Correction done by adding the outflow node to the puddle to extend it
        It will then create a larger puddle that may be merged with others

        Parameters
        ----------
        puddleloop : set
            id of puddles to be extended.

        Returns
        -------
        None.

        """
        print(len(puddleloop), "boucles")
        for ipuddle in puddleloop:
            seed = self.puddledict[ipuddle]['nodes']
            ioutflow = self.puddledict[ipuddle]['outflow']
            seed.add(ioutflow)
            nodes, outlet, outflow = self.floodPuddle(seed)
            self.puddledict[ipuddle] = {'nodes': nodes,
                                        'outlet': outlet,
                                        'outflow': outflow}

    def mergePuddlesAndAssignFlow(self):
        """
        Loop on the puddle merging and the flow computation
        In each iteration, puddles are merged and the flow is directed
        If loops are detected, we iterate again

        Returns
        -------
        None.

        """
        self.mergePuddles()
        self.assignFlowDirection()
        puddleloop = self.detectLoops()
        while puddleloop:
            print("Correcting loops")
            self.extendPuddles(puddleloop)
            self.mergePuddles()
            self.assignFlowDirection()
            puddleloop = self.detectLoops()
    
    def makeSingleFlow(self):
        """
        Remove multiple direction flow by allowing no more than one descendant stream
        If a node has several downstreams, a new node is created at the same place
        and defined as a spring.

        Returns
        -------
        None.

        """
        nodelist = [i for i in self.nodedict if i >= 0 and self.nodedict[i]['type'] not in (df.PEAK, df.JUNCTION) and not self.isSpring(i)]
        for inode in nodelist:
            thalweglist = self.downstreamThalwegs(inode)
            if len(thalweglist) > 1: # we have multiple downstream flows
                slopelist = [self.computeStartSlope(it) for it in thalweglist]
                #smallslope = max(slopelist) # downward slope is negative
                # if all defluents go outside, we arbitrarily take the first one
                #if smallslope == df.NOTIN:
                #    continue
                bigslope = min(slopelist)
                i = slopelist.index(bigslope)
                del thalweglist[i] # this thalweg stays with inode
                # remove thalweglist from inode
                node = self.nodedict[inode]
                oldlist = [it for it in node['thalweg'] if abs(it) not in thalweglist]
                newlist = [it for it in node['thalweg'] if abs(it) in thalweglist]
                node['thalweg'] = oldlist
                # create a new node
                newkey = self.nodekey
                newnode = {'ij': node['ij'],
                           'z': node['z'],
                           'thalweg': newlist, # à mettre dans le bon sens
                           #'ridge': [inode], # pour compatibilité
                           'ridge': [], # pour compatibilité
                           'type': df.SPRING
                    }
                self.nodedict[newkey] = newnode
                # change the start/end node in newlist's thalwegs
                for it in newlist:
                    t = self.thalwegdict[abs(it)]
                    if t['start'] == inode:
                        t['start'] = newkey
                    elif t['end'] == inode:
                        t['end'] = newkey
                self.nodekey += 1

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

def getStreamsInBasin(self, thalweg, minorder = 1):
    thalwegstack = [thalweg]
    thalweglist = [thalweg]
    while thalwegstack:
        it = thalwegstack.pop()
        inode = self.toNode(it)
        order = self.nodedict[inode]['order']
        if self.thalwegdict[it]['order'] >= minorder and order > self.thalwegdict[it]['order']:
            thalweglist.append(it)
        inode = self.fromNode(it)
        upthalweg = self.upstreamThalwegs(inode)
        thalwegstack += upthalweg
    return thalweglist
        
def detectPuddleSinks(self):
    for ipuddle in self.puddledict:
        nodes = self.puddledict[ipuddle]['nodes']
        for inode in nodes:
            thalweg = self.nodedict[inode]['thalweg']
            b = True
            for t in thalweg:
                b = self.toNode(t) == inode
                if not b:
                    break
            if b and not self.isSink(inode):
                print(inode, "est un puit dans", ipuddle)
