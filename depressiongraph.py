# coding=utf-8
"""
depressiongraph.py:
    Class defining a graph of depressions.
    
    author: Eric Guilbert
"""
import math
import numpy as np
import deffunction as df
from celltype import CellType
from extendedsurfacenetwork import ExtendedSurfaceNetwork

class DepressionGraph:
    """
    A class defining a set of depressions connected into a graph
    A depression is built as an aggregation of dales coming from the ESN
    A depression is delineated by a set of ridges
    A depression is characterised by a pit and an outlet. 
    The pit is the lowest pit of the depression and the outlet is the 
    lowest saddle on the border of the depression.
    In the depression graph, nodes are depressions 
    and adjacency is defined by shared saddles so that an edge is defined
    by a saddle and two depressions sharing the saddle. Since a saddle can
    be connected to more than two depressions, the same saddle can be found
    in several edges. Two depressions can share several saddles so just the
    depressions is not enough to define the edge.
    """
    
    def __init__(self, network, culvert: "bool" = False):
        """
        

        Parameters
        ----------
        network : ExtendedSurfaceNetwork
            The ESN has been previously computed with its dales. The first step
            of the constructor is to compute the SN ascending cells that 
            represent the elementary depressions. After that, it computes the
            parameters of each depression and the edges of the graph.

        Returns
        -------
        None.

        """
        self.esn = network
        self.depdict = {}
        self.saddledict = {}
        self.edgedict = {}
        pitlist = [i for i in network.nodedict if network.nodedict[i]['type'] & CellType.PIT_FLAG]
        #pitlist.remove(-1)
        saddledict = {}
        for ipit in pitlist:
            # merge all the dales around a pit to form a depression 
            # take all the downslope lines (thalwegs or surface watercourses) connected to that pit
            thalweglist = network.getDownslopeLinesAtNode(ipit)
            culvertlist = []
            if culvert:
                # get the list of all underground watercourses connected to the pit
                line = network.nodedict[ipit]['line']
                culvertlist = [i for i in line if i > 0 and network.linedict[i]['type'] & CellType.UNDERGROUND_FLAG]
            thalwegset = {-j for j in thalweglist}
            stack = set(thalwegset)
            # go upward to get thalwegs above
            while stack:
                it = stack.pop()
                istart = network.linedict[it]['start']
                start = network.nodedict[istart]
                thalweglist = network.getDownslopeLinesAtNode(istart)
                if culvert:
                    line = network.nodedict[istart]['line']
                    culvertlist.extend([i for i in line if i > 0 and network.linedict[abs(i)]['type'] & CellType.UNDERGROUND_FLAG])
                tlist = [-i for i in thalweglist if i < 0] # thalwegs arriving at istart
                
                if start['type'] == CellType.SADDLE and tlist: # there is a thalweg arriving at the saddle
                    # we check if the thalwegs in tlist are in the same depression or not
                    # two thalwegs are in the same depression if they are not separated by a ridge
                    itpos = start['line'].index(it)
                    inextridge = network.nextRidgeIndexAtNode(istart, itpos)
                    ipreviousridge = network.previousRidgeIndexAtNode(istart, itpos)
                    other = []
                    if inextridge > ipreviousridge:
                        other = start['line'][ipreviousridge + 1: inextridge]
                    else:
                        other = start['line'][ipreviousridge +1:] + start['line'][:inextridge]
                    tlist = [t for t in tlist if -t in other]
                stack.update(tlist)
                thalwegset.update(tlist)
            # get ridge ring containing thalwegs
            plusridgeset, minusridgeset = network.getDivideSetBoundingThalwegs(thalwegset)
            minusridgeset = {-t for t in minusridgeset}
            self.depdict[ipit] = {
                'pit': ipit,
                'dale': thalwegset,
                'edge': set(),
                'outlet': None,
                'borderplus': plusridgeset,
                'borderminus': minusridgeset
            }
            self.computeOutlet(ipit)
            if culvert and culvertlist:
                outlet = self.depdict[ipit]['outlet']
                out = network.linedict[culvertlist[0]]['start']
                z = network.nodedict[out]['z']
                if z <= network.nodedict[outlet]['z']:
                    self.depdict[ipit]['outlet'] = out
                else:
                    self.depdict[ipit]['outlet'] = out
                    print("culvert higher than outlet", ipit, outlet, z)
            border = plusridgeset | minusridgeset
            nodeset = {network.linedict[ir]['start'] for ir in border}
            saddleset = {i for i in nodeset if network.nodedict[i]['type'] == CellType.SADDLE}
            for i in saddleset:
                if i in saddledict:
                    saddledict[i].append(ipit)
                else:
                    saddledict[i] = [ipit]
        # an edge connects two depressions and has the common saddle for attribute
        for i in saddledict:
            tmplist = sorted(saddledict[i], reverse = True)
            while tmplist:
                idep = tmplist.pop()
                for jdep in tmplist:
                    # if two deps are connected by several saddles, the edge is the lowest saddle
                    if (idep,jdep) in self.edgedict:
                        n1 = self.edgedict[idep,jdep]
                        z1 = network.nodedict[n1]['z']
                        z2 = network.nodedict[i]['z']
                        if z2 < z1:
                            self.edgedict[idep,jdep] = i
                        elif z2 == z1:
                            self.edgedict[idep, jdep] = min(i, n1)
                    else:
                        self.edgedict[idep,jdep] = i
        for ie in self.edgedict:
            idep, jdep = ie
            self.depdict[idep]['edge'].add(ie)
            self.depdict[jdep]['edge'].add(ie)
            i = self.edgedict[ie]
            if i in self.saddledict:
                self.saddledict[i].append(ie)
            else:
                self.saddledict[i] = [ie]

    def valid(self):
        linedict = self.esn.linedict
        for idep, dep in self.depdict.items():
            errorlist = [ir for ir in dep['borderplus'] if linedict[ir]['type'] != CellType.RIDGE]
            if errorlist:
                print(idep, errorlist)
        for idep, dep in self.depdict.items():
            errorlist = [ir for ir in dep['borderplus'] if linedict[ir]['type'] != CellType.RIDGE]
            if errorlist:
                print(idep, errorlist)

    def minDepth(self, iddep):
        """
        retuns the depth of a dep measured between its pit and its outlet

        Parameters
        ----------
        iddep : integer
            the id of the dep.

        Returns
        -------
        float
            the depth of the dep.

        """
        dep = self.depdict[iddep]
        node = self.esn.nodedict
        zmax = node[dep['outlet']]['z']
        zmin = node[dep['pit']]['z']
        return zmax - zmin
        
    def maxDepth(self, iddep):
        """
        returns the depth of a dep measured between its pit and the highest
        point on its boundary

        Parameters
        ----------
        iddep : integer
            the id of the dep.

        Returns
        -------
        float
            the depth of the dep.

        """
        dep = self.depdict[iddep]
        network = self.esn
        zmin = network.nodedict[dep['pit']]['z']
        ridges = dep['borderplus'] | dep['borderminus']
        nodeset = set()
        for ir in ridges:
            inode = network.ridgedict[ir]['end']
            nodeset.add(inode)
        # get the lowest node
        zmax = network.terrain.minz
        for inode in nodeset:
            z = network.nodedict[inode]['z']
            if z > zmax:
                zmax = z
        return zmax - zmin
        
    def getBorder(self, iddep):
        """
        Returns the ridges located on the border of the depression

        Parameters
        ----------
        iddep : integer
            Depression id.

        Returns
        -------
        set of ridge ids
            Ridges are returned in a set, meaning they are not ordered. 

        """
        dep = self.depdict[iddep]
        borderminus = dep['borderminus']
        return dep['borderplus'].union(borderminus)

    def getSprings(self, iddep):
        """
        Returns the set of nodes in the depression that are springs. A spring
        is a node from which all flows start and none arrives. Springs can be
        on the border or inside the depression

        Parameters
        ----------
        iddep : integer
            id of the depression.

        Returns
        -------
        node : set of integers
            The ids of the nodes that are springs.

        """
        network = self.esn
        # get all the ends of the thalwegs
        dep = self.depdict[iddep]
        node = {network.fromNode(t) for t in dep['dale']}
        node = {i for i in node if network.isSpring(i)}
        return node
            
    def getLowerDep(self, iddep):
        outlet = self.depdict[iddep]['outlet']
        edges = self.saddledict[outlet]
        dep = {i for edge in edges for i in edge if i != iddep}
        lowerdep = [d for d in dep if self.depdict[d]['outlet'] != outlet]
        if not lowerdep:
            return None
        return lowerdep[0]
        
    def computeOutlet(self, iddep):
        """
        The outlet of a depression is the lowest saddle on its border

        Parameters
        ----------
        iddep : integer
            id of the depression.

        Returns
        -------
        idmin : integer
            id of the outlet.

        """
        network = self.esn
        # get all the saddles located on the prominence border
        ridges = self.getBorder(iddep)
        saddleset = set()
        for ir in ridges:
            if ir < 0:
                print("erreur", ir, iddep)
            isaddle = network.linedict[ir]['start']
            saddleset.add(isaddle)
        # get the lowest saddle
        idmin = -1
        zmin = network.terrain.maxz
        #ijmin = (network.terrain.m,network.terrain.n)
        for isaddle in saddleset:
            saddle = network.nodedict[isaddle]
            z = saddle['z']
            if z < zmin:
                idmin = isaddle
                zmin = z
                #ijmin = saddle['ij']
            elif z == zmin:
                #if saddle['ij'][0] < ijmin[0]:
                #    idmin = isaddle
                #    ijmin = saddle['ij']
                #elif saddle['ij'][0] == ijmin[0] and saddle['ij'][1] < ijmin[1]:
                #    idmin = isaddle
                #    ijmin = saddle['ij']
                idmin = min(idmin, isaddle)
                
        if idmin == -1:
            print("error outlet", iddep)
        self.depdict[iddep]['outlet'] = idmin
         
    def mergeDepressions(self, mindep, idep):
        """
        Merge two depressions into one

        Parameters
        ----------
        mindep : integer
            id of the depression that is kept.
        idep : integer
            id of the depression merged to mindep.

        Returns
        -------
        None.

        """
        network = self.esn
        # get all other edges of idep
        edge = min(mindep, idep), max(mindep, idep)
        otheredges = set(self.depdict[idep]['edge'])
        otheredges.remove(edge)
        # update other edges by replacing other dep values with mindep
        otherdep = -1
        for e in otheredges:
            # get the other dep connected to mindep
            if e[0] == idep:
                otherdep = e[1]
            else:
                otherdep = e[0]
            newe = min(mindep, otherdep), max(mindep, otherdep)
            olds = self.edgedict[e] # saddle between idep and otherdep
            # newe replaces e
            # if newe already exists, check which saddle to keep
            if newe in self.edgedict:
                j1 = self.edgedict[newe] # saddle between mindep and otherdep
                z1 = network.nodedict[j1]['z']
                z2 = network.nodedict[olds]['z']
                if z2 < z1 or (z2 == z1 and olds < j1):
                    self.edgedict[newe] = olds # the edge is the older saddle
                    self.saddledict[olds].append(newe)
                    self.saddledict[j1].remove(newe)
            else:
                self.edgedict[newe] = olds
                self.saddledict[olds].append(newe)
            self.saddledict[olds].remove(e)
            # the old edge is removed
            del self.edgedict[e]
            # replace the edge in the dep
            self.depdict[otherdep]['edge'].remove(e)
            self.depdict[otherdep]['edge'].add(newe)
            # and add it in mindep
            self.depdict[mindep]['edge'].add(newe)

        self.depdict[mindep]['edge'].remove(edge)
        saddle = self.edgedict[edge]
        self.saddledict[saddle].remove(edge)
        del self.edgedict[edge]

        # update dep border and outlet
        superdep = self.depdict[mindep]
        borderplus = superdep['borderplus']
        borderminus = superdep['borderminus']
        # add all the borders
        dep = self.depdict[idep]
        borderplus.update(dep['borderplus'])
        borderminus.update(dep['borderminus'])
        # remove those which also appear on the other side
        newborderminus = borderminus.difference(borderplus)
        newborderplus = borderplus.difference(borderminus)
        superdep['borderplus'] = newborderplus
        superdep['borderminus'] = newborderminus
        # update the dale set
        superdep['dale'].update(dep['dale'])
        
    def buildPuddleGraph(self):
        """
        Merges depressions of the graph sharing the same outlet
        The process is repeated until no more outlet is shared
        The resulting depressions correspond to puddles where the water can
        accumulate and leave through the outlet towards another depression

        Returns
        -------
        None.

        """
        network = self.esn
        # order saddles from the lowest to the highest
        sortedsaddles = [i for i in self.saddledict]
        nsaddle = len(sortedsaddles)
        nmax = nsaddle
        nsup = nsaddle + 1
        sortedsaddles.sort(key = lambda i: network.nodedict[i]['z'])
        while nsup > nsaddle:
            nsup = nsaddle
            if nsup/nmax < 0.8: # try to optimise by not sorting too often
                sortedsaddles = [i for i in self.saddledict]
                sortedsaddles.sort(key = lambda i: network.nodedict[i]['z'])
                nmax = nsup
            print("number of saddles", nsaddle)
            nsaddle = len(sortedsaddles)
            for isaddle in sortedsaddles:
                if not isaddle in self.saddledict:
                    nsaddle -= 1
                    continue
                edges = self.saddledict[isaddle]
                if not edges:
                    continue
                # get the depressions connected to this saddle
                deps = {i for edge in edges for i in edge}
                # keep those whose outlet is the saddle
                deps = {i for i in deps if self.depdict[i]['outlet'] == isaddle}
                if len(deps) <=1:
                    continue
                # two depressions have the same outlet
                # keep the lowest one
                zmin = network.terrain.maxz
                mindep = -1
                # get the lowest depression
                for dep in deps:
                    ipit = self.depdict[dep]['pit']
                    z = network.nodedict[ipit]['z']
                    if z < zmin:
                        zmin = z
                        mindep = dep
#                if len(edges) != 1:
#                    print(": Edges at saddle", isaddle, "are", edges)
                if mindep == -1:
                    continue
                deps.remove(mindep)
                for dep in deps:
                    self.mergeDepressions(mindep, dep)
                    del self.depdict[dep]
                    nsaddle -= 1
                self.computeOutlet(mindep)
            klist = list(self.saddledict.keys())
            for k in klist:
                if not self.saddledict[k]:
                    del self.saddledict[k]
#            sortedsaddles = [i for i in self.saddledict]
#            nsaddle = len(sortedsaddles)

#        klist = list(self.edgedict.keys())
#        for k in klist:
#            if not self.edgedict[k]:
#                del self.edgedict[k]
        
    def buildPuddleGraph2(self):
        """
        Merges depressions of the graph sharing the same outlet
        The process is repeated until no more outlet is shared
        The resulting depressions correspond to puddles where the water can
        accumulate and leave through the outlet towards another depression

        Returns
        -------
        None.

        """
        network = self.esn
        # order saddles from the highest to the lowest
        sortedsaddles = [i for i in self.saddledict]
        nsaddle = len(sortedsaddles)
        nsup = nsaddle + 1
        while nsup > nsaddle:
            nsup = nsaddle
            print("number of saddles", nsaddle)
            # sorting should be done before
            sortedsaddles.sort(key = lambda i: network.nodedict[i]['z'], reverse = True)
            while sortedsaddles:
                isaddle = sortedsaddles.pop() # à réécrire pour éviter de refaire un sort
                edges = self.saddledict[isaddle]
                if not edges:
                    continue
                # get the depressions connected to this saddle
                deps = {i for edge in edges for i in edge}
                # keep those whose outlet is the saddle
                deps = {i for i in deps if self.depdict[i]['outlet'] == isaddle}
                if len(deps) <=1:
                    continue
                # two depressions have the same outlet
                # keep the lowest one
                zmin = network.terrain.maxz
                mindep = -1
                # get the lowest depression
                for dep in deps:
                    ipit = self.depdict[dep]['pit']
                    z = network.nodedict[ipit]['z']
                    if z < zmin:
                        zmin = z
                        mindep = dep
#                if len(edges) != 1:
#                    print(": Edges at saddle", isaddle, "are", edges)
                deps.remove(mindep)
                for dep in deps:
                    self.mergeDepressions(mindep, dep)
                    del self.depdict[dep]
                self.computeOutlet(mindep)
            klist = list(self.saddledict.keys())
            for k in klist:
                if not self.saddledict[k]:
                    del self.saddledict[k]
            sortedsaddles = [i for i in self.saddledict]
            nsaddle = len(sortedsaddles)

        klist = list(self.edgedict.keys())
        for k in klist:
            if not self.edgedict[k]:
                del self.edgedict[k]

    def mergeFlatDepressions(self, precision):
        """
        À faire: agrégation des dépressions en fonction d’une précision
        La précision peut être une différence de hauteur, une surface ou un volume

        Parameters
        ----------
        precision : TYPE
            DESCRIPTION.

        Returns
        -------
        None.

        """
        network = self.esn
        # get depressions flatter than precision
        sorteddeps = [i for i in self.depdict.keys() if self.maxDepth(i) < precision]
        #sorteddeps = sorted(self.depdict.keys(), key = lambda i: self.minDepth(i), reverse = True)
        ndep = len(sorteddeps)
        nsup = ndep + 1
        #while nsup > ndep:
        removed = set()
        sorteddeps.sort(key = lambda i: self.maxDepth(i), reverse = True)
        print(len(sorteddeps))
        while sorteddeps:
            nsup = ndep
            mindep = sorteddeps.pop() # this is the flatest depression
            if mindep == 545479:
                print(mindep, self.depdict[mindep]['edge'])
            if mindep in removed:
                continue
            if self.maxDepth(mindep) >= precision:
                continue

            dep = self.depdict[mindep]
            # get edges of dep
            edges = dep['edge']
            # get neighbouring depressions
            deps = set(i for e in edges for i in e)
            # bogue mindep n’est pas dans deps
            deps.remove(mindep)
            if not deps:
                continue
            # get the flatest depression
            idep = min(deps, key = self.maxDepth)
            if network.nodedict[dep['pit']]['z'] > network.nodedict[self.depdict[idep]['pit']]['z']:
                (idep, mindep) = (mindep, idep)
            self.mergeDepressions(mindep, idep)
            #print(len(sorteddeps), mindep, idep)
            if mindep == 545479 or idep == 545479:
                print(mindep, idep)
            removed.add(idep)
            del self.depdict[idep]
            self.computeOutlet(mindep)
            
            klist = list(self.saddledict.keys())
            for k in klist:
                if not self.saddledict[k]:
                    del self.saddledict[k]

        klist = list(self.edgedict.keys())
        for k in klist:
            if not self.edgedict[k]:
                del self.edgedict[k]
            
    def mergeFlatPuddles(self, precision):
        """
        À faire: agrégation des dépressions en fonction d’une précision
        La précision peut être une différence de hauteur, une surface ou un volume

        Parameters
        ----------
        precision : TYPE
            DESCRIPTION.

        Returns
        -------
        None.

        """
        network = self.esn
        # get saddles ordered in z
        saddles = set(self.edgedict.values())
        sortedsaddles = sorted(saddles, key = lambda i: network.nodedict[i]['z'], reverse = True)
        print(len(sortedsaddles))
        while sortedsaddles:
            isaddle = sortedsaddles.pop() # get the lowest saddle
            zmax = network.nodedict[isaddle]['z']
            edges = self.saddledict[isaddle] # get edges connected to this saddle
            # get the depressions connected to the edge with a flat bottom
            deps = set(i for e in edges for i in e if zmax - network.nodedict[i]['z'] < precision)
            if not deps:
                continue
            # get all the edges that connect these deps
            edges = [e for e in edges if e[0] in deps and e[1] in deps]
            # sort deps according to depth
            sorteddeps = sorted(deps, key = lambda i: network.nodedict[i]['z'])
            # merge the deps along these edges starting from the highest
            idep = sorteddeps.pop()
            while sorteddeps:
                tmpdeps = sorteddeps[:]
                while tmpdeps:
                    mindep = tmpdeps.pop()
                    if (mindep, idep) in edges or (idep, mindep) in edges:
                        print(len(sorteddeps), mindep, idep, isaddle)
                        self.mergeDepressions(mindep, idep)
                        tmpdeps = []
                        del self.depdict[idep]
                        idep = sorteddeps.pop()
            self.computeOutlet(mindep)
            
            klist = list(self.saddledict.keys())
            for k in klist:
                if not self.saddledict[k]:
                    del self.saddledict[k]

        klist = list(self.edgedict.keys())
        for k in klist:
            if not self.edgedict[k]:
                del self.edgedict[k]
            
    def computeBorder(self, depid):
        """
        Compute the border of a depression. Can be used to update the value.
        Ridges in the border are stored in two sets depending on their orientation

        Parameters
        ----------
        depid : integer
            Depression id.

        Returns
        -------
        None.

        """
        network = self.esn
        dep = self.depdict[depid]
        dales = dep['dale']
        stack = set(dep['lower'])
        while stack:
            idep = stack.pop()
            dep = self.depdict[idep]
            dales.update(dep['dale'])
            stack.update(dep['lower'])
        # dales contains all the dales inside the depression
        # get all the ridges in the depression 
        # and keep only those on the border of the depression
        ridges = []
        for idale in dales:
            dale = network.daledict[idale]
            ridges.extend(dale['boundary'])
            for hole in dale['hole']:
                ridges.extend(hole)
        # now we have all the ridges located in the depression
        # to keep those on the border, we keep only those which apear once
        # others appear twice as positive and negative
        ridgeplus = {t for t in ridges if t > 0}
        ridgeminus = {t for t in ridges if t < 0}
        borderplus = {t for t in ridgeplus if -t not in ridgeminus}
        borderminus = {-t for t in ridgeminus if -t not in ridgeplus}
        self.depdict[depid]['borderplus'] = borderplus
        self.depdict[depid]['borderminus'] = borderminus

    def getThalwegsBelowOutlet(self, depid):
        """
        Return thalwegs in a depression located below the outlet.
        A thalweg is below the outlet if its highest point is strictly lower
        than the outlet or if it is connected to the outlet.

        Parameters
        ----------
        depid : integer
            id of the depression.

        Returns
        -------
        thalwegs : set
            The set of thalweg ids.

        """
        network = self.esn
        dep = self.depdict[depid]
        ioutlet = dep['outlet']
        outlet = network.nodedict[ioutlet]
        z = outlet['z']
        dales = dep['dale']
        thalwegs = network.getDownslopeLinesAtNode(ioutlet)
        thalwegs = {abs(it) for it in thalwegs if it in dales}
        for it in dales:
            t = network.linedict[it]
            istart = t['start']
            start = network.nodedict[istart]
            if start['z'] <= z:
                thalwegs.add(it)
        return thalwegs

    def getThalwegsInsidePuddles(self):
        """
        Return the set of thalwegs below the outlet in each depression
        The result is returned in a dictionary where the key is the depression id

        Returns
        -------
        thalwegdict : dict
            A dictionary where the keys are the depression ids and the values 
            are sets of thalwegs below the outlet in each depression.

        """
        thalwegdict = {}
        for i in self.depdict:
            thalwegs = self.getThalwegsBelowOutlet(i)
            thalwegdict[i] = (self.depdict[i]['outlet'], thalwegs)
        return thalwegdict
    
    def createDepMatrix(self):
        network = self.esn
        terrain = network.terrain
        ridgedict = network.linedict
        boundary = set(terrain.boundary)
        self.depmatrix = terrain.nodata*np.ones((terrain.m, terrain.n), dtype=int)
        
        for idep in self.depdict:
            if idep == -1:
                continue
            dep = self.depdict[idep]
            ridges = [r for r in dep['borderminus']]
            ridges.extend(dep['borderplus'])
            borderset = set()
            for r in ridges:
                borderset.update(ridgedict[r]['polyline'])
            pointset = set()
            for t in dep['dale']:
                line = network.linedict[t]['polyline'][1:]
                pointset.update(line)
            pointset.difference_update(borderset)
            stack = set(pointset)
            while stack:
                i,j = stack.pop()
                neighbours = [(i+di, j+dj) for di,dj in df.ldr4]
                neighbours = [p for p in neighbours if p not in pointset and p not in borderset and p not in boundary]
                stack.update(neighbours)
                pointset.update(neighbours)
            for (i,j) in pointset:
                self.depmatrix[i,j] = idep
        
    def createPuddleMatrix(self):
        network = self.esn
        terrain = network.terrain
        zdict = {ip: network.nodedict[self.depdict[ip]['outlet']]['z'] for ip in self.depdict}
        self.puddlematrix = self.depmatrix.copy()
        it = np.nditer(self.puddlematrix, flags = ['multi_index'], op_flags = ['readwrite'])
        for x in it:
            (i, j) = it.multi_index
            ip = x.item()
            if ip !=terrain.nodata and terrain.dtm[i,j] > zdict[ip]:
                x[...] = terrain.nodata
