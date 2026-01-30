# -*- coding: utf-8 -*-
"""
Created on Tue Oct  1 23:44:12 2023

@author: Yassmine Zada
"""

from time import perf_counter_ns
from extendedsurfacenetwork import ExtendedSurfaceNetwork
import geopandas as gpd
from shapely.geometry import LineString
from shapely.wkt import loads
import deffunction as df
from math import sqrt
from pyproj import Transformer
from terrain import bresenham


class TopoHydroNetwork(ExtendedSurfaceNetwork):

        
    def buildTopoHydroFromTerrain(self, terrain, culvert, sd = True):
        """
        Builds a surface network on a given terrain. Start by detecting saddles.
        Follows with thalwegs. Then integrates culverts and corrects thalweg. 
        Follows by ridges. Ridges are constrained to avoid 
        conflicts with thalwegs.

        Parameters
        ----------
        terrain : Terrain
            Instance of the class Terrain on which the network is built.
        correction : boolean
            True if a correction is done to stay closer to the gradient
            False if the descent direction is not corrected as in D8
        fichier :
            path to culvert file
        
        Returns
        -------
        None.

        """
        self.culvertdict = {}
        self.terrain = terrain
        self.culvert = culvert
        
        self.buildFlowThalwegsFromTerrain(terrain, correction = False)
        start = perf_counter_ns()
        self.addBifurcation()
        end = perf_counter_ns()
        print("Bifurcation computation time:", end - start)
        start = perf_counter_ns()
        self.computeCulverts()
        end = perf_counter_ns()
        print("Culverts computation time:", end - start)
        start = perf_counter_ns()
        self.correctThalwegs()
        end = perf_counter_ns()
        print("Correction Thalweg time:", end - start)
        start = perf_counter_ns()
        self.addNewThalwegs(correction=True)
        end = perf_counter_ns()
        print("New Thalweg computation time:", end - start)
        self.orderThalwegsAroundNodes()
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
        start = perf_counter_ns()
        divide = True
        print("Build ridges")
        self.singleflowdirection = sd
        smooth = False
        if smooth:
            df.distancefactor = 2/3
            correction = True
        else:
            df.distancefactor = sqrt(2)/2
            correction = False
        self.instantiateRidgesAtSaddles(divide)
        # erase line property in saddles
        for ipt in self.nodedict:
            if 'line' in self.nodedict[ipt]:
                del self.nodedict[ipt]['line']
        self.instantiateRidgesAtPitsAndConfluencesAndTransfluences()
        self.traceRidgesT(divide, correction)
        self.orderRidgesAroundNodesC()
        end = perf_counter_ns()
        print("Ridge computation time:", end - start)
        start = perf_counter_ns()
        self.computeThalwegDales()
        end = perf_counter_ns()
        print("Dale computation time:", end - start)
        start = perf_counter_ns()
        self.flowcomputedict()
        self.addPuddlesC()
        self.mergePuddlesAndAssignFlowC()
        end = perf_counter_ns()
        print("Flow direction computation time:", end - start)
        start = perf_counter_ns()
        # self.flowcomputedict()
        self.accumulateFlowC()
        end = perf_counter_ns()
        print("Flow accumulation computation time:", end - start)
        start = perf_counter_ns()
        self.reassignmentnode()
        end = perf_counter_ns()
        print("Reassignment nodes time:", end - start)
    
    def flowcomputedict(self):
        """
        calculates the flow dictionary, including surface flows and flows through culverts

        Parameters
        ----------
        None. 
        
        Returns
        -------
        None.

        """
        self.flowdict = self.thalwegdict.copy()
        max_thalweg_key = max(self.thalwegdict.keys())
        for key in self.flowdict:
            self.flowdict[key]['type'] = 'thalweg'
        for key, culvert_info in self.culvertdict.items():
            max_thalweg_key += 1
            self.flowdict[max_thalweg_key] = culvert_info
            self.flowdict[max_thalweg_key]['type'] = 'culvert'
            self.flowdict[max_thalweg_key]['dale'] = -999
            
    def reassignmentnode(self) :
        """
        reassigns the node type according to the flow direction 

        Parameters
        ----------
        None. 
        
        Returns
        -------
        None.

        """
        for ipt, pt in self.nodedict.items():
            downstream = self.downstreamThalwegsC(ipt)
            downlength = len(downstream)
            upstream = self.upstreamThalwegsC(ipt)
            uplength = len(upstream)
            if uplength == 0 and downlength > 0 :
                self.nodedict[ipt]['fluence'] = 1 ###Spring
            elif uplength > 0 and downlength > 1 :
                self.nodedict[ipt]['fluence'] = 2 ###Diffluence
            elif uplength > 1 and downlength == 1 :
                self.nodedict[ipt]['fluence'] = 3 ###Confluence
            elif uplength == 1 and downlength == 1 :
                self.nodedict[ipt]['fluence'] = 4 ###Transfluence
            elif uplength == 1 and downlength == 1 :
                self.nodedict[ipt]['fluence'] = 5 ###Sink (Doesn't exist in our case)
            else :
                self.nodedict[ipt]['fluence'] = 6 ###Other (Peak and RJunction)
        
        
    def computeBifurcation(self) :
        """
        Calculates connection points between thalwegs and culverts 

        Parameters
        ----------
        fichier :
            path to culvert file
            
        Returns
        -------
        None.

        """
        
        lignes_ponceaux = self.culvert[0]
        ponceaux = self.terrain.pixelCulverts(lignes_ponceaux)   ##Pixelation of culverts
        ponceaux_avec_z = ponceaux[0]
        ponceaux_avec_z_prolonges = ponceaux[1]

        T = self.InOutCulverts(ponceaux_avec_z, ponceaux_avec_z_prolonges)  ##Calculation of culvert inlets and outlets
        entrees_candidates = T[1]
        sorties_candidates = T[2]
        
        T1 = self.entrees_sorties_sansZ(entrees_candidates, sorties_candidates)
        entrees_candidatess = T1[0]
        sorties_candidatess = T1[1]

        T2 = self.entrees_sorties_finaux(entrees_candidatess, sorties_candidatess, sorties_candidates)
        bifurcations_entrees = T2[1]
        bifurcations_sorties = T2[3]
        
        bifurcations =  []
        for i in range(len(bifurcations_entrees)) :
            bifurcations.append(bifurcations_entrees[i])
            bifurcations.append(bifurcations_sorties[i])
        
        return bifurcations_entrees, bifurcations_sorties, bifurcations

    def addBifurcation(self):
        """
        Adds connection points between thalwegs and culverts to the dictionary
        of nodes and assigns them a specific type 

        Parameters
        ----------
        fichier :
            path to culvert file
            
        Returns
        -------
        None.

        """
        
        # matrix dimension 
        A = self.terrain.dtm
        pixelclass = self.terrain.pixelclass
        
        self.bifurcationdict = {}  # dictionary of bifurcations
        bifurcationidx = {}  # dictionary indexing the coordinates of the bifurcations
        
        B = self.computeBifurcation()
        bifurcations = B[2]
        cles = list(self.nodedict.keys())  # retrieve dictionary keys
        nodekey = cles[-1]    # retrieve the index of the last point in nodedict

        count = 0
        for point in bifurcations :
            i = point[0]
            j = point[1]
            z = A[i, j]
        
            self.bifurcationdict[count] = {
                'ij': (i, j),
                'z': z,
                'thalweg': [],
                'ridge': [],
                'culvert': [],
                'type': df.BIFURCATION
                }
            bifurcationidx[i, j] = count
            if pixelclass[i, j] == df.SADDLE or pixelclass[i, j] == df.PIT or pixelclass[i, j] == df.CONFLUENCE :
                key = self.nodeidx[(i,j)]   # retrieve the index of this node
                self.nodedict[key]['culvert'] = []   # add a “culvert” field to that node

            else :
                
                self.nodeidx[(i,j)] = self.nodekey  # add to nodeidx the following indices at the bifurcations

                # add the bifurcations to the nodedict
                self.nodedict[self.nodekey] = {
                    'ij': (i, j),
                    'z': z,
                    'thalweg': [],
                    'ridge': [],
                    'culvert': [],  # Add the “culvert” field here and initialise it with an empty list.
                    'type': df.BIFURCATION
                }
                self.terrain.pixelclass2[i, j] = df.BIFURCATION   #The reason why we don't add it in pixelclass is to preserve the fact that it is a thalweg dividing in two
                # increment the index
                self.nodekey += 1
            count += 1
  

    def correctThalwegs(self) :
        """
        Allows you to correct thalwegs that have been split in two by culverts 
        by inserting a node and then having two thalwegs instead of one 

        Parameters
        ----------
        fichier :
            path to culvert file
    

        Returns
        -------
        None.

        """
        
        cles = list(self.thalwegdict.keys())  # retrieve dictionary keys
        cle = cles[-1]    # retrieve the index of the last point in nodedict
        nextthalweg = cle + 1  
        
        for k in self.bifurcationdict:
            
            # retrieve the bifurcations indices
            i = self.bifurcationdict[k]['ij'][0]
            j = self.bifurcationdict[k]['ij'][1]
            
            if self.terrain.pixelclass[i,j] == df.THALWEG :
                
                idthalweg = self.floorkey[i,j]
                endindex = self.nodeidx[(i,j)] 
                self.nodedict[endindex]['thalweg'].append(-idthalweg)
                self.nodedict[endindex]['thalweg'].append(nextthalweg)

                # split the existing thalweg in two
                t = self.thalwegdict[idthalweg]
                ptlist = t['polyline']
                iend = t['end']
                icfl = ptlist.index((i,j))
                ptlist1 = ptlist[:icfl+1]
                ptlist2 = ptlist[icfl:]
                # modify the first thalweg to stop at the confluence
                t['polyline'] = ptlist1
                t['end'] = endindex   #bifurcation serait la fin
                # create a new thalweg
                self.thalwegdict[nextthalweg] = {
                    'start' : endindex,
                    'polyline' : ptlist2,
                    'end' : iend,
                    'lefthill' : [],
                    'righthill' : [],
                }
                
                # we update the node at the end of the new thalweg
                ndt = self.nodedict[iend]['thalweg']
                ndt.remove(-idthalweg)
                ndt.append(-nextthalweg)
                # we update floorkey with nextthalweg
                for k in ptlist2[1:-1]:
                    self.floorkey[k] = nextthalweg
                nextthalweg += 1
                self.terrain.pixelclass[i,j] = df.BIFURCATION

    
    def computeCulverts(self) :
        """
        Calculates the new list of culverts with the new entry and exit points

        Parameters
        ----------
        fichier :
            path to culvert file

        Returns
        -------
        None.

        """
        
        culvertkey = 1
        # floorkey={}
        B = self.computeBifurcation()
        dict_ponceaux = self.culvert[1]
        b1 = B[0]
        b2 = B[1]
        L = self.newCulverts(b1, b2)
        culverts = L[0]
        for culvertt in culverts :
            polylist = []    # initiate the list of polylists to be added to the culvert dictionary
            # extract the i,j from the inlets and outlets
            i_in = culvertt[0][0]      
            j_in = culvertt[0][1]
            i_out = culvertt[-1][0] 
            j_out = culvertt[-1][1]
            # floorkey[i_in, j_in] = culvertkey
            
            # loop to fill the polylist with its corresponding values
            for point in culvertt :
                i = point[0]
                j = point[1]
                # floorkey[i,j] = culvertkey
                polylist.append((i,j))
            
            # retrieve the indices of the inputs as startindex and the indices of the outputs as endindex
            startindex = self.nodeidx[(i_in,j_in)]
            endindex = self.nodeidx[(i_out,j_out)]
            
            # store culvert indices in nodedict for points that allow water to pass through culverts
            self.nodedict[startindex]['culvert'].append(culvertkey)
            self.nodedict[endindex]['culvert'].append(-culvertkey)
            
            # create the culvert dictionary
            self.culvertdict[culvertkey]={
                'start' : startindex,
                'polyline' : polylist,
                'end' : endindex,
            }
            culvertkey += 1   
       # for key in self.culvertdict.keys() & dict_ponceaux.keys():
            #self.culvertdict[key]['diametre'] = dict_ponceaux[key]['DIAMETRE']
        
    
    def addNewThalwegs(self, correction=True):
        
        bifurcation_nn_conn_th = []
        
        for k in self.bifurcationdict:
            
            # retrieve the bifurcations indices
            i = self.bifurcationdict[k]['ij'][0]
            j = self.bifurcationdict[k]['ij'][1]
            if self.terrain.pixelclass2[i,j] == df.BIFURCATION:
                key = self.nodeidx[(i, j)]
                if not self.nodedict[key]['thalweg'] :
                    bifurcation_nn_conn_th.append(key+1)

        print("bifurcation_nn_conn_th = ", bifurcation_nn_conn_th)
       
        cles = list(self.nodedict.keys())  # retrieve dictionary keys
        nodekey = cles[-1]    # retrieve the index of the last point in nodedict
        
        cles = list(self.thalwegdict.keys())  # retrieve dictionary keys
        cle = cles[-1]    # retrieve the index of the last point in nodedict
        nextthalwegkey = cle + 1
        for indice in bifurcation_nn_conn_th:
            pt = self.nodedict[indice-1]
            (i,j) = self.nodedict[indice-1]['ij']
            polylist = [(i,j)]
            keepgoing = True # no node yet
            outbound = False # not out of the terrain
            # endindex = -1
            # check if we are out of the domain or on a node
            if self.terrain.isInside(i,j):
                if self.terrain.pixelclass[i,j] == df.SLOPE:
                    keepgoing = True
                    self.terrain.pixelclass[i,j] = df.THALWEG
                    self.floorkey[i,j] = nextthalwegkey
                else:
                    keepgoing = False
            else:
                keepgoing = False
                outbound = True # we’ll connect the thalweg to the virtual pit
            # if not a node and inbound, we look for the next point
            self.terrain.pixelclass[i,j] = df.BIFURCATION ##ADDED 20231202
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
                        det = 1/(x1*y2 - x2*y1)
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
                            else: lag = 0
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
               
                # check again if the point is a node or out
                if self.terrain.isInside(i,j):
                    if self.terrain.pixelclass[i,j] == df.SLOPE:
                        keepgoing = True
                        self.terrain.pixelclass[i,j] = df.THALWEG
                        self.floorkey[i,j] = nextthalwegkey
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
                    # endindex = nodekey + 1
                    # self.nodeidx[(i,j)] = endindex
                    self.nodedict[self.nodeidx[(i,j)]] = {
                        'ij' : (i,j),
                        'z': self.terrain.dtm[i,j],
                        'thalweg' : [-idthalweg, nextthalwegkey],
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
                    self.thalwegdict[nextthalwegkey] = {
                        'start' : endindex,
                        'polyline' : ptlist2,
                        'end' : iend,
                        'lefthill' : [],
                        'righthill' : [],
                    }
                    # we update the node at the end of the new thalweg
                    ndt = self.nodedict[iend]['thalweg']
                    ndt.remove(-idthalweg)
                    ndt.append(-nextthalwegkey)
                    # we update floorkey with nextthalweg
                    for k in ptlist2[1:-1]:
                        self.floorkey[k] = nextthalwegkey
                    nextthalwegkey += 1
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
            pt['thalweg'].append(nextthalwegkey)
            self.nodedict[endindex]['thalweg'].append(-nextthalwegkey)
            self.thalwegdict[nextthalwegkey]={
                'start' : indice-1,
                'polyline' : polylist,
                'end' : endindex,
                'lefthill' : [],
                'righthill' : [],
            }
            nextthalwegkey += 1
            
    def getThalwegNeighbourNodesC(self, ipt):
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
        if 'culvert' in pt:
            for it in pt['culvert'] :
                if it > 0 and it is not None :
                    lowerneighbour.add(self.culvertdict[it]['end'])
                elif it is not None:
                    higherneighbour.add(self.culvertdict[-it]['start']) 
            
        return lowerneighbour, higherneighbour
    
    def accumulateFlowC(self):
        """
        computes the flow accumulation at the mouth of each thalweg. The flow
        accumulation is defined by the area of all dales upstream of the mouth.
        If a node has severalboolean downstream thalwegs, the accumulation can be sent
        to the thalwegs in proportion of the slope or to only one thalweg.

        Parameters
        ----------
        sd : 
            True if the network does not have diffluences (single direction).
            In that case, all the flow follows the steepest slope.
            Otherwise, the flow to each thalweg depends on the slope

        Returns
        -------
        None.

        """
        sd = self.singleflowdirection
        visited = {}
        for it in self.flowdict:
            t = self.flowdict[it]
            idale = t['dale']
            visited[it] = len(self.upstreamThalwegsC(self.fromNodeC(it)))
            if idale != -1 and idale != -999:
                t['accumulation'] = self.daledict[idale]['area']
            else:
                t['accumulation'] = 0

        streamlist = [self.downstreamThalwegsC(ipt) for ipt in self.nodedict if self.isSpringC(ipt)]
        thalweglist = [it for sublist in streamlist for it in sublist]
        for ist in thalweglist:
            stack = [ist]
            while stack:
                it = stack.pop()
                inode = self.toNodeC(it)
                downstream = self.downstreamThalwegsC(inode)
                if downstream:
                    n = len(downstream)
                    slopelist = [min(0, self.computeStreamSlopeC(idt))
                                 for idt in downstream]
                    if df.NOTIN in slopelist:
                        for i in range(n):
                            if slopelist[i] != df.NOTIN:
                                slopelist[i] = 0
                    accu = [0]*n
                    accutot = self.flowdict[it]['accumulation']
                    if sd:
                        maxslope = min(slopelist)
                        imax = slopelist.index(maxslope)
                        for i in range(n):
                            lt = downstream[i]
                            if self.flowdict[lt]['type'] == 'thalweg':
                                nd = self.flowdict[lt]['start']
                                
                            if self.flowdict[lt]['type'] == 'culvert':
                                accu[i] = accutot
                            # elif i == imax and 'culvert' not in self.nodedict[nd]:
                            elif i == imax :
                                accu[imax] = accutot
                            else:
                                accu[i] = 0
                                jt = downstream[i]
                                self.flowdict[jt]['headwater'] = True
                    for i in range(n):
                        idt = downstream[i]
                        self.flowdict[idt]['accumulation'] += accu[i]
                        visited[idt] -= 1
                        if visited[idt] == 0:
                            stack.append(idt)
                            
    def instantiateRidgesAtPitsAndConfluencesAndTransfluences(self):
        """
        Instantiates ridge objects between two thalwegs at pits of degree three
        or more. 

        Returns
        -------
        None.

        """
        # key for new ridges
        ridgestart = max(self.ridgedict) + 1  # REMARQUE
        ridgekey = ridgestart
        divideside = {}

        # we check z>nodata to avoid starting a ridge in a hole
        # we start ridges between two consecutive thalwegs from pits and confluences
        # pits with one thalweg are excluded
#        pitlist = [i for i in self.nodedict \
#                   if ((self.nodedict[i]['type'] == df.PIT and len(self.nodedict[i]['thalweg'])>1)\
#                   or self.nodedict[i]['type'] == df.CONFLUENCE) and self.nodedict[i]['z']>self.terrain.nodata]
        pitlist = [i for i in self.nodedict if (self.nodedict[i]['type'] == df.PIT and i != -1) or
                    self.nodedict[i]['type'] == df.CONFLUENCE or self.nodedict[i]['type'] == df.BIFURCATION]  # ADDDDD BIFURCATION
       

        for ipt in pitlist:
            # if ipt!=testipt:
            #    continue
            pt = self.nodedict[ipt]

            # we get the list of thalwegs
            thalweglist = pt['thalweg']
            #print('Node ', ipt, thalweglist)

            (i, j) = pt['ij']
            # orderedneighbours is the neighbours of (i,j) in clockwise order
            orderedneighbours = [(i+di, j+dj)
                                 for (di, dj) in self.terrain.neighbour[i, j]]
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
                if t > 0:
                    pleft = self.thalwegdict[t]['polyline'][1]
                else:
                    pleft = self.thalwegdict[-t]['polyline'][-2]
                if tright > 0:
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
                    slopeneighbours = orderedneighbours[ileft:] + \
                        orderedneighbours[:iright+1]
                elif nt > 0:
                    slopeneighbours = orderedneighbours[ileft:iright+1]
                else:
                    slopeneighbours = orderedneighbours[ileft:] + \
                        orderedneighbours[:iright]
                slopeneighbours = [
                    v for v in slopeneighbours if self.terrain.isInside(v[0], v[1])]
                # the ridge in this hill must go through one of these pixels

                # take the steepest slope in pixels of slopeneighbours
                z = self.terrain.dtm[i, j]
                slope = [2*(self.terrain.dtm[v] - z)/(abs(v[0] - i) +
                                                      abs(v[1] - j) + 1) for v in slopeneighbours]
                maxslope = max(slope)
                maxcount = slope.count(maxslope)
                point = None
                if maxcount > 1:
                    maxcandidates = [slopeneighbours[i]
                                     for i in range(len(slope)) if slope[i] == maxslope]
                    maxcandidates.sort(
                        key=lambda x: df.order8.index((x[0]-i, x[1]-j)))
                    point = maxcandidates[-1]
                else:
                    maxpos = slope.index(maxslope)
                    point = slopeneighbours[maxpos]
                polylist = [(i, j), point]
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
                ileft = len(innerneighbours) - 1 - \
                    innerneighbours[::-1].index(pleft)
            else:
                ileft = innerneighbours.index(pleft)
            #print('tleft, tright', t, tright, ihill)
            #print('ileft, iright', ileft, iright)
            # slopeneighbours is the list of pixels between the thalwegs
            slopeneighbours = []
            if ileft >= iright:
                slopeneighbours = innerneighbours[ileft:iright-1:-1]
            else:
                slopeneighbours = innerneighbours[iright:] + \
                    innerneighbours[:ileft+1]
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
            (i, j) = point

            # get the first pixel outside the domain
            qleft = self.thalwegdict[-tleft]['polyline'][-1]
            qright = self.thalwegdict[-t]['polyline'][-1]
            ileft = outerneighbours.index(qleft)
            iright = outerneighbours.index(qright)
            slopeneighbours = []
            if ileft > iright:
                slopeneighbours = outerneighbours[ileft:iright-1:-1]
            else:
                slopeneighbours = outerneighbours[iright:] + \
                    outerneighbours[:ileft+1]
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

    def traceRidgesT(self, divide = False, correction = True):
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
        
        pixelclass = self.terrain.pixelclass
        pixelclass2 = self.terrain.pixelclass2
        

#        testir = [10363]
        for ihill in self.hilldict:
            floorkey={} # coordinate index recording ridge id
            junctionidx = {} # coordinate index to junctions
            ridgeindangling = {} # coordinate index for ridges on a dangling thalweg
            # for each coordinate, records the ridge id and the side of the thalweg to check
            # if another rige comes there, a junction should be inserted            junctionidx = {} # coordinate index to junctions
            ridgeindangling = {} # coordinate index for ridges on a dangling thalweg

            
            if ihill % 1000 == 0:
                print(ihill, end=" ")
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
                # print("id id thalweg :", dtoverlap)
                thalwegside = 0 # -1/0/1 if left/undefined/right side
                # on which side of the dangling line we come from
                if dtoverlap > 0:
                    if divide and 'startside' in self.ridgedict[ir]: 
                        thalwegside = self.ridgedict[ir]['startside']
                    else:
                        tmppoly = self.thalwegdict[dtoverlap]['polyline']
                        kt = tmppoly.index((i,j))
                        if kt == 0:
                            #print('kt = 0, thalweg', dtoverlap, 'ridge', ir)
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
                oldij = (i,j)
                while keepgoing:
                    if nextrounds:
                        ldr = self.terrain.neighbour[i,j]
                        v = self.terrain.getNeighbourHeights(i, j)
                        z = self.terrain.dtm[i,j]
    
                        #print(plgn)                    
                        [validridge, kmax, inlist] = self.computeAscent(i, j, v, z, ldr, plgn, danglingconfluence, confluencewedge, thalwegside, dtoverlap, danglingmap)

                        if not validridge:
                            print('validridge', ir, ihill)
                            exit(0)
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
                        if pixelclass[i,j] == df.CONFLUENCE or pixelclass2[i,j] == df.BIFURCATION :
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
                        if pixelclass[i,j] == df.CONFLUENCE or pixelclass2[i,j] == df.BIFURCATION :
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
                                # ajout 7/9/22
                                p4 = polylist[-2]
                                wedge2 = (p1, (i,j), p2)
                                wedge3 = (p1, (i,j), p3)
                                side23 = df.sideOfWedge(wedge2, p3)
                                side24 = df.sideOfWedge(wedge2, p4)
                                side32 = df.sideOfWedge(wedge3, p2)
                                side34 = df.sideOfWedge(wedge3, p4)
                                side2 = side23 * side24
                                side3 = side32 * side34
                                if side2 < 0:
                                    confluencewedge = wedge2
                                    thalwegside = side24
                                elif side3 < 0:
                                    confluencewedge = wedge3
                                    thalwegside = side34
                                elif side3 > 0:
                                    confluencewedge = wedge2
                                    thalwegside = side24
                                elif side2 > 0:
                                    confluencewedge = wedge3
                                    thalwegside = side34
                                # remplace
                                #thalwegside = df.sideOfWedge(confluencewedge, p1)
                                # fin ajout 7/9/22
                            # if we know on which side to stay
                            else:
                                jthalweg = ithalweg + thalwegside
                                if jthalweg == len(lthalweg):
                                    jthalweg = 0
                                kthalweg = ithalweg - thalwegside
                                if kthalweg == len(lthalweg):
                                    kthalweg = 0
                                it2 = df.getSecondPointPosition(lthalweg[jthalweg])
                                it3 = df.getSecondPointPosition(lthalweg[kthalweg])
                                p2 = self.thalwegdict[abs(lthalweg[jthalweg])]['polyline'][it2]
                                p3 = self.thalwegdict[abs(lthalweg[kthalweg])]['polyline'][it3]
                                confluencewedge = (p1,(i,j),p2)
                                side23 = df.sideOfWedge(confluencewedge, p3)
                                thalwegside = -side23
                                dtoverlap = lthalweg[jthalweg]
                            danglingconfluence = True
                        # if we are not at a confluence (regular thalweg pixel)
                        else:
                            # correction 8/9/22
                            if oldoverlap > 0 and dtoverlap != oldoverlap and pixelclass[i,j] != df.SADDLE :
                                tmppoly = self.thalwegdict[dtoverlap]['polyline']
                                kt = tmppoly.index((i,j))
                                if kt != 0:
                                    thalwegside = df.sideOfLine(tmppoly, kt, oldij)
                                else:
                                    thalwegside = 0
                            # fin correction
                            # if we don't know on which side to stay, find it out
                            if thalwegside == 0:
                                tmppoly = self.thalwegdict[dtoverlap]['polyline']
                                kt = tmppoly.index((i,j))
                                if kt != 0:
                                    #print("ça va planter", kt, i,j, ir, dtoverlap, oldij)
                                    thalwegside = df.sideOfLine(tmppoly, kt, oldij)
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
                    if pixelclass[i,j] == df.SLOPE:
                        # nothing special, keep going
                        keepgoing = True
                        pixelclass[i,j] = df.RIDGE
                        floorkey[i,j] = ir
                    # if divide is false (default) confluences are not considered 
                    # as end points since ridges only start from saddles
                    # if divide is true, ridges can start and end at confluences
                    # elif pixelclass[i,j] == df.THALWEG or (pixelclass[i,j] == df.CONFLUENCE and not divide) or (pixelclass[i,j] == df.BIFURCATION and not divide):
                    elif pixelclass[i,j] == df.THALWEG or (pixelclass[i,j] == df.CONFLUENCE and not divide) or (pixelclass2[i,j] == df.BIFURCATION and not divide):
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
                    elif (divide and pixelclass[i,j] == df.CONFLUENCE) or (divide and pixelclass2[i,j] == df.BIFURCATION) :
                    # elif divide and pixelclass[i,j] == df.CONFLUENCE :
                        if (i,j) in floorkey:
                            junction = False
                            keepgoing = False
                            if dtoverlap > 0:
                                if (i,j) in ridgeindangling:
                                    ridgeindangling[(i,j)].append((ir, thalwegside))
                                else:
                                    ridgeindangling[(i,j)] = [(ir, thalwegside)]
                        else:
                            floorkey[i,j] = ir
                            if dtoverlap > 0:
                                ridgeindangling[(i,j)] = [(ir, thalwegside)]
                            keepgoing = False
                            junction = False
                    elif pixelclass[i,j] == df.RIDGE:
                        keepgoing = False
                        junction = True
                    else: # ridge terminates at a saddle or peak
                        keepgoing = False
                        #outbound = False
                # end of while loop reached when keepgoing is false
            
            
                if validridge:
#                    idend = -1
#                    if pixelclass[i,j] == CONFLUENCE:
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
                        if pixelclass[i,j] == df.RIDGE:
                            pixelclass[i,j] = df.JUNCTION
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
                        #if pixelclass[i,j] == df.CONFLUENCE:
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
                        
    def orderRidgesAroundNodesC(self):
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
                (i, j) = pt['ij']
                neighbours = [(i+di, j+dj) for (di, dj) in df.ldr8]
                # arrange the lines in the same order as the neighbours
                try:
                    ridgelist.sort(key=lambda it: neighbours.index(
                        self.ridgedict[abs(it)]['polyline'][df.getSecondPointPosition(it)]))
                except (IndexError, ValueError) as err:
                    print('Exception raised in orderRidgesAroundNodes',
                          err, ipt, pt, neighbours)
                # check if two ridges overlap at the node
                for i in range(nr):
                    ir1 = ridgelist[i-1]
                    ir2 = ridgelist[i]
                    absir1 = abs(ir1)
                    absir2 = abs(ir2)
                    hill1 = self.ridgedict[absir1]['hill']
                    hill2 = self.ridgedict[absir2]['hill']
                    secondpoint = self.ridgedict[absir1]['polyline'][df.getSecondPointPosition(
                        ir1)]
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
                                print(
                                    'Exception raised in orderRidgesAroundNodes', err)
                                print('if hill1 == hill2', ipt, it, ir1, ir2)
                            if ts1 and ts2 :
                                if ts1 == ts2:
                                    print("Error at node", ipt, "ridges",
                                          ir1, "and", ir2, "on the same side")
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
                                    
    def insertIntoHills(self, i_in, j_in) :
        ###1###CORRECTION OF THALWEGS THAT HAVE BEEN CORRECTED
        key = self.nodeidx[(i_in, j_in)]
        if len(self.nodedict[key]['thalweg']) == 2 :
            if self.nodedict[key]['thalweg'][0] < 0 :
                thalwegsearch = - self.nodedict[key]['thalweg'][0]
                print("thalwegsearch :",thalwegsearch)
                nextthalwegsearch = self.nodedict[key]['thalweg'][1]
                print("nextthalwegsearch :",nextthalwegsearch)
            elif self.nodedict[key]['thalweg'][0] > 0 :
                thalwegsearch = - self.nodedict[key]['thalweg'][1]
                print("thalwegsearch :",thalwegsearch)
                nextthalwegsearch = self.nodedict[key]['thalweg'][0]
                print("nextthalwegsearch :",nextthalwegsearch)
            if self.thalwegdict[thalwegsearch]['righthill'] :
                righthill = self.thalwegdict[thalwegsearch]['righthill']
                print("righthill :", righthill)
                if thalwegsearch in self.hilldict[righthill]['boundary']:
                    boundary = self.hilldict[righthill]['boundary']
                    indice_boundary = boundary.index(thalwegsearch)
                    boundary.insert(indice_boundary + 1, nextthalwegsearch)
                    self.thalwegdict[nextthalwegsearch]['righthill'] = righthill
                if thalwegsearch in self.hilldict[righthill]['hole']:
                    boundary = self.hilldict[righthill]['hole']
                    indice_boundary = boundary.index(thalwegsearch)
                    boundary.insert(indice_boundary + 1, nextthalwegsearch)
                    self.thalwegdict[nextthalwegsearch]['righthill'] = righthill
            if self.thalwegdict[thalwegsearch]['lefthill'] :
                lefthill = self.thalwegdict[thalwegsearch]['lefthill']
                print("lefthill :", lefthill)
                if thalwegsearch in self.hilldict[lefthill]['boundary']:
                    boundary = self.hilldict[lefthill]['boundary']
                    indice_boundary = boundary.index(thalwegsearch)
                    boundary.insert(indice_boundary + 1, nextthalwegsearch)
                    self.thalwegdict[nextthalwegsearch]['lefthill'] = lefthill
                if thalwegsearch in self.hilldict[lefthill]['hole']:
                    boundary = self.hilldict[lefthill]['hole']
                    indice_boundary = boundary.index(thalwegsearch)
                    boundary.insert(indice_boundary + 1, nextthalwegsearch)
                    self.thalwegdict[nextthalwegsearch]['lefthill'] = lefthill
                
    def correctHills(self, terrain, fichier):
        B = self.computeBifurcation(terrain, fichier)
        b1 = B[0]
        b2 = B[1]
        L = self.newCulverts(b1, b2)
        culverts = L[0]
        self.orderThalwegsAroundNodes()
        for culvert in culverts :
            i_in = culvert[0][0]      
            j_in = culvert[0][1]
            i_out = culvert[-1][0] 
            j_out = culvert[-1][1]
            if self.terrain.pixelclass2[i_in,j_in] == df.BIFURCATION:
               self.insertIntoHills(i_in, j_in)

            if self.terrain.pixelclass2[i_out,j_out] == df.BIFURCATION:
               self.insertIntoHills(i_out, j_out)
               key = self.nodeidx[(i_out, j_out)]
               if len(self.nodedict[key]['thalweg']) == 1 :
                   nextthalwegsearch = self.nodedict[key]['thalweg'][0]
                   intersection = self.thalwegdict[nextthalwegsearch]['end']
    
    
    def InOutCulverts(self, ponceaux_avec_z, ponceaux_avec_z_prolonges) :
        """
        Updates the list of culverts by connecting them to the network's thalwegs

        Parameters
        ----------
        ponceaux_avec_z :
            list of pixelated culverts with height z
        ponceaux_avec_z_prolonges :
            list of extended pixelated culverts with height z

        Returns
        -------
        ponceaux_update :
            an updated list of culverts
        entrees_candidates :
            candidate culvert entrances
        sorties_candidates :
            candidate culvert exits
        min_locaux :
            list of local minimums
        max_locaux :
            list of local maximums

        """
        
        entrees_candidates = []   # list of candidate inlet at the culvert level
        ind_dict = []   # 
        sorties_candidates = []    # list of candidate outlets at culvert level
        ind_dict2 = [] # list representing the keys of the dictionary of candidate outputs
        ponceaux_update = []
        l = 1
    
        for ponceau in ponceaux_avec_z :
            ## Culvert index 
            ind = ponceaux_avec_z.index(ponceau)   # index of the culvert in the list of culverts_with_z
            # r = 0
            ## Determine the orientation of the elevation profile along the y-axis.
            if ponceau[0][1] < ponceau[-1][1]:
                note = "The elevation profile is oriented from left to right along the y-axis"
            elif ponceau[0][1] > ponceau[-1][1]:
                note = "The elevation profile is oriented from right to left along the y-axis"
            
            ## Create a list of z-elevations for each culvert
            values = [item[2] for item in ponceau]
            
            ###CANDIDATE INLETS
            ## Determine all local maximums and minimums at each culvert.
            min_locaux = []
            for i in range(1, len(values)-1):     
                if values[i-1] > values[i] and values[i+1] > values[i]:         
                    min_locaux.append(ponceau[i])
            max_locaux = []
            for i in range(1, len(values)-1):     
                if values[i-1] < values[i] and values[i+1] < values[i]:         
                    max_locaux.append(ponceau[i])
            
            ## Define rejected minima that are on the road: by detecting minima that are between very close 
            # elevation maxima (even the local minimum will have the same elevation)
            min_rejetes = []
            for min_point in min_locaux:
                for max_point1, max_point2 in zip(max_locaux[:-1], max_locaux[1:]):
                    # dist = df.calculate_distance(max_point1[0],max_point1[1],max_point2[0],max_point2[1])
                    z_1 = max_point1[2]
                    z_2 = max_point2[2]
                    # if dist < 5 :
                    if z_1 - z_2 < 0.5 :
                        if df.entre(max_point1, max_point2, min_point):
                           min_rejetes.append(min_point)
            
            ## Remove these rejected minimums from the list of local minimums
            min_locaux1 = [x for x in min_locaux if x not in min_rejetes]
            
            ## Separate the different possible cases
            min_locaux = []
            ## If there is only one local minimum in the list of local minima after rejecting those that are on the road
            if len(min_locaux1) == 1 :
                min_locaux = min_locaux1
            ## If there is no local minimum: in this case, the culvert must be extended by 3 or 5 metres (5 metres is acceptable given that for burning, they extend it by 5 metres)
            elif len(min_locaux1) == 0 :
                ## To simplify matters, the culverts have already been extended. All that needs to be done is to detect the culvert index and search for the extended one.
                ponceau = ponceaux_avec_z_prolonges[ind]
                # t = self.terrain.extendCulvert(ponceau, 5)
                # ponceau = t[1]
                # r = 1   ## This is just to show that the culvert has been extended in the elevation profile diagram
                # note2 = "The culvert was extended to find an outlet"
                # ind = ponceaux_avec_z_prolonges.index(ponceau)   ## assign the culvert index with the index in the list of extended culverts
                ## Determine local minimums again by eliminating those that are on the road (already detected)
                values = [item[2] for item in ponceau] ## Re-extract the z values of the culverts, given that we have an extended culvert
                min_locaux1 = []   # clear the list of local minimums
                ## Redefining local minimums with the new extended culvert
                for i in range(1, len(values)-1):     
                    if values[i-1] > values[i] and values[i+1] > values[i]:         
                        min_locaux1.append(ponceau[i]) 
                min_locaux1 = [x for x in min_locaux1 if x not in min_rejetes]  # reject local minimums that are the same because they are on the road
                if len(min_locaux1) == 1 :
                    min_locaux = min_locaux1
                elif len(min_locaux1) == 0 :
                    print("WARNING : Our method does not always work, and there is no candidate inlet even when extending!")
                elif len(min_locaux1) >= 2 :
                    values_z = [item[2] for item in min_locaux1] 
                    max_z = max(values_z)
                    for element in min_locaux1 :
                        if element[2] == max_z :
                            min_locaux.append(element)
            ## If there is more than one local minimum: choose the highest one
            elif len(min_locaux1) >= 2 :
                values_z = [item[2] for item in min_locaux1] 
                max_z = max(values_z)
                for element in min_locaux1 :
                    if element[2] == max_z :
                        min_locaux.append(element)
            ## List of candidate inlets at the culvert level
            entrees_candidates.append(min_locaux)
            ind_dict.append(ind+1)
            
            
            ###CANDIDATE OUTLETS
            liste_cand = []
            liste_cand_g = []
            cc = 0
            # cross the half-culvert to find the other local minimum to be considered as the outlet
            for k in range(int(len(ponceau)/2),len(ponceau)):   
                liste_cand.append(ponceau[k][2])   # store the z's in a list
                liste_cand_g.append(ponceau[k])   # store the culvert point in a list
            if min_locaux[0] in liste_cand_g :   ## if our candidate inlet belongs to the list of points belonging to this half-culvert
                # Then reset our lists and cross the culvert to the other side
                cc = 1  ## It's just to differentiate which half-culvert was used
                liste_cand = []  
                liste_cand_g = []
                for k in range(int(len(ponceau)/2)):
                    liste_cand.append(ponceau[k][2])
                    liste_cand_g.append(ponceau[k])
                    ## This part has made it possible to define the half-culvert end that does not contain our inlet
            # select the lowest points relative to the candidate inlet from among those determined in the list
            elements_inf = [x for x in liste_cand if x < min_locaux[0][2]]  
            sortie_cand =  []
            ## If the list of points below the candidate inlet does indeed contain elements and the candidate inlet is not on the same side of the culvert as these points
            if elements_inf and cc == 0 :
                sortie_z = elements_inf[0]   # select item 0 from the list of lower points to choose the point that is directly below the inlet
            # If the list of points below the candidate inlet contains elements and the candidate inlet is on the same side of the culvert as these points
            elif elements_inf and cc == 1 :
                sortie_z = elements_inf[-1]  # choose the last element of the inf elements since it is the one that directly corresponds to the lowest point of the inlet
            # Extend the culvert if no suitable outlets are found
            # If there is no lower element, the culvert must be extended
            elif not elements_inf :   
                # ind = ponceaux_avec_z.index(ponceau)
                ponceau = ponceaux_avec_z_prolonges[ind]
                # t = self.terrain.extendCulvert(ponceau, 5)
                # ponceau = t[1]
                # r = 1  ## This is just to show that the culvert has been extended in the elevation profile diagram
                # note2 = "The culvert was extended to find an outlet"
                # ind = ponceaux_avec_z_prolonges.index(ponceau)
                # print("liste_cand =", liste_cand)
                # print("liste_cand_g =", liste_cand_g)
                ##Redefine the candidate outlet in the same way with the extended culvert
                cc = 0
                for k in range(int(len(ponceau)/2),len(ponceau)):
                    liste_cand.append(ponceau[k][2])
                    liste_cand_g.append(ponceau[k])
                if min_locaux[0] in liste_cand_g :
                    cc = 1
                    liste_cand = []
                    liste_cand_g = []
                    for k in range(int(len(ponceau)/2)):
                        liste_cand.append(ponceau[k][2])
                        liste_cand_g.append(ponceau[k])
                elements_inf = [x for x in liste_cand if x < min_locaux[0][2]]
                # print("elements_inf = ", elements_inf)
                if elements_inf and cc == 0:
                    sortie_z = elements_inf[0]
                elif elements_inf and cc == 1 :
                    sortie_z = elements_inf[-1]
            ## Fill in the list of candidate outlets
            for kk in ponceau :
                if kk[2] == sortie_z:
                    sortie_cand.append(kk)
                    
            ## List of candidate outlets at the culvert level
            sorties_candidates.append(sortie_cand)
            ind_dict2.append(ind+1)
            
            ponceaux_update.append(ponceau)
            
            ## elevation profiles
            i_values = [entry[0] for entry in ponceau]
            j_values = [entry[1] for entry in ponceau]
            z_values = [entry[2] for entry in ponceau]
            
            i_min = [entry[0] for entry in min_locaux]
            j_min = [entry[1] for entry in min_locaux]
            z_min = [entry[2] for entry in min_locaux]
            
            i_max = [entry[0] for entry in max_locaux]
            j_max = [entry[1] for entry in max_locaux]
            z_max = [entry[2] for entry in max_locaux]
            
            i_initial = i_values[0]
            j_initial = j_values[0]
            
            distance = []
            distance_min = []
            distance_max = []
            for i, j in zip(i_values, j_values):
                dist =+ df.calculate_distance(i, j, i_initial, j_initial)
                distance.append(dist)
            for i, j in zip(i_min, j_min):
                dist =+ df.calculate_distance(i, j, i_initial, j_initial)
                distance_min.append(dist)
            for i, j in zip(i_max, j_max):
                dist =+ df.calculate_distance(i, j, i_initial, j_initial)
                distance_max.append(dist)
            i_min_s = [entry[0] for entry in sortie_cand]
            j_min_s = [entry[1] for entry in sortie_cand]
            z_min_s = [entry[2] for entry in sortie_cand]
            distance_min_s = []
            for i, j in zip(i_min_s, j_min_s):
                dist =+ df.calculate_distance(i, j, i_initial, j_initial)
                distance_min_s.append(dist)
            
        return ponceaux_update, entrees_candidates, sorties_candidates, min_locaux, max_locaux
                   
            
    def entrees_sorties_sansZ(self, entrees_candidates, sorties_candidates) :
        """
        Calculates candidate inputs and candidate outputs without height

        Parameters
        ----------
        entrees_candidates :
            list of candidate entry points
        sorties_candidates :
            list of candidate entry points

        Returns
        -------
        entrees_candidatess :
            list of candidate entry points without height z
        sorties_candidatess :
            list of candidate entry points without height z
        
        """ 
        
        entrees_candidatess = []
        for point in entrees_candidates :
            for x, y, _ in point :
                ponceau = [x, y]
                entrees_candidatess.append(ponceau)
    
        sorties_candidatess = []
        for point in sorties_candidates :
            for x, y, _ in point :
                ponceau = [x, y]
                sorties_candidatess.append(ponceau)
                
        return entrees_candidatess, sorties_candidatess
    
    def entrees_sorties_finaux(self, entrees_candidatess, sorties_candidatess, sorties_candidates) : 
        """
        Calculates candidate inputs and candidate outputs without height

        Parameters
        ----------
        entrees_candidatess :
            list of candidate entry points without height z
        sorties_candidatess :
            list of candidate entry points without height z
        sorties_candidates :
            list of candidate entry points

        Returns
        -------
        dictionnaire_entrees :
            dictionary of candidate inputs
        transfluences_entrees 
        dictionnaire_sorties :
            dictionary of candidate exits
        transfluences_sorties
        
        """ 
        
        dictionnaire_entrees = {}
        transfluences_entrees = []
          
        for point in entrees_candidatess :
            ind = entrees_candidatess.index(point)
            k = point[0]
            m = point[1]
            neighbors = [(k-1, m-1), (k, m-1), (k+1, m-1), (k-1, m), (k+1, m), (k-1, m+1), (k, m+1), (k+1, m+1)]      
            if self.terrain.pixelclass[k, m] == df.THALWEG or self.terrain.pixelclass[k, m] == df.SADDLE or self.terrain.pixelclass[k, m] == df.CONFLUENCE or self.terrain.pixelclass[k, m] == df.PIT :
                intersect = [k, m]
                transfluences_entrees.append(intersect)
                dictionnaire_entrees[ind] = intersect
            else :
                cand = []
                z_cand = []
                ind_s = entrees_candidatess.index(point)
                point_s = sorties_candidates[ind_s]
                z_s = point_s[0][2]
                for nk, nm in neighbors :
                    if self.terrain.pixelclass[nk, nm] == df.THALWEG or self.terrain.pixelclass[nk, nm] == df.SADDLE or self.terrain.pixelclass[nk, nm] == df.CONFLUENCE or self.terrain.pixelclass[nk, nm] == df.PIT :
                        nz = self.terrain.getHeight(nk,nm)
                        if nz > z_s :
                            z_cand.append(nz)
                            pt_cand = [nk, nm, nz]
                            cand.append(pt_cand)
                nz_min = min(z_cand)
                for pt in cand :
                    if pt[2] == nz_min :
                        intersect = [pt[0], pt[1]]
                        transfluences_entrees.append(intersect)
                        dictionnaire_entrees[ind] = intersect
                    else :
                        pass ## provide for a case where there is no inlet
    
        dictionnaire_sorties = {}                
        transfluences_sorties = []
          
        for point in sorties_candidatess :
            ind = sorties_candidatess.index(point)
            k = point[0]
            m = point[1]
            neighbors = [(k-1, m-1), (k, m-1), (k+1, m-1), (k-1, m), (k+1, m), (k-1, m+1), (k, m+1), (k+1, m+1)]      
            if self.terrain.pixelclass[k, m] == df.THALWEG or self.terrain.pixelclass[k, m] == df.SADDLE or self.terrain.pixelclass[k, m] == df.CONFLUENCE or self.terrain.pixelclass[k, m] == df.PIT :
                intersect = [k, m]
                transfluences_sorties.append(intersect)
                dictionnaire_sorties[ind] = intersect
            else :
                cand = []
                z_cand = []
                for nk, nm in neighbors :
                    if self.terrain.pixelclass[nk, nm] == df.THALWEG or self.terrain.pixelclass[nk, nm] == df.SADDLE or self.terrain.pixelclass[nk, nm] == df.CONFLUENCE or self.terrain.pixelclass[nk, nm] == df.PIT :
                        nz = self.terrain.getHeight(nk,nm)
                        z_cand.append(nz)
                        pt_cand = [nk, nm, nz]
                        cand.append(pt_cand)
                if z_cand :
                    nz_min = min(z_cand)
                    for pt in cand :
                        if pt[2] == nz_min :
                            intersect = [pt[0], pt[1]]
                            transfluences_sorties.append(intersect)
                            dictionnaire_sorties[ind] = intersect
                else :
                    intersect = [k, m]
                    transfluences_sorties.append(intersect) # if there is no outlet, stop at the outlet already detected
                    dictionnaire_sorties[ind] = intersect
    
        return dictionnaire_entrees, transfluences_entrees, dictionnaire_sorties, transfluences_sorties
    
    def newCulverts(self, transfluences_entrees, transfluences_sorties) :
        """
        Updates the list of culverts by connecting them to the network's thalwegs

        Parameters
        ----------
        transfluences_entrees :
            list of entry points
        transfluences_sorties :
            list of exit points

        Returns
        -------
        new_ponceaux_pixel :
            an updated list of culverts
        ponceaux_dictionnaire :
            an updated dictionary of culverts
        
        """
        
        new_ponceaux = []
        for i in range(len(transfluences_entrees)):
            ponceau = [transfluences_entrees[i], transfluences_sorties[i]]
            new_ponceaux.append(ponceau)
    
        new_ponceaux_pixel = []
        for ponceau in new_ponceaux :
            coords = ponceau
            # print("coords= ", coords)
            ponceaux = bresenham(coords)
            ponceaux_list = list(ponceaux)
            new_ponceaux_pixel.append(ponceaux_list)
        keys = []
        for i in range(len(new_ponceaux_pixel)):
            keys.append(i+1)
        pairs = zip(keys, new_ponceaux_pixel) 
        ponceaux_dictionnaire = dict(pairs)
    
        return new_ponceaux_pixel, ponceaux_dictionnaire


def readBDculverts(fichier, terrain) :
    """
    Reads culvert lines

    Parameters
    ----------
    fichier :
        path to culvert file
    terrain : 
        terrain

    Returns
    -------
    lignes_ponceaux : list
        a list storing culvert lines 
    dict_ponceaux : dictionary  
        a dictionary storing culvert lines

    """
    
    # Reading the BD file of the culverts already converted to GPKG
    data = gpd.read_file(fichier)
    
    culvert_crs = data.crs


    # Creation of a culvert dictionary
    ponceaux_dict = data.to_dict()  
            
    # Transformer for reprojecting data to the same coordinate system as the DTM
    transformer = Transformer.from_crs(culvert_crs, terrain.crs, always_xy=True)
    
    # Creating LineString culvert lines and storing them in a list
    lignes_ponceaux = []
    dict_ponceaux = {}
    j = 1
    for i, ligne in ponceaux_dict['geometry'].items():
        # print('ligne = ', ligne)
        line_origine = loads(str(ligne))
        coordonnees_ligne = list(line_origine.coords)
        # Creation of culvert lines while converting coordinates from one system to another
        line_cible = LineString(transformer.transform(x, y) for x, y in coordonnees_ligne)
        # line_cible = LineString((x, y) for x, y in coordonnees_ligne)
        # print(line_cible)
        lignes_ponceaux.append(line_cible)
        j += 1
    # print("lignes_ponceaux = ", lignes_ponceaux)
    # print(dict_ponceaux)
   
    return lignes_ponceaux, dict_ponceaux
            
