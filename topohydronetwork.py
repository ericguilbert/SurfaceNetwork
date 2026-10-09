# -*- coding: utf-8 -*-
"""
Created on Tue Oct  1 23:44:12 2023

@author: Yassmine Zada
"""

from time import perf_counter_ns
from math import sqrt

from extendedsurfacenetwork import ExtendedSurfaceNetwork
from celltype import CellType
import deffunction as df
from terrain import bresenham

import geopandas as gpd
from shapely.geometry import LineString
from shapely.wkt import loads
from pyproj import Transformer


class TopoHydroNetwork(ExtendedSurfaceNetwork):

        
    def buildTopoHydroFromTerrain(self, terrain, culvert, sd = True, roadheight = 0.5):
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
        self.terrain = terrain
#        self.culvert = culvert
        self.roadheight = roadheight
        
        self.buildFlowThalwegsFromTerrain(terrain, correction = False)
        start = perf_counter_ns()
        inlet, outlet = self.computeCulvertInletsAndOutlets(culvert)
        self.createCulverts(inlet, outlet)
        end = perf_counter_ns()
        print("Culvert computation time:", end - start)
        start = perf_counter_ns()
        self.correctThalwegs()
        end = perf_counter_ns()
        print("Thalweg split time:", end - start)
        start = perf_counter_ns()
        self.addNewWatercourses(correction=True)
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
        return
#        self.reassignmentnode()
    
    def isDownslopeLine(self, cellvalue: "CellType") -> bool:
        """
        returns True if the line is a surface watercourse.
        The method overrides the method that returns True for thalwegs

        Parameters
        ----------
        cellvalue : "CellType"
            A celltype value.

        Returns
        -------
        bool
            True if the cell value is a surface watercourse.

        """
        return cellvalue & CellType.SURFACE_FLAG and cellvalue & CellType.LINE

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
        
        
    def computeCulvertInletsAndOutlets(self, culvert) :
        """
        Calculates connection points between thalwegs and culverts 

        Parameters
        ----------
        culvert :
            list of pixelated culverts. Each culvert contains two pixel lines,
            the pixelated original culvert and the pixelated extended culvert
            
        Returns
        -------
        node:

        """
        
        ponceaux_avec_z = culvert[0]
        ponceaux_avec_z_prolonges = culvert[1]

        # position inlets and outlets along culverts        
        entrees_candidates, sorties_candidates = self.InOutCulverts(ponceaux_avec_z, ponceaux_avec_z_prolonges)  ##Calculation of culvert inlets and outlets
        
        # remove the z coordinates to inlet and outlet nodes
        T1 = self.entrees_sorties_sansZ(entrees_candidates, sorties_candidates)
        entrees_candidatess = T1[0]
        sorties_candidatess = T1[1]

        # relocate inlets and outlets on thalwegs and watercourse nodes
        T2 = self.entrees_sorties_finaux(entrees_candidatess, sorties_candidatess, sorties_candidates)
        inlet = T2[0]
        outlet = T2[1]
        
        # node =  []
        # for i in range(len(inlet)) :
        #     node.append(inlet[i])
        #     node.append(outlet[i])
        
        return inlet, outlet

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
        
        lignes_ponceaux = self.culvert
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

    def createCulverts(self, inlet, outlet):
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
        
        
        # self.nodekey contains the index of the next node
        # self.linekey contains the index of the next line

        # compute the line of each culvert
        culvertline = self.newCulverts(inlet, outlet)
        # avant, L[0] liste des ponceaux, L[1] dictionnaire des ponceaux
        # maintenant, on ne retourne que la liste
        
        for iline, line in enumerate(culvertline):
            point = inlet[iline] # take the starting point of the culvert
            i = point[0]
            j = point[1]
            z = A[i, j]
        
            startkey = self.nodekey
            if (i,j) in self.nodeidx:
                # modifier le nœud existant
                startkey = self.nodeidx[(i,j)]   # retrieve the index of this node
                self.nodedict[startkey]['line'].append(self.linekey)   # add the culvert to the line list
                self.nodedict[startkey]['type'] = self.nodedict[startkey]['type'] | CellType.CULVERTNODE
            else :
                # créer un nouveau nœud
                self.nodeidx[(i,j)] = startkey  # add to nodeidx the following indices at the bifurcations

                # add the culvert to the nodedict
                self.nodedict[startkey] = {
                    'ij': (i, j),
                    'z': z,
                    'line': [self.linekey], # the only line at this node is the culvert.
                    'type': CellType.CULVERTNODE
                }
                self.nodekey += 1
            pixelclass[i,j] = pixelclass[i,j] | CellType.CULVERTNODE_FLAG
                # pixelclass[i, j] = df.BIFURCATION   #The reason why we don't add it in pixelclass is to preserve the fact that it is a thalweg dividing in two

            point = outlet[iline] # take the end point of the culvert
            i = point[0]
            j = point[1]
            z = A[i, j]
            # créer le ponceau à ajouter dans linedict
        
            key = self.nodekey    
            if (i,j) in self.nodeidx:
                # modifier le nœud existant
                key = self.nodeidx[(i,j)]   # retrieve the index of this node
                self.nodedict[key]['line'].append(-self.linekey) # negative because it is an end node
                self.nodedict[key]['type'] = self.nodedict[key]['type'] | CellType.CULVERTNODE
            else :
                # créer un nouveau nœud
                self.nodeidx[(i,j)] = self.nodekey  # add to nodeidx the following indices at the bifurcations

                # add the culvert to the nodedict
                self.nodedict[self.nodekey] = {
                    'ij': (i, j),
                    'z': z,
                    'line': [-self.linekey], # negative because it is the end node
                    'type': CellType.CULVERTNODE
                }
                self.nodekey += 1
            pixelclass[i,j] = pixelclass[i,j] | CellType.CULVERTNODE_FLAG
            
            self.linedict[self.linekey] = {
                    'start': startkey,
                    'polyline': line,
                    'end': key,
                    'type': CellType.CULVERT
                }
            self.linekey += 1

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
        
        
        nodedict = self.nodedict
        culvertnode = [i for i in nodedict if nodedict[i]['type'] & CellType.CULVERTNODE_FLAG]
        for k in culvertnode:
            
            # retrieve the bifurcations indices
            node = nodedict[k]
            i = node['ij'][0]
            j = node['ij'][1]
            
            if CellType.isThalweg(self.terrain.pixelclass[i,j]):
                
                idthalweg = self.floorkey[i,j]
                node['line'].append(-idthalweg)
                node['line'].append(self.linekey)
                node['type'] = node['type'] | CellType.TJUNCTION

                # split the existing thalweg in two
                t = self.linedict[idthalweg]
                ptlist = t['polyline']
                iend = t['end']
                icfl = ptlist.index((i,j))
                ptlist1 = ptlist[:icfl+1]
                ptlist2 = ptlist[icfl:]
                # modify the first thalweg to stop at the confluence
                t['polyline'] = ptlist1
                t['end'] = k   # the culvert node is the end of the thalweg
                # create a new thalweg
                self.linedict[self.linekey] = {
                    'start' : k,
                    'polyline' : ptlist2,
                    'end' : iend,
                    'lefthill' : [],
                    'righthill' : [],
                    'type' : CellType.THALWEG
                }
                
                # we update the node at the end of the new thalweg
                ndt = self.nodedict[iend]['line']
                ndt.remove(-idthalweg)
                ndt.append(-self.linekey)
                # we update floorkey with nextthalweg
                for k in ptlist2[1:-1]:
                    self.floorkey[k] = self.linekey
                self.linekey += 1
                self.terrain.pixelclass[i,j] = node['type']

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
        
    
    def addNewWatercourses(self, correction=True):
        
        culvert = [k for k in self.linedict if self.linedict[k]['type'] & CellType.CULVERT_FLAG]
        
        endnode = []
        for ic in culvert:
            # culvert start node is on a thalweg, no need to add a watercourse
            # retrieve the end of the culvert
            inode = self.linedict[ic]['end']
            if self.nodedict[inode]['line'] == [-ic]:
                endnode.append(inode)

        print("list of isolated culvert nodes = ", endnode)
        
        for ic in culvert:
            inode = self.linedict[ic]['end']
            pt = self.nodedict[inode]
            if not (pt['line'] == [-ic]):
                continue
            (i,j) = pt['ij']
            polylist = [(i,j)]
            self.traceLineDown(inode, polylist, correction, CellType.SURFACE)
            pt['type'] |= CellType.TJUNCTION_FLAG

                                    
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
        Compute inlets and outlets of culverts and return the list of culverts
        with candidate inlets and outlets and min and max along culverts

        Parameters
        ----------
        ponceaux_avec_z :
            list of pixelated culverts with height z
        ponceaux_avec_z_prolonges :
            list of extended pixelated culverts with height z

        Returns
        -------
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
        l = 1
    
        for ind, ponceau in enumerate(ponceaux_avec_z):
            ## Culvert index 
#            ind = ponceaux_avec_z.index(ponceau)   # index of the culvert in the list of culverts_with_z
            # r = 0
            ## Determine the orientation of the elevation profile along the y-axis.
            if ponceau[0][1] < ponceau[-1][1]:
                #note = "Le profil d'élévation est orienté de gauche à droite selon l'axe y"
                note = "The elevation profile is oriented from left to right along the y-axis"
            elif ponceau[0][1] > ponceau[-1][1]:
                #note = "Le profil d'élévation est orienté de droite à gauche selon l'axe y"
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
                    if z_1 - z_2 < self.roadheight :
                        if df.between(max_point1, max_point2, min_point):
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
            # plt.figure(figsize=(8, 8))
            # plt.plot(distance , z_values, marker='o', label='Pixels du ponceau')
            # plt.scatter(distance_min , z_min, marker='o', s=200, facecolors = 'none', edgecolors='red', linewidths=2, label='Entrées Candidates')
            # plt.scatter(distance_max , z_max, marker='o', s=200, facecolors = 'none', edgecolors='orange', linewidths=2, label='Maximums (Route)')
            # plt.scatter(distance_min_s , z_min_s, marker='o', s=200, facecolors = 'none', edgecolors='green', linewidths=2, label='Sorties Candidates')
            # plt.xlabel('Distance along the culvert')
            # plt.ylabel('Elevation Z')
            # plt.title(f"Elevation profile along the culvert {l}")
            # l = l+1
            # plt.grid(True)
            # plt.legend()
            # plt.show()
            
        return entrees_candidates, sorties_candidates
                   
            
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
        transfluences_entrees 
        transfluences_sorties
        
        """ 
        
        dictionnaire_entrees = {}
        transfluences_entrees = []
          
        pixelclass = self.terrain.pixelclass
        for point in entrees_candidatess :
            ind = entrees_candidatess.index(point)
            k = point[0]
            m = point[1]
            neighbors = [(k-1, m-1), (k, m-1), (k+1, m-1), (k-1, m), (k+1, m), (k-1, m+1), (k, m+1), (k+1, m+1)]      
            if CellType.isThalweg(pixelclass[k, m]) or CellType.isCourseNode(pixelclass[k, m]) :
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
                    if CellType.isThalweg(pixelclass[nk, nm]) or CellType.isCourseNode(pixelclass[nk, nm]) :
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
            if CellType.isThalweg(pixelclass[k, m]) or CellType.isCourseNode(pixelclass[k, m]) :
                intersect = [k, m]
                transfluences_sorties.append(intersect)
                dictionnaire_sorties[ind] = intersect
            else :
                cand = []
                z_cand = []
                for nk, nm in neighbors :
                    if CellType.isThalweg(pixelclass[nk, nm]) or CellType.isCourseNode(pixelclass[nk, nm]) :
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
    
        return transfluences_entrees, transfluences_sorties
    
    def newCulverts(self, inlet, outlet) :
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
        for i in range(len(inlet)):
            ponceau = [inlet[i], outlet[i]]
            new_ponceaux.append(ponceau)
    
        new_ponceaux_pixel = []
        for ponceau in new_ponceaux :
            coords = ponceau
            # print("coords= ", coords)
            ponceaux = bresenham(coords)
#            ponceaux_list = list(ponceaux)
#            new_ponceaux_pixel.append(ponceaux_list)
            new_ponceaux_pixel.append(ponceaux)
    
        return new_ponceaux_pixel


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
    print(culvert_crs)
    print(terrain.crs)

    # Creation of a culvert dictionary
    ponceaux_dict = data.to_dict()  
            
    # Transformer for reprojecting data to the same coordinate system as the DTM
    transformer = Transformer.from_crs(culvert_crs, terrain.crs, always_xy=True)
    
    # Creating LineString culvert lines and storing them in a list
    lignes_ponceaux = []
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
   
    return lignes_ponceaux
            
