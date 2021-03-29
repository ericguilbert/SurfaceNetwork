# coding=utf-8
"""
terrain.py:
	Class storing a DTM from a raster grid, used for the computation of a surface network
	It includes an array storing the triangulation of the DTM
	
	author: Eric Guilbert
"""
from iomodule import readRasterDTM
from collections import deque

import numpy as np
import deffunction as df

class Terrain:
    # class describing a terrain defined by:
    # m, n: number of rows and columns of the terrain
    # x_size, y_size: raster resolution
    # upper_left_x, upper_left_y: coordinates of top left corner
    # dtm: an mxn array of elevations
    # pixelclass: an mxn array with a morphometric class for each pixel
    # neighbour: an mxn array with the list of neighbours of each pixel
    # diag: an m-1xn-1 array storing the orientation of each diagonal

    # constructor loading a raster from a file
    def __init__(self, chemin):
        """
        Constructor defining a new terrain from a raster file.
        The constructor loads the heights in an array called dtm
        and initialises other arrays: 
            pixelclass stores the morphometric class of each pixel,
            diag indicates the orientation of the diagonal, can be ±1.
            neighbour for each pixel contains the list of neighbouring pixels according to the triangulation
            boundary is a list of pixels located at the border outside the terrain
            inbound is the list of pixels located on the inner side of the border
            maxz and minz are scalars
        The constructor is not able yet to handle holes.

        Parameters
        ----------
        chemin : String
            path to the file containing the raster terrain.

        Returns
        -------
        None.

        """
        readRasterDTM(chemin, self) # reads the file and fill dtm array and box
        print(self.m, 'rows', self.n, 'columns')
        print('nodata =', self.nodata)
        outside = False
        # fill the pixelclass array
        # at the beginning all points are either slope or notin class
        self.pixelclass = np.ones((self.m, self.n))
        for i in range(self.m):
            for j in range(self.n):
                if self.dtm[i, j] > self.nodata:
                    self.pixelclass[i, j] = df.SLOPE
                else:
                    self.pixelclass[i, j] = df.NOTIN # outside the terrain
                    outside = True


        # if there are no pixel outside, the boundary is the border of the raster
        if not outside:
            self.boundary = [(i, -1) for i in range(-1, self.m)] + [(self.m, j) for j in range(-1, self.n)] \
                          + [(i, self.n) for i in range(self.m, -1, -1)] + [(-1, j) for j in range(self.n, -1, -1)]
            self.inbound = [(i, 0) for i in range(0, self.m)] + [(self.m - 1, j) for j in range(0, self.n)] \
                         + [(i, self.n - 1) for i in range(self.m - 1, 0, -1)] + [(0, j) for j in range(self.n - 1, 0, -1)]
        else:
            self.boundary, self.inbound = self.detectBoundary()

        # initialises other attributes
        self.diag = np.zeros((self.m - 1, self.n - 1))
        self.neighbour = np.zeros((self.m, self.n), dtype=np.ndarray)
        self.maxz = np.max(self.dtm)

    def getHeight(self, ij):
        """
        Return the height of the pixel (i,j). No domain check is done

        Parameters
        ----------
        ij : pair of integer
            coordinates of the pixel.

        Returns
        -------
        float
            height of the pixel.

        """
        return self.dtm[ij]

    def fromIndexToCoordinates(self, i, j):
        """
        Converts (i,j) matrix index into (x,y) coordinates. Coordinates are
        given at the centre of the pixel. While (i,j) are (row, column), (x,y) are
        (abscissa, ordinate)

        Parameters
        ----------
        i : integer
            row index.
        j : integer
            column index.

        Returns
        -------
        x : float
            x coordinate.
        y : float
            y coordinate.

        """
        x = j * self.x_size + self.upper_left_x + self.x_size / 2
        y = i * self.y_size + self.upper_left_y + self.y_size / 2
        return x, y

    # transform raster point coordinates into index
    def fromCoordinatesToIndex(self, x, y):
        """
        Converts (x,y) coordinates into (i,j) matrix index. Coordinates are
        given at the centre of the pixel. While (i,j) are (row, column), (x,y) are
        (abscissa, ordinate)

        Parameters
        ----------
        x : float
            x coordinate.
        y : float
            y coordinate.

        Returns
        -------
        i : integer
            row index.
        j : integer
            column index.

        """
        j = (x - self.upper_left_x) / self.x_size - 0.5
        i = (y - self.upper_left_y) / self.y_size - 0.5
        return int(i), int(j)

    # check if a point is inside or outside the terrain
    def isInside(self, i, j):
        """
        Checks if a pixel is inside the terrain or not. It can take into account
        holes if they have been marked, which is not done yet. Coordinates are
        row and column indexes and are supposed to be integer. If they are not,
        the result can be erroneous along the borders.

        Parameters
        ----------
        i : integer
            row index.
        j : integer
            column index.

        Returns
        -------
        boolean
            Returns True if the point is in the terrain.

        """
        return self.isInDomain(i, j) and self.pixelclass[i, j] > df.OUT

    def isInSubPixel(self, i, j):
        """
        Checks if a pixel is inside the terrain or not. It can take into account
        holes if they have been marked, which is not done yet. Coordinates are
        row and column indexes and can be given at subpixel accuracy (with float values). 
        This was mainly done for dealing with bilinear surfaces.

        Parameters
        ----------
        i : float
            row index.
        j : float
            column index.

        Returns
        -------
        boolean
            Returns True if the point is in the terrain.

        """
        return self.isInSubPixelDomain(i, j) and self.pixelclass[int(i), int(j)] > df.OUT

    def isInDomain(self, i, j):
        """
        Check if a pixel is inside the raster or not. This test is part of the
        check if a pixel is part of the terrain. Coordinates are integer values.
        If they are not, the result can be erroneous along the borders.

        Parameters
        ----------
        i : integer
            row index.
        j : integer
            column index.

        Returns
        -------
        boolean
            Returns True if the point is within the boundaries of the raster.

        """
        return (i >= 0) and (i < self.m) and (j >= 0) and (j < self.n)

    def isInSubPixelDomain(self, i, j):
        """
        Check if a pixel is inside the raster or not. This test is part of the
        check if a pixel is part of the terrain. Coordinates are float.
        This was mainly done for dealing with bilinear surfaces.

        Parameters
        ----------
        i : float
            row index.
        j : float
            column index.

        Returns
        -------
        boolean
            Returns True if the point is within the boundaries of the raster.

        """
        return (i >= 0) and (i <= self.m - 1) and (j >= 0) and (j <= self.n - 1)

    # build neighbour array
    def buildNeighbourArray(self):
        """
        Builds the list of neighbours of each vertex. Neighbours are the four
        vertices connected by the grid edges and the vertices connected by a 
        diagonal. Diagonals must have been computed before hand.

        Returns
        -------
        None.

        """
        m1 = self.m - 1
        n1 = self.n - 1
        dr = ((-1, -1), (-1, 0), (-1, 1), (0, 1), (1, 1), (1, 0), (1, -1), (0, -1))

        for i in range(self.m):
            for j in range(self.n):
                edge = [0, 1, 0, 1, 0, 1, 0, 1]
                if (i > 0) and (j > 0) and self.diag[i - 1, j - 1] == 1:
                    edge[0] = 1
                if (i > 0) and (j < n1) and self.diag[i - 1, j] == -1:
                    edge[2] = 1
                if (i < m1) and (j < n1) and self.diag[i, j] == 1:
                    edge[4] = 1
                if (i < m1) and (j > 0) and self.diag[i, j - 1] == -1:
                    edge[6] = 1
                neigh = []
                for ie, e in enumerate(edge):
                    if (e == 1):
                        neigh.append(dr[ie])
                self.neighbour[i, j] = neigh

    # define the boundary line surrounding the terrain
    # we need both an outer boundary and an inner boundary for sorting thalwegs
    def detectBoundary(self):
        """
        Compute the boundary of the terrain. Two boundaries are computed:
            the outer boundary where all pixels are outside the terrain but adjacent
            to it and the inner boundary where all pixels are inside but adjacent
            to the outside. The method should handle holes but this has not been
            done yet, which means that they are considered as low flat areas.

        Returns
        -------
        boundary : list of integer pairs
            list of pixels forming the outer boundary, some pixels can be outside the raster
        inbound : list of integer pairs
            list of pixels forming the inner boundary, all pixels are inside the terrain.

        """
        # computation of the outer boundary
        instart, jnstart = (0, 0)
        while self.pixelclass[instart, jnstart] == df.NOTIN:
            self.pixelclass[instart, jnstart] = df.OUT
            instart += 1
            if instart == self.m:
                instart = 0
                jnstart += 1
        # instart,jnstart is the first pixel inside the terrain
        istart = instart - 1
        jstart = jnstart # istart, jstart is on the outer boundary
        invldr = ((1, 0), (0, 1), (-1, 0), (0, -1))
        ldr = df.ldr4
        boundary = self.walkBoundary(istart, jstart, invldr)
        boundary.reverse()  # because we turn in the opposite direction around the global virtual pit

        # fill all outer pixels
        outstack = boundary[:]
        while outstack:
            (i, j) = outstack.pop()
            if self.isInDomain(i, j) and self.pixelclass[i, j] == df.NOTIN:
                self.pixelclass[i, j] = df.OUT
                pxlaround = [(i + di, j + dj) for (di, dj) in df.ldr4]
                outstack.extend(pxlaround)
        # all pixels in holes are considered inside
        # to be modified to handle holes
        for i in range(self.m):
            for j in range(self.n):
                if self.pixelclass[i, j] == df.NOTIN:
                    self.pixelclass[i, j] = df.SLOPE

        # computation of the inner boundary
        # istart,jstart is the first pixel inside the terrain
        istart = instart
        jstart = jnstart
        inbound = self.walkInbound(istart, jstart, ldr)
        inbound.reverse()

        return boundary, inbound

    def walkBoundary(self, istart, jstart, ldr):
        """
        This method starts with a point located on the outer boundary and
        walks along the border until coming back to the first point.

        Parameters
        ----------
        istart : integer
            row index of first point.
        jstart : integer
            column index of first point.
        ldr : list of pairs of ±1/0
            list of neighbour shifts, to be checked, ldr4 should work.

        Returns
        -------
        boundary : list of pairs of integer
            list of pixels located on the outer boundary.

        """
        g = deque(ldr)
        boundary = []
        i = istart
        j = jstart
        b = False
        #counter = 0
        # istart, jstart is outside the terrain
        while not b:
            boundary.append((i, j))
            pxlaround = [(i + di, j + dj) for (di, dj) in g]
            nexti = 4
            for ipxl, pxl in enumerate(pxlaround):
                if not self.isInside(pxl[0], pxl[1]):
                    nexti = ipxl
                    break
            if nexti == 4:
                break # there’s a problem?
            (i, j) = pxlaround[nexti]
            g.rotate(-nexti + 1)
            b = ((i, j) == (istart, jstart))
        return boundary

    def walkInbound(self, istart, jstart, ldr):
        """
        This method starts with a point located on the inner boundary and
        walks along the border until coming back to the first point.

        Parameters
        ----------
        istart : integer
            row index of first point.
        jstart : integer
            column index of first point.
        ldr : list of pairs of ±1/0
            list of neighbour shifts, to be checked, ldr4 should work.

        Returns
        -------
        inbound : list of pairs of integer
            list of pixels located on the inner boundary.

        """
        inbound = []
        g = deque(ldr)
        i = istart
        j = jstart
        b = False
        # istart, jstart is inside the terrain
        while not b:
            inbound.append((i, j))
            pxlaround = [(i + di, j + dj) for (di, dj) in g]
            #        print('pxl', i, j, self.pixelclass[i,j], 'neighb', pxlaround)
            nexti = 4
            for ipxl, pxl in enumerate(pxlaround):
                if self.isInside(pxl[0], pxl[1]):
                    nexti = ipxl
                    break
            if nexti == 4:
                break
            (i, j) = pxlaround[nexti]
            g.rotate(-nexti + 1)
            b = ((i, j) == (istart, jstart))
        return inbound

    def diagonalisation(self):
        """
        Define a diagonal in each cell that connects to the lowest point of the cell.
        This was used before we maximise the number of saddles.

        Returns
        -------
        None.

        """
        # we triangulate the domain by cutting each grid cell in two
        # the diagonal connects the lowest point to limit the number of pits
        # and facilitate the flow

        # the matrix diag stores the orientation of each diagonal in a cell,
        # -1 diagonal is SW-NE
        # +1 diagonal is NW-SE
        # 0  no diagonal (which is possible but would be a bad idea unless under constraint)
        # diag = ones((m-1,n-1)) # nw-se orientation
        # diag = -1*ones((m-1,n-1)) # sw-ne orientation
        # orientation towards the lowest point
        for i in range(self.m - 1):
            for j in range(self.n - 1):
                min1 = min(self.dtm[i, j], self.dtm[i + 1, j + 1])
                min2 = min(self.dtm[i + 1, j], self.dtm[i, j + 1])
                if (min1 < min2):
                    self.diag[i, j] = 1
                else:
                    self.diag[i, j] = -1

        # neighbour is a matrix where each cell contains the list of neighbours of a point
        self.buildNeighbourArray()

    # raster triangulation with consideration of saddles and thalwegs
    def saddleAndThalwegAwareDiagonalisation(self):
        """
        Flip the diagonals towards the highest point of the cell. If saddles and
        thalwegs have been computed before (which should be the case), diagonals
        connecting to saddles and followed by thalwegs are not flipped.

        Returns
        -------
        None.

        """
        # we modify the existing triangulation to compute ridges
        # the diagonal connects the highest point in each cell
        # but leaves thalwegs and saddles unmodified

        # the matrix diag stores the orientation of each diagonal in a cell,
        # -1 diagonal is SW-NE between (i+1,j) and (i,j+1)
        # +1 diagonal is NW-SE between (i,j) and (i+1,j+1)
        # 0  no diagonal (which is possible but would be a bad idea unless under constraint)

        # set the diagonals around passes according to ridges and thalwegs

        # go through all diagonals
        for i in range(self.m - 1):
            for j in range(self.n - 1):
                # check if one point is a pass or two points are thalwegs or confluences
                # in that case, we don't modify the diagonal
                #                if (self.diag[i,j] == 1 and \
                #                ( (self.pixelclass[i,j] != SADDLE and self.pixelclass[i+1,j+1] != SADDLE) and (self.pixelclass[i,j] == SLOPE or self.pixelclass[i+1,j+1] == SLOPE)))\
                #                or (self.diag[i,j] == -1 and \
                #                ( (self.pixelclass[i+1,j] != SADDLE and self.pixelclass[i,j+1] != SADDLE) and (self.pixelclass[i+1,j] == SLOPE or self.pixelclass[i,j+1] == SLOPE))):
                if self.pixelclass[i, j] == df.SADDLE or self.pixelclass[i, j + 1] == df.SADDLE or self.pixelclass[
                    i + 1, j] == df.SADDLE or self.pixelclass[i + 1, j + 1] == df.SADDLE:
                    continue
                if (self.pixelclass[i, j] != df.SLOPE and self.pixelclass[i + 1, j + 1] != df.SLOPE) \
                        or (self.pixelclass[i + 1, j] != df.SLOPE and self.pixelclass[i, j + 1] != df.SLOPE):
                    continue
                max1 = max(self.dtm[i, j], self.dtm[i + 1, j + 1])
                max2 = max(self.dtm[i + 1, j], self.dtm[i, j + 1])
                if (max1 > max2):
                    self.diag[i, j] = 1
                else:
                    self.diag[i, j] = -1

        # neighbour is a matrix where each cell contains the list of neighbours of a point
        self.buildNeighbourArray()

    # raster triangulation with consideration of saddles
    def saddleAndFlowAwareDiagonalisation(self, passdict):
        """
        Orientation of the diagonals towards the lowest point of each cell. Saddles
        that have been computed before are preserved by aligning the diagonals 
        with their critical lines.

        Parameters
        ----------
        passdict : dict
            Dictionary containing the saddles. Critical lines are defined (their
            first segment) for each saddle.

        Returns
        -------
        None.

        """
        # we triangulate the domain by cutting each grid cell in two
        # first, diagonals are set to preserve critical lines around saddles
        # second, diagonals connect the lowest point to limit the number of pits
        # and facilitate the flow
        # this is good for extracting thalwegs and fit them on a drainage network
        # but it leads to a very large number of peaks, 

        # the matrix diag stores the orientation of each diagonal in a cell,
        # -1 diagonal is SW-NE
        # +1 diagonal is NW-SE
        # 0  no diagonal (which is possible but would be a bad idea unless under constraint)
        # diag = ones((m-1,n-1)) # nw-se orientation
        # diag = -1*ones((m-1,n-1)) # sw-ne orientation

        # set the diagonals around passes according to ridges and thalwegs
        for ipss, pss in passdict.items():
            linelist = pss['line']
            for l in linelist:
                p1 = l[0]
                p2 = l[1]
                # if the line is diagonal
                if p1[0] != p2[0] and p1[1] != p2[1]:
                    # check which diagonal it is
                    if p2[0] < p1[0]:
                        if p2[1] < p1[1]:
                            self.diag[p2] = 1
                        else:
                            self.diag[p2[0], p1[1]] = -1
                    else:
                        if p2[1] < p1[1]:
                            self.diag[p1[0], p2[1]] = -1
                        else:
                            self.diag[p1] = 1

        # orientation towards the lowest point
        for i in range(self.m - 1):
            for j in range(self.n - 1):
                if self.diag[i, j] == 0:
                    min1 = min(self.dtm[i, j], self.dtm[i + 1, j + 1])
                    min2 = min(self.dtm[i + 1, j], self.dtm[i, j + 1])
                    if (min1 < min2):
                        self.diag[i, j] = 1
                    else:
                        self.diag[i, j] = -1

        # neighbour is a matrix where each cell contains the list of neighbours of a point
        self.buildNeighbourArray()

    # get all 8 neighbour points according to ldr
    def getAllNeighbourHeights(self, i, j, ldr=df.ldr8):
        """
        Returns the height of all the neighbours of a pixel. Neighbouring pixels
        are indicated by ldr. By default all eight neighbours are returned.

        Parameters
        ----------
        i : integer
            DESCRIPTION.
        j : integer
            DESCRIPTION.
        ldr : list of pairs of integer, optional
            Contains the vectors to each neighbour. The default is df.ldr8.

        Returns
        -------
        v : list of floats
            Heights of neighbouring points.

        """
        # ldr = ((-1,-1), (-1,0), (-1,1), (0,1), (1,1), (1,0), (1,-1), (0,-1))
        l = len(ldr)
        v = [0] * l
        for k in range(l):
            dr = ldr[k]
            ik = i + dr[0]
            jk = j + dr[1]
            if self.isInside(ik, jk):
                #            if ((ik>=0) and (ik<self.m) and (jk>=0) and (jk<self.n)):
                v[k] = self.dtm[ik, jk]
            else:
                v[k] = self.nodata
        return v

    # get neighbour points
    def getNeighbourHeights(self, i, j):
        """
        Returns the height of all the neighbours of a pixel. Neighbouring pixels
        are given by the neighbour array and depend on the triangulation.

        Parameters
        ----------
        i : integer
            DESCRIPTION.
        j : integer
            DESCRIPTION.

        Returns
        -------
        v : list of floats
            Heights of neighbouring points.

        """
        ldr = self.neighbour[i, j]
        l = len(ldr)
        v = [0] * l
        for k in range(l):
            dr = ldr[k]
            ik = i + dr[0]
            jk = j + dr[1]
            if self.isInside(ik, jk):
                #            if ((ik>=0) and (ik<self.m) and (jk>=0) and (jk<self.n)):
                v[k] = self.dtm[ik, jk]
            else:
                v[k] = self.nodata
        return v

    def checkSaddle(self, i, j, ldr):
        """
        Checks if a pixel is a saddle by comparing it with its neighbours given in ldr.
        The method can be sped up. The method return an object containing the list of
        critical lines (their first segment only)

        Parameters
        ----------
        i : integer
            row index.
        j : integer
            column index.
        ldr : list of pair of integers
            List of vectors to neighbouring points.

        Returns
        -------
        criticalline : list of numbers
            The list contains the first segment of each critical line. Each line is
            defined by the coordinates of the first point (i,j), the coordinates of 
            the second point, ±1 for the type of line and the height of the second point.

        """
        z = self.dtm[i, j]
        nv = len(ldr)
        v = self.getAllNeighbourHeights(i, j, ldr)  # elevation of 8 neighbours
        signe = [0] * nv  # sign of elevation difference with centre
        Nc = 0
        criticalline = []
        #       print(v, nv)

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
            # to be optimised: replace the del by a new list
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
            # a critical line is defined by the central pixel, the second pixel, its type and elevation
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
        return criticalline

    # computation of saddles directly from the raster
    def computeSaddles(self):
        """
        Detection of all saddles on the raster. Comparison is done on all eight neighbours.
        The objective is to detect as many saddles as possible, they will be sorted after.
        In the pixelclass, all saddles are marked TOBEPROCESSED for that.

        Returns
        -------
        saddledict : dictionary
            Contains all the saddles with their critical lines.
        saddleidx : dictionary
            Hash table storing saddle coordinates for indexing.

        """
        # matrix dimensions
        A = self.dtm
        pixelclass = self.pixelclass
        m = self.m
        n = self.n

        # these are counters for counting saddle points
        count_p = 0  # number of saddles (counting multiplicity)
        count_n = 0  # number of saddles (without multiplicity)

        saddledict = {}  # dictionary of saddles
        saddleidx = {}  # dictionary indexing saddle coordinates

        # 8 directions around a pixel, goes clockwise
        ldr = ((-1, -1), (-1, 0), (-1, 1), (0, 1), (1, 1), (1, 0), (1, -1), (0, -1))
        nv = len(ldr)

        # here starts the computation of saddle points
        # we go through each point and check its elevation against its neighbours
        for i in range(m):
            for j in range(n):
                if pixelclass[i, j] > df.OUT:
                    z = A[i, j]
                    v = self.getAllNeighbourHeights(i, j)  # elevation of 8 neighbours
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
                        pixelclass[i, j] = df.TOBEPROCESSED
                        saddleidx[i, j] = count_n
                        count_n += 1
                        count_p += (Nc / 2 - 1)
        return saddledict, saddleidx
