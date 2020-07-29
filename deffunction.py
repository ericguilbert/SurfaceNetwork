# coding=utf-8
"""
deffunction.py:
    Several functions used to compute different characteristics from the terrain
    or the surface network
    
    author: Eric Guilbert
"""

import math
from matplotlib.pylab import * # check if used

from collections import deque

# pixel codes
TOBEPROCESSED = 999
NOTIN = -999
OUT = -99
PIT = -3
CONFLUENCE = -2
THALWEG = -1
SADDLE = 0
RIDGE = 1
RIDWEG = -5  # overlapping ridge and thalweg, not used
JUNCTION = 2
PEAK = 3
SLOPE = 99

ldr8 = ((-1, -1), (-1, 0), (-1, 1), (0, 1), (1, 1), (1, 0), (1, -1), (0, -1))
ldr4 = ((-1, 0), (0, 1), (1, 0), (0, -1))
order4 = ((-1, 0), (0, -1), (0, 1), (1, 0)) # lexicographic order of neighbours

# step for gradient descent on bilinear patches
step = 0.1

# epsilon defining numerical zero
epsilon = 1.e-02

def gradLength(z, v, ldr):
    """
    gradLength returns rounded gradient norms from a point to a series of neighbouring points
    distance is 3/2 along the diagonal because we want the closest direction to gradient.
    For one triangle
    The steepest descent is a vector going from one vertex to the opposite edge. The position
    on the opposite edge has an abscissa in [0;1], 0 if the vector goes along the pixel edge p1,
    1 if it goes along the diagonal to vertex p2. If it is below ½, the vector is closer to p1,
    if it is above, it is closer to p2. Since grad z = (z1-z0, z2-z1) if p0 is the local origin,
    the gradient will be closer to p2 if (z2-z0)>3/2(z2-z1)

    Parameters
    ----------
    z : float
        Height of the central point.
    v : float list
        Elevation of neighbouring points.
    ldr : pair list
        Direction of neighbouring points. Should be ldr8 for eight directions

    Returns
    -------
    dz : float list
        Contains the gradient norm (the slope) in each direction.
    """
    dz = []
    for i, p in enumerate(v):
        dr = ldr[i]
        if dr[0] == 0 or dr[1] == 0:
            dz.append(p - z) # the norm is simply the height difference
        else:
            dz.append((p - z) * 2 / 3) # 2/3 of the height difference
    return dz


def lexicoSign(di, dj):
    """
    Computes the lexicographic sign based on the coordinate difference.
    Coordinates correspond to raster coordinates: (1,1) in top left corner
    For two points p1 and p2, we note (di, dj) = p2 - p1.
    If p1 is closer to the raster top, the result is positive.
    If p1 is closer to the raster left, the result is positive.

    Parameters
    ----------
    di : integer
        coordinate difference in i.
    dj : integer
        coordinate difference in j.

    Returns
    -------
    integer
        Either 1 (positive) or -1 (negative).

    """
    if di > 0:
        return 1
    if di < 0:
        return -1
    if dj > 0:
        return 1
    return -1


# function computing if two segments cross at their middle without overlapping
def crossingSegment(segment1, segment2):
    """
    The function checks if the two segments cross in a single poirt at their middle.
    It is used to check if critical lines starting at two neighbouring saddles are
    in conflict or not. Since segment coordinates are integer values, two segments
    can only cross, be equal (but this case cannot occur for two different saddles)
    or be disjoint (and parallel).

    Parameters
    ----------
    segment1 : a pair of pairs
        contains the coordinates of its two points.
    segment2 : a pair of pairs
        contains the coordinates of ist two points.

    Returns
    -------
    boolean
        True if the segments cross at their middle.
        False otherwise

    """
    p1 = segment1[0]
    p2 = segment1[1]
    q1 = segment2[0]
    q2 = segment2[1]
    # first check if they have the same middle point
    if (p1[0] + p2[0] == q1[0] + q2[0]) and (p1[1] + p2[1] == q1[1] + q2[1]):
        # the cross product is non zero if they cross, zero if they overlap
        return bool((p2[0] - p1[0]) * (q2[1] - q1[1]) -
                    (q2[0] - q1[0]) * (p2[1] - p1[1]))
    return False


# function returning the angle of a vector from vector (1,0)
def anglex(vctr):
    """
    Function returning the angle between a vector and the vector (1,0).

    Parameters
    ----------
    vctr : a pair of integer values
        vctr is a vector defined by its (i, j) coordinates.

    Returns
    -------
    float
        an angle between -pi and +pi.

    """
    if vctr == (0, 0):
        return 0
    lenvector = math.hypot(vctr[0], vctr[1])
    nvctr = (vctr[0] / lenvector, vctr[1] / lenvector)
    angl = math.atan2(nvctr[1], nvctr[0])  # angle is in ]-pi, pi]
    return angl


# function returning the anglex of the segment of a polyline starting at orgn
def polylineangle(polyline, orgn):
    """
    Returns the angle of the first segment of a polyline with the i-axis.
    The first segment depends on the orientation considered, which is indicated
    by passing the point to take as the starting point.

    Parameters
    ----------
    polyline : list of pairs
        polyline defined by a list of (i,j) coordinates.
    orgn : a pair of integer
        the coordinate of the first point (either the first or the last of the polyline).

    Returns
    -------
    float
        an angle between -pi and +pi.

    """
    if orgn == polyline[0]:
        secondpt = polyline[1]
        return anglex((secondpt[0] - orgn[0], secondpt[1] - orgn[1]))
    secondpt = polyline[-2]
    return anglex((secondpt[0] - orgn[0], secondpt[1] - orgn[1]))


def getBoundaryIndex(polyline, boundary, inbound):
    """
    The function returns the position of the last two points inside the terrain boundaries:
    It returns the position of the last point along the outer boundary and the 
    position of the point before the last along the inner boundary
    This is used as a key to sort thalwegs around the virtual pit.

    Parameters
    ----------
    polyline : list of pairs
        polyline defined by a list of (i,j) coordinates.
    boundary : list of pairs
        list of points defining the outer boundary (all points are outside the terrain).
    inbound : list of pairs
        list of points defining the inner boundary (all points are inside the terrain).

    Returns
    -------
    pair of integer
        The index of the points inside the boundary and the inner boundary.

    """
    try:
        return boundary.index(polyline[-1]), inbound.index(polyline[-2])
    except Exception as e:
        print('getBoundaryIndex', polyline[-2:])
        print('Exception message: {!s}'.format(e))


def getNeighbourPoints(nodelist):
    """
    Get all the neighbours around a list of points. Neighbours are unsorted. It
    takes all eight neighbours of each point and remove those already in the list.

    Parameters
    ----------
    nodelist : list of pairs
        list of points around which we compute the neighbours.

    Returns
    -------
    neighbourlist : list of pairs
        Coordinates of the points that are outside but adjacent to the set of points.

    """
    #ldr8 = ((-1, -1), (-1, 0), (-1, 1), (0, 1), (1, 1), (1, 0), (1, -1), (0, -1))
    neighbourset = set()
    for pxl in nodelist:
        (i, j) = pxl
        neighbours = set([(i + di, j + dj) for (di, dj) in ldr8])
        neighbourset |= neighbours
    neighbourset -= set(nodelist)
    neighbourlist = list(neighbourset)
    return neighbourlist


def buildPixelRing(pixellist, nd=4):
    """
    Function sorting a list of points into a ring. If needed, a point can be 
    duplicated meaning that there is an overlap. There must not be any gap in 
    the point list.

    Parameters
    ----------
    pixellist : list of pairs
        coordinates of the points.
    nd : integer, optional
        number of neighbours to consider. Not implemented, we always take four neighbours.

    Returns
    -------
    orderedlist : list of pairs
        the sorted points starting from the top left corner and ending there.

    """
    #ldr4 = ((-1, 0), (0, 1), (1, 0), (0, -1))
    g = deque(ldr4)
    maxrun = 2 * len(pixellist) + 1
    # take the pixel with min i and min j (top left corner)
    startpxlid = min(enumerate(pixellist), key=lambda a: a[1])[0]
    # walk along all pixels until coming back to start
    startpxl = pixellist[startpxlid]
    orderedlist = [startpxl] # stores the result
    nextpxl = (startpxl[0], startpxl[1] + 1)
    if nextpxl not in pixellist:
        print('problem in buildPixelRing')
    pxl = 0
    counter = 0
    # print('nextpxl', nextpxl)
    # print('pixels', pixellist)
    while nextpxl != startpxl and counter <= maxrun:
        orderedlist.append(nextpxl)
        pxl = nextpxl
        i = pxl[0]
        j = pxl[1]
        # coordinates of neighbouring pixels
        pxlaround = [(i + di, j + dj) for (di, dj) in g]
        # print('pxlaround', pxlaround)
        nextipxl = 3
        for ipxl, pxl in enumerate(pxlaround[:-1]):
            if pxl in pixellist:
                nextipxl = ipxl
        nextpxl = pxlaround[nextipxl]
        g.rotate(-nextipxl + 1) # permutation to start with the right neighbour
        # for ipxl, pxl in enumerate(pxlaround):
        #     if pxl in pixellist:
        #         g.rotate(-ipxl+1)
        #         nextpxl = pxl
        #         break
        counter += 1
    #     print('nextpxl', nextpxl)
    # print('orderedlist', orderedlist)
    # counter is there for information: if too many rounds, may be a problem
    if counter > maxrun: 
        print('counter max', maxrun)
    return orderedlist


def getSecondPointPosition(it):
    """
    Function that returns the position of the second point of a line starting from a saddle.
    Since the thalweg id is signed, we return either 1 or -2 depending on
    the id sign.
    
    Parameters
    ----------
    it : integer
        index of the thalweg.

    Returns
    -------
    integer
        index of the second point, either -2 if it is the one before the last or
        1 if it is the second point of the line.

    """
    if it < 0:
        return -2
    return 1

def walkOuterBoundary(image, border, code):
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
    None. image is modified with outer boundary pixels marked with -code

    """
    for i,j in border:
        image[i,j] = code
    (istart, jstart) = min(border)
    istart -= 1
    ldr = list(ldr4)
    #ldr.reverse()
    g = deque(ldr)
    g.rotate(1)
    i = istart
    j = jstart
    b = False
    outcode = -code
    counter = 0
    boundary = [(i,j)]
    # istart, jstart is outside the terrain
    while not b:
#        print(i,j)
        image[i,j] = outcode
        pxlaround = [(i + di, j + dj) for (di, dj) in g]
        candidatepixel = []
        nexti = 4
        for ipxl, (pxli, pxlj) in enumerate(pxlaround):
            if image[pxli, pxlj] != code and image[pxli, pxlj] != outcode:
                candidatepixel.append((pxli, pxlj))
                nexti = ipxl
        if candidatepixel:
            (i,j) = candidatepixel[-1]
            boundary.append((i,j))
        else:
            (i,j) = boundary.pop()
        g.rotate(-nexti + 1)
        b = ((i, j) == (istart, jstart))
        if counter == 2:
            image[istart, jstart] = 0
        counter += 1
#    return boundary

def walkInnerBoundary(image, border, code):
    if border[0] == border[-1]:
        border.pop()
    for i,j in border:
        image[i,j] = code

    for j, pj in enumerate(border):
        pi = border[j-1]
        dp = (pj[0] - pi[0], pj[1] - pi[1])
        discr = dp[0]*dp[1]
        pin = (0,0)
        if discr > 0:
            pin = (pj[0]-dp[1], pj[1])
        elif discr < 0:
            pin = (pj[0], pj[1]-dp[1])
        else:
            pk = border[j-2]
            dpk = (pk[0] - pi[0], pk[1] - pi[1])
            if not (dp[0]*dpk[0]+dp[1]*dpk[1] > 0 and dp[0]*dpk[1]-dp[1]*dpk[0] > 0):
                pin = (pi[0]-dp[1], pi[1]+dp[0])
#                print(pk, pi, pj, dp, dpk, pin)
        if (image[pin] != code):
            image[pin] = -code
        
    vi = [i[0] for i in border]
    vj = [i[1] for i in border]
    mini = min(vi) + 1
    maxi = max(vi)
    minj = min(vj) + 1
    maxj = max(vj)
    for i in range(mini, maxi):
        for j in range(minj, maxj):
            if image[i,j] == 0:
                tmp = [(i + di, j + dj) for (di, dj) in ldr4]
                lst = [image[k] for k in tmp]
                nb = lst.count(-code)
                if nb > 1:
                    image[i,j] = -code

def sideOfLine(line, i, pt):
    p = line[i]
    q = line[i-1]
    r = line[i+1]
    u = (pt[0] - p[0], pt[1] - p[1])
    v1 = (q[0] - p[0], q[1] - p[1])
    v2 = (r[0] - p[0], r[1] - p[1])
    
    sign1 = v1[0]*u[1] - v1[1]*u[0]
    sign2 = v2[0]*u[1] - v2[1]*u[0]
    
    if sign1 > 0 and sign2 < 0:
        return 1
    if sign1 < 0 and sign2 > 0:
        return -1
    
    if v1[0]*u[0] + v1[1]*u[1] > 0:
        if sign1 > 0:
            return 1
        elif sign1 < 0: 
            return -1
        return 0
    if sign2 < 0:
        return 1
    elif sign2 > 0:
        return -1
    return 0

def sideOfWedge(wedge, pt):
    r = wedge[0]
    p = wedge[1]
    q = wedge[2]
    u = (pt[0] - p[0], pt[1] - p[1])
    v1 = (q[0] - p[0], q[1] - p[1])
    v2 = (r[0] - p[0], r[1] - p[1])
    
    sign1 = v1[0]*u[1] - v1[1]*u[0]
    sign2 = v2[0]*u[1] - v2[1]*u[0]
    
    #print(u, v1, v2, sign1, sign2)
    if sign1 > 0 and sign2 < 0:
        return 1
    if sign1 < 0 and sign2 > 0:
        return -1
    
    if v1[0]*u[0] + v1[1]*u[1] > 0:
        if sign1 > 0:
            return 1
        elif sign1 < 0: 
            return -1
        return 0
    if sign2 < 0:
        return 1
    elif sign2 > 0:
        return -1
    return 0