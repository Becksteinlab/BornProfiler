# BornProfiler -- dealing with boxes
# Copyright (c) 2010 Oliver Beckstein
"""
Simulation boxes --- :mod:`bornprofiler.boxes`
==============================================

Helper functions and classes to process and analyze various simulation
boxes. (Not used?)

"""

import numpy

def normal(a,b):
    n = numpy.cross(a,b)
    return n/numpy.linalg.norm(n)

# corners of a cube
cube = numpy.array([[0,0,0],[1,0,0],[0,1,0],[0,0,1],[1,1,0],[1,0,1],[0,1,1],[1,1,1]])

class Plane(object):
    """Oriented plane."""
    def __init__(self, point, normal):
        self.point = numpy.asarray(point)
        self.normal = numpy.asarray(normal)
    def cmp(self, x):
        """Return +1 if point is above plane, -1 below, 0 in plane."""
        return numpy.sign(numpy.dot(self.normal, numpy.asarray(x) - self.point))
    def __repr__(self):
        return "Plane(%(point)r,%(normal)r)" % vars(self)

class Unitcell(object):
    """Base class for simulation cells"""

class Orthorhombic(Unitcell):
    def __init__(self, center, lengths):
        self.center = numpy.asarray(center)
        self.lengths = numpy.asarray(lengths)
        self.corners = self.center + (cube - 0.5)*self.lengths
        # corners MUST be in the right order, i.e. ordered by cartesian coordinates
        # and x varying fastest, then y, then z
        # See :data:`cube` as example
        x = self.corners
        e = x[1:4] - x[0]  # basis vectors, origin x[0]
        # normal for each face; opposite faces are parallel but oriented
        # oppositely
        normals = [-normal(e[2], e[0]), -normal(e[0], e[1]), -normal(e[1],e[2]),
                    normal(e[2], e[0]),  normal(e[0], e[1]),  normal(e[1],e[2])]
        # point on each face
        points = [x[0], x[0], x[0], x[7], x[7], x[7]]
        self.faces = [Plane(p,n) for p,n in zip(points,normals)]

    def sorted_corners(self):
        # sort corners by distance from first element
        distances = numpy.array([numpy.linalg.norm(x) for x in (self.corners - self.corners[0])])
        return self.corners[distances.argsort()]

    def _point_cmp(self,x):
        return numpy.array([f.cmp(x) for f in self.faces])

    def point_inside(self, x):
        """Check if point *x* is inside the cube.

        This is a simple algorithm that just checks if *x* is "below"
        each bounding plane. For a convex hull this is equivalent to
        checking that the point is inside.

        It also counts points on faces as "inside".
        """
        return self.point_cmp(x) != 1

    def point_cmp(self, x):
        """Compare position of point *x* to the box.

        :Returns: -1 if x is in inside, 0 if x is on a face, 1 if it is outside.
        """
        pos = self._point_cmp(x)
        if numpy.all(pos == -1):
            return -1   # all inside
        elif numpy.all(numpy.logical_or(pos == -1, pos == 0)):
            return 0    # at least one face AND "inside"
        return 1        # outside

    def contains(self, other):
        """Return ``True`` if box *other* is contained in this box."""
        return self > other

    def inside(self, other):
        """Return ``True`` if this box is contained in box *other*"""
        return self < other

    def __cmp__(self, y):
        """x.__cmp__(y) <==> cmp(x,y)

        cmp(x, y) -> integer

        Return -1 if x<y (i.e. x fits into y), 0 if x==y, and +1 if x>y.

        Only works for parallel edges at the moment!!
        """
        try:
            # primitive check... just make sure every corner is in the cube
            return numpy.all([self.point_inside(c) for c in y.corners])
        except AttributeError:
            raise TypeError("Other object must contain a list of points in 'corners'")

    def __repr__(self):
        return "Orthorhombic(%(center)r,%(lengths)r)" % vars(self)
