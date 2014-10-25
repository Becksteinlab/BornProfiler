import numpy
from scipy.interpolate import interp1d


class curve(object):
    """Class for encapsulating the information stored in a group of BP dat file.
    Allows for addition, subtraction, and averaging with appropriate 
    interpolation."""
    def __init__(self,*args,**kwargs):
        self.datarray = args[0]
        self.spacial = self.datarray[:,0:3]
        self.values = numpy.array([self.datarray[:,-1]]).T
        self.z = self.spacial[:,-1]
        self.minz = numpy.amin(self.z)
        self.maxz = numpy.amax(self.z)
        self.zdensity = self.calc_zdensity()
        try:
            self.spacing = kwargs['spacing']
        except:
            self.spacing = 1/self.zdensity

    def calc_zdensity(self):
        distance = self.maxz - self.minz
        density = (self.spacial.shape[0] - 1)/distance
        return density

    def z_interpolator(self):
        return interp1d(self.z.flatten(),self.values.flatten())

    def __add__(self,curve2):
        spacing = max(self.spacing,curve2.spacing)
        valid_range = [max(self.minz,curve2.minz),min(self.maxz,curve2.maxz)]
        points = numpy.arange(valid_range[0],valid_range[1],spacing)
        numpoints = points.shape[0]
        result = self.z_interpolator()(points) + curve2.z_interpolator()(points)
        spacial_array = numpy.hstack((numpy.zeros((numpoints,2)),points.reshape(numpoints,1)))
        complete_array = numpy.hstack((spacial_array,result.reshape(numpoints,1)))
        return curve(complete_array,spacing=spacing)

    def __sub__(self,curve2):
        spacing = max(self.spacing,curve2.spacing)
        valid_range = [max(self.minz,curve2.minz),min(self.maxz,curve2.maxz)]
        points = numpy.arange(valid_range[0],valid_range[1],spacing)
        numpoints = points.shape[0]
        result = self.z_interpolator()(points) - curve2.z_interpolator()(points)
        spacial_array = numpy.hstack((numpy.zeros((numpoints,2)),points.reshape(numpoints,1)))
        complete_array = numpy.hstack((spacial_array,result.reshape(numpoints,1)))
        return curve(complete_array,spacing=spacing)

    def __div__(self,number):
        divided_vals = self.values/number
        complete_array = numpy.hstack((self.spacial,divided_vals))
        return curve(complete_array)
