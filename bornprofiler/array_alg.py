import numpy
from scipy.interpolate import interp1d


class curve(object):
    """Class for encapsulating the information stored in a group of BP dat file.
    Allows for addition, subtraction, and averaging with appropriate 
    interpolation. Interpolation spacing and axis are provided as kwargs while
    array containing [[x,y,z,E],...] is provided as first arg."""
    def __init__(self,*args,**kwargs):
        self.datarray = args[0]
        self.numpoints = self.datarray.shape[0]
        if self.numpoints < 2:
            raise Exception("At least two points needed to define curve")
        self.spacial = self.datarray[:,0:3]
        self.values = numpy.array([self.datarray[:,-1]]).T
        try:
            self.axis = kwargs['axis']
        except:
            self.axis = 2
        self.axis_points = self.spacial[:,self.axis]
        self.min_axis = numpy.amin(self.axis_points)
        self.max_axis = numpy.amax(self.axis_points)
        self.axis_density = self.calc_axis_density()
        try:
            self.spacing = kwargs['spacing']
        except:
            self.spacing = 1/self.axis_density

    def calc_axis_density(self):
        distance = self.max_axis - self.min_axis
        density = (float(self.numpoints) - 1)/distance
        return density

    def axis_interpolator(self):
        return interp1d(self.axis_points.flatten(),self.values.flatten())

    def __add__(self,curve2):
        if self.axis != curve2.axis:
            raise Exception("Curves to be added do not have the same interpolation axis")
        spacing = max(self.spacing,curve2.spacing)
        valid_min,valid_max = [max(self.min_axis,curve2.min_axis),min(self.max_axis,curve2.max_axis)]
        if valid_min >= valid_max:
            raise Exception("Curves to be added do not overlap on their interpolation axis")
        points = numpy.arange(valid_min,valid_max,spacing)
        numpoints = points.shape[0]
        if numpoints < 2:
            raise Exception("Curves to be added do not share sufficient space on interpolation axis to allow for addition")
        result = self.axis_interpolator()(points) + curve2.axis_interpolator()(points)
        output_array = numpy.zeros((numpoints,4))
        output_array[:,self.axis] = points
        output_array[:,3] = result 
        return curve(output_array,spacing=spacing)

    def __sub__(self,curve2):
        if self.axis != curve2.axis:
            raise Exception("Curves to be subtracted do not have the same interpolation axis")
        spacing = max(self.spacing,curve2.spacing)
        valid_min,valid_max = [max(self.min_axis,curve2.min_axis),min(self.max_axis,curve2.max_axis)]
        if valid_min >= valid_max:
            raise Exception("Curves to be subtracted do not overlap on their interpolation axis")
        points = numpy.arange(valid_min,valid_max,spacing)
        numpoints = points.shape[0]
        if numpoints < 2:
            raise Exception("Curves to be added do not share sufficient space on interpolation axis to allow for addition")
        result = self.axis_interpolator()(points) - curve2.axis_interpolator()(points)
        output_array = numpy.zeros((numpoints,4))
        output_array[:,self.axis] = points
        output_array[:,3] = result
        return curve(output_array,spacing=spacing)


    def __div__(self,number):
        divided_vals = self.values/number
        complete_array = numpy.hstack((self.spacial,divided_vals))
        return curve(complete_array)

    def __mul__(self,number):
        multiplied_vals = self.values * number
        complete_array = numpy.hstack((self.spacial,multiplied_vals))
        return curve(complete_array)
