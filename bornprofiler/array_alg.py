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
    
    def curve_overlap(self,curve2):
        if self.axis != curve2.axis:
            raise Exception("Curves do not have the same interpolation axis")
        spacing = max(self.spacing,curve2.spacing)
        valid_min,valid_max = [max(self.min_axis,curve2.min_axis),min(self.max_axis,curve2.max_axis)]
        if valid_min >= valid_max:
            raise Exception("Curves do not overlap on their interpolation axis")
        points = numpy.arange(valid_min,valid_max,spacing)
        numpoints = points.shape[0]
        if numpoints < 2:
            raise Exception("Curves do not share sufficient space on interpolation")
        return [points,numpoints,spacing]


    def __add__(self,addend):
        if type(addend) == curve:
            points, numpoints, spacing = self.curve_overlap(addend)
            result = self.axis_interpolator()(points) + addend.axis_interpolator()(points)
            output_array = numpy.zeros((numpoints,4))
            output_array[:,self.axis] = points
        else:
            points = self.spacial
            numpoints = self.numpoints
            spacing = self.spacing
            result = self.values.flatten() + addend
            output_array = numpy.zeros((numpoints,4))
            output_array[:,(0,2)] = points
        output_array[:,3] = result 
        return curve(output_array,spacing=spacing)

    def __sub__(self,subtrahend):
        if type(subtrahend) == curve:
            points, numpoints, spacing = self.curve_overlap(subtrahend)
            result = self.axis_interpolator()(points) - subtrahend.axis_interpolator()(points)
            output_array = numpy.zeros((numpoints,4))
            output_array[:,self.axis] = points
        else:
            points = self.spacial
            numpoints = self.numpoints
            spacing = self.spacing
            result = self.values.flatten() - subtrahend
            output_array = numpy.zeros((numpoints,4))
            output_array[:,(0,2)] = points

        output_array[:,3] = result
        return curve(output_array,spacing=spacing)


    def __div__(self,divisor):
        if type(divisor) == curve:
            points, numpoints, spacing = self.curve_overlap(divisor)
            result = self.axis_interpolator()(points)/divisor.axis_interpolator()(points)
            output_array = numpy.zeros((numpoints,4))
            output_array[:,self.axis] = points
        else:
            points = self.spacial
            numpoints = self.numpoints
            spacing = self.spacing
            result = self.values.flatten()/divisor
            output_array = numpy.zeros((numpoints,4))
            output_array[:,(0,2)] = points
        output_array[:,3] = result
        return curve(output_array)

    def __mul__(self,multiplier):
        if type(multiplier) == curve:
            points, numpoints, spacing = self.curve_overlap(multiplier)
            result = self.axis_interpolator()(points)*multiplier.axis_interpolator()(points)
            output_array = numpy.zeros((numpoints,4))
            output_array[:,self.axis] = points
        else:
            points = self.spacial
            numpoints = self.numpoints
            spacing = self.spacing
            result = self.values.flatten()*multiplier
            output_array = numpy.zeros((numpoints,4))
            output_array[:,(0,2)] = points
        output_array[:,3] = result
        return curve(output_array)
