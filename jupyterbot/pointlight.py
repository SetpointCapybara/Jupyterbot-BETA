from utils import *


class PointLight:


    #######################################
    # Attributes
    #######################################

    @property
    def name(self):
        """The object name."""
        return self._name

    @property
    def color(self):
        """The object color."""
        return self._color

    @property
    def intensity(self):
        """The light intensity."""
        return self._intensity

    @property
    def position(self):
        """The light position."""
        return self._position

    @property
    def max_distance(self):
        """The light maximum distance (meters)."""
        return self._max_distance

    #######################################
    # Constructor
    #######################################

    def __init__(self, name="genLight", color="white", intensity=1, position=[0,0,0], max_distance=0):

        # Error handling

        if not (str(type(name)) == "<class 'str'>"):
            raise Exception("The parameter 'name' should be a string")

        if not Utils.is_a_color(color):
            raise Exception("The parameter 'color' should be a color")

        if not Utils.is_a_number(intensity) or intensity < 0:
            raise Exception("The parameter 'intensity' should be a nonnegative float")

        if not Utils.is_a_vector(position, 3):
            raise Exception("The parameter 'position' should be a 3D vector")

        if not Utils.is_a_number(max_distance) or max_distance < 0:
            raise Exception("The parameter '_max_distance' should be a nonnegative float")

        # end error handling

        self._name = name
        self._color = color
        self._intensity = intensity
        self._position = position
        self._max_distance = max_distance



    #######################################
    # Std. Print
    #######################################

    def __repr__(self):

        string = "Point light with name '" + self.name + "': \n\n"
        string += " Color: " + str(self.color) + "\n"
        string += " Position: " + str(self.position) + "\n"
        string += " Intensity: " + str(self.intensity) + "\n"
        string += " Maximum distance: " + str(self.max_distance) + "\n"

        return string

    #######################################
    # Methods
    #######################################

    def gen_code(self):
        """Generate code for injection."""
        return '''const var_''' + self._name + ''' = new PointLight(\"'''+self.color+'''\", '''+str(self.intensity)+''','''+str(self.max_distance)+''');
		var_''' + self._name + '''.position.set('''+str(self.position[0])+''','''+str(self.position[1])+''',''' +str(self.position[2])+''');
		scene.add(var_''' + self._name + ''');
		//USER INPUT GOES HERE
		'''
