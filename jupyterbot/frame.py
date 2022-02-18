from utils import *
from ball import *
from meshmaterial import *
import numpy as np



class Frame:
    """
  A ball object.

  Parameters
  ----------
  htm : 4x4 numpy array or 4x4 nested list
      The object's configuration
      (default: the same as the current HTM).

  name : string
      The object's name
      (default: 'genBall').

  radius : positive float
      The object's radius, in meters
      (default: 1).

  mass : positive float
      The object's mass, in kg
      (default: 1).

  color : string
      The object's color, in HTML compatible string
      (default: "red").

  mesh_material: 'MeshMaterial' object
      The object mesh material. If set to 'None', the default is used.
      (default: None).
  """

    #######################################
    # Attributes
    #######################################

    @property
    def size(self):
        """The axis size, in meters."""
        return self._size

    @property
    def name(self):
        """Name of the object."""
        return self._name

    @property
    def htm(self):
        """Object pose. A 4x4 homogeneous transformation matrix written is scenario coordinates."""
        return np.array(self._htm)

    @property
    def axis_color(self):
        """"""
        return self._axis_color

    @property
    def axis_name(self):
        """"""
        return self._axis_name

    #######################################
    # Constructor
    #######################################

    def __init__(self, htm=np.identity(4), name="genFrame", size=0.5, axis_color=["red","lime","blue"], axis_names=['x', 'y', 'z']):

        # Error handling
        if not Utils.is_a_matrix(htm, 4, 4):
            raise Exception("The parameter 'htm' should be a 4x4 homogeneous transformation matrix")

        if not Utils.is_a_number(size) or size <= 0:
            raise Exception("The parameter 'size' should be a positive float")

        if not (str(type(name)) == "<class 'str'>"):
            raise Exception("The parameter 'name' should be a string")

        if not (str(type(axis_color)) == "<class 'list'>") or (not (len(axis_color) == 3)):
            raise Exception("The parameter 'list' should be a a list of 3 HTML-compatible color strings")

        for color in axis_color:
            if not Utils.is_a_color(color):
                raise Exception("The parameter 'list' should be a a list of 3 HTML-compatible color strings")


        # end error handling


        self._htm = np.array(htm)
        self._name = name
        self._size = size
        self._axis_names = axis_names
        self._axis_color = axis_color
        self._ball = Ball(name="dummy_ball_"+name, htm=htm, radius=0.0001, mesh_material=MeshMaterial(opacity=0))

        # Set initial total configuration
        self.set_ani_frame(self._htm)

    #######################################
    # Std. Print
    #######################################

    def __repr__(self):

        string = "Axes with name '" + self.name + "': \n\n"
        string += " HTM: \n" + str(self.htm) + "\n"

        return string

    #######################################
    # Methods
    #######################################

    def add_ani_frame(self, time, htm=None):
        """
    Add a single configuration to the object's animation queue.

    Parameters
    ----------
    time: positive float
        The timestamp of the animation frame, in milliseconds.
    htm : 4x4 numpy array or 4x4 nested list
        The object's configuration
        (default: the same as the current HTM).

    Returns
    -------
    None
    """
        if htm is None:
            htm = self._htm

        # Error handling
        if not Utils.is_a_matrix(htm, 4, 4):
            raise Exception("The parameter 'htm' should be a 4x4 homogeneous transformation matrix")

        if not Utils.is_a_number(time) or time < 0:
            raise Exception("The parameter 'time' should be a positive float")
        # end error handling

        self._ball.add_ani_frame(time, htm)

    # Set config. Restart animation queue
    def set_ani_frame(self, htm=None):
        """
    Reset object's animation queue and add a single configuration to the
    object's animation queue.

    Parameters
    ----------
    htm : 4x4 numpy array or 4x4 nested list
        The object's configuration
        (default: the same as the current HTM).

    Returns
    -------
    None
    """
        if htm is None:
            htm = self._htm

        # Error handling
        if not Utils.is_a_matrix(htm, 4, 4):
            raise Exception("The optional parameter 'htm' should be a 4x4 homogeneous transformation matrix")
        # end error handling

        self._ball.set_ani_frame(htm)

    def gen_code(self):
        """Generate code for injection."""


        string = self._ball.gen_code()
        string = self.code.replace("//USER INPUT GOES HERE", "")
        string += '''var var_axes_''' + self._name + ''' = new AxesHelper('''+str(self.size)+''');
        var_dummy_ball_''' + self._name + '''.shape.add(var_axes_''' + self._name + ''');
        var_axes_''' + self._name + '''.setColors(\''''+self.axis_color[0]+'''\',\''''+self.axis_color[1]+'''\',\''''+self.axis_color[2]+'''\');
		sceneElements.push(var_dummy_ball_''' + self._name + ''');
		//USER INPUT GOES HERE
		'''
        return string


