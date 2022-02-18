from utils import *
import numpy as np
from meshmaterial import *


class Ball:
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
      The object's color, a HTML - compatible string
      (default: "red").

  opacity : float between 0 and 1
      The opacity. 1 = fully opaque, and 0 = transparent.

  mesh_material: 'MeshMaterial' object
      The object mesh material. If set to 'None', the default is used.
      (default: None).
  """

    #######################################
    # Attributes
    #######################################

    @property
    def radius(self):
        """The ball radius, in meters."""
        return self._radius

    @property
    def name(self):
        """Name of the object."""
        return self._name

    @property
    def htm(self):
        """Object pose. A 4x4 homogeneous transformation matrix written is scenario coordinates."""
        return np.array(self._htm)

    @property
    def mass(self):
        """Mass of the object, in kg."""
        return self._mass

    @property
    def color(self):
        """Color of the object, HTML compatible"""
        return self.mesh_material.color

    @property
    def mesh_material(self):
        """Mesh properties of the object"""
        return self._mesh_material

    #######################################
    # Constructor
    #######################################

    def __init__(self, htm=np.identity(4), name="genBall", radius=1, mass=1, color="red", opacity=1, mesh_material=None):

        # Error handling
        if not Utils.is_a_matrix(htm, 4, 4):
            raise Exception("The parameter 'htm' should be a 4x4 homogeneous transformation matrix")

        if not Utils.is_a_number(mass) or mass < 0:
            raise Exception("The parameter 'mass' should be a positive float")

        if not Utils.is_a_number(radius) or radius < 0:
            raise Exception("The parameter 'radius' should be a positive float")

        if not (str(type(name)) == "<class 'str'>"):
            raise Exception("The parameter 'name' should be a string")

        if not Utils.is_a_color(color):
            raise Exception("The parameter 'color' should be a color")

        if not ((mesh_material is None) or (Utils.get_jupyterbot_type(mesh_material) == "jupyterbot.MeshMaterial")):
            raise Exception("The parameter 'mesh_material' should be either 'None' or a 'MeshMaterial' object")

        if (not Utils.is_a_number(opacity)) or opacity < 0 or opacity > 1:
            raise Exception("The parameter 'opacity' should be a float between 0 and 1")
        # end error handling


        self._radius = radius
        self._htm = np.array(htm)
        self._name = name
        self._str_name = "'var_" + name + "'"
        self._mass = 1

        if mesh_material is None:
            self._mesh_material = MeshMaterial(color=color, opacity=1)
        else:
            self._mesh_material = mesh_material

        # Set initial total configuration
        self.set_ani_frame(self._htm)

    #######################################
    # Std. Print
    #######################################

    def __repr__(self):

        string = "Ball with name '" + self.name + "': \n\n"
        string += " Radius (m): " + str(self.radius) + "\n"
        string += " Color: " + str(self.color) + "\n"
        string += " Mass (kg): " + str(self.mass) + "\n"
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

        f = [htm[0][0], htm[0][1], htm[0][2], htm[0][3],
             htm[1][0], htm[1][1], htm[1][2], htm[1][3],
             htm[2][0], htm[2][1], htm[2][2], htm[2][3],
             0, 0, 0, 1, time]

        self._htm = htm
        self._frames.append(f)

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

        self._frames = []
        self.code = ''
        self.add_ani_frame(0,htm)

    def gen_code(self):
        """Generate code for injection."""
        self._str_frames = str(self._frames)
        string = self.mesh_material.gen_code(self._name) + "\n"
        string += '''const var_''' + self._name + ''' = new Ball(''' + self._str_name + ''',''' + str(
            self._radius) + ''',''' + self._str_frames + ''', material_'''+self._name+''');
		sceneElements.push(var_''' + self._name + ''');
		//USER INPUT GOES HERE
		'''
        return string

    # Compute inertia matrix with respect to the inertia frame
    def inertia_matrix(self, htm=None):
        """
    The 3D inertia matrix of the object, written in the 'inertia frame', that is, a frame attached to the center of mass of the object which has the same orientation as the scenario frame.

    Parameters
    ----------
    htm : 4x4 numpy array or 4x4 nested list
        The object's configuration for which the inertia matrix will be computed
        (default: the same as the current HTM).

    Returns
    -------
     inertia_matrix : 3x3 numpy array
        The 3D inertia matrix.
    """

        if htm is None:
            htm = self._htm

        # Error handling
        if not Utils.is_a_matrix(htm, 4, 4):
            raise Exception("The optional parameter 'htm' should be a 4x4 homogeneous transformation matrix")
        # end error handling

        I = (2 / 5) * self._mass * (self._radius * self._radius)

        return I * np.identity(3)

    def copy(self):
        """Return a deep copy of the object, without copying the animation frames."""
        return Ball(self.htm, self.name + "_copy", self.radius, self.mass, self.color)

    # Compute the h projection of a point into an object
    def h_projection(self, point, h=0.001, htm=None):
        """
    The h projection of a point in the object, that is, the 
    h-closest point in the object to a point 'point'.

    Parameters
    ----------
    point : 3D vector
        The point for which the projection will be computed.

    h : positive float
        Smoothing parameter, in meters
        (defalt: 0.001 m)

    htm : 4x4 numpy array or 4x4 nested list
        The object's configuration
        (default: the same as the current HTM).            

    Returns
    -------
     projPoint : 3D vector
        The h-projection of the point 'point' in the object.

     d : positive float
        The h-distance between the object and 'point'.     
    """

        if htm is None:
            htm = self._htm

        # Error handling
        if not Utils.is_a_matrix(htm, 4, 4):
            raise Exception("The optional parameter 'htm' should be a 4x4 homogeneous transformation matrix")

        if not Utils.is_a_vector(point, 3):
            raise Exception("The parameter 'point' should be a 3D vector")

        if not Utils.is_a_number(h) or h <= 0:
            raise Exception("The optional parameter 'h' should be a positive number")
        # end error handling
        tpoint = np.transpose(htm[0:3, 0:3]) @ (point - htm[0:3, 3])

        delta = 0.001
        r = np.linalg.norm(tpoint)

        drf = Utils.fun_Int(r + delta, h, self.radius)
        drb = Utils.fun_Int(r - delta, h, self.radius)

        dr = (drf - drb) / (2 * delta)
        d = 0.5 * (drf + drb)
        ppoint = tpoint - (dr / (r + 0.00001)) * tpoint

        return htm[0:3, 0:3] @ ppoint + htm[0:3, 3], d
