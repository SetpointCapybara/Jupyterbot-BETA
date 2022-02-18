import numpy as np
from utils import *
from model3d import *
from meshmaterial import *

class Link:

    #######################################
    # Attributes
    #######################################

    @property
    def theta(self):
        """The 'theta' parameter of the Denavit-Hartenberg convention (in rad)"""
        return self._theta

    @property
    def d(self):
        """The 'd' parameter of the Denavit-Hartenberg convention (in meters)"""
        return self._d

    @property
    def a(self):
        """The 'a' parameter of the Denavit-Hartenberg convention (in meters)"""
        return self._a

    @property
    def alpha(self):
        """The 'alpha' parameter of the Denavit-Hartenberg convention (in rad)"""
        return self._alpha

    @property
    def joint_number(self):
        """The joint number in the kinematic chain."""
        return self._joint_number

    @property
    def joint_type(self):
        """The joint type (0=revolute, 1=prismatic)."""
        return self._joint_type

    @property
    def mass(self):
        """The link's mass, in kg."""
        return self._mass

    @property
    def center_shift(self):
        """The link's inertia matrix, in kg m²."""
        return self._center_shift

    @property
    def inertia_matrix(self):
        """The link's inertia matrix, in kg m²."""
        return self._inertia_matrix

    @property
    def col_objects(self):
        """Collection of objects that compose the collision model of this link."""
        return self._col_objects

    @property
    def model_3d(self):
        """The 3d model of the object"""
        return self._model_3d

    #######################################
    # Constructor
    #######################################

    def __init__(self, joint_number, theta, d, alpha, a, joint_type, model_3d,
                 mass=1, center_shift = [0,0,0], inertia_matrix=np.identity(3)):

        # Error handling
        if str(type(joint_number)) != "<class 'int'>" or joint_number < 0:
            raise Exception("The 'joint_number' parameter should be a nonnegative integer")

        if not Utils.is_a_number(theta):
            raise Exception("The 'theta' parameter should be a float")

        if not Utils.is_a_number(d):
            raise Exception("The 'd' parameter should be a float")

        if not Utils.is_a_number(alpha):
            raise Exception("The 'alpha' parameter should be a float")

        if not Utils.is_a_number(a):
            raise Exception("The 'a' parameter should be a float")

        if joint_type != 0 and joint_type != 1:
            raise Exception("The 'joint_type' parameter should be either '0' (rotative) or '1' (prismatic)")

        if not (Utils.get_jupyterbot_type(model_3d) == "jupyterbot.Model3D"):
            raise Exception("The parameter 'model_3d' should be a 'Model3D' object")

        if not Utils.is_a_number(mass) or mass < 0:
            raise Exception("The parameter 'mass' should be a positive float")

        if not Utils.is_a_pd_matrix(np.array(inertia_matrix), 3):
            raise Exception("The parameter 'inertia_matrix' should be a symmetric positive definite 3x3 matrix")

        if not Utils.is_a_vector(center_shift, 3):
            raise Exception("The parameter 'center_shift' should be a 3d vector")

        # Code
        self._joint_number = joint_number
        self._theta = theta
        self._d = d
        self._alpha = alpha
        self._a = a
        self._joint_type = joint_type
        self._col_objects = []
        self._model_3d = model_3d
        self._center_shift = center_shift

        self._mass = mass
        self._inertia_matrix = np.array(inertia_matrix)





    #######################################
    # Std. Print
    #######################################

    def __repr__(self):

        string = "Link " + str(self._joint_number) + "' "

        if self._joint_type == 0:
            string += "with rotative joint:\n\n"
            string += " θ (rad) : [variable] \n"
            string += " d (m)   : " + str(self._d) + " \n"
        if self._joint_type == 1:
            string += "with prismatric joint:\n\n"
            string += " θ (rad) : " + str(self._theta) + " \n"
            string += " d (m)   : [variable] \n"

        string += " α (rad) : " + str(self._alpha) + " \n"
        string += " a (m)   : " + str(self._a) + " \n"
        string += " Link mass (kg): " + str(self._mass) + "\n"
        string += " Link inertia matrix (kg*m²): \n " + str(self._inertia_matrix) + "\n"


        return string

    #######################################
    # Methods
    #######################################

    def attach_col_object(self, obj, htm):
        """
    Attach an object (ball, box or cylinder) into the link as a collision
    object.
    
    Parameters
    ----------
    obj: object
        Object to be verified.

    htm : 4x4 numpy array or 4x4 nested list
        The transformation between the link's HTM and the object's HTM

    """

        # Error handling
        if not Utils.is_a_matrix(htm, 4, 4):
            raise Exception("The parameter 'htm' should be a 4x4 homogeneous transformation matrix")

        if not (Utils.is_a_simple_object(obj)):
            raise Exception("The parameter 'obj' must be either a box, a ball or a cylinder")

        self._col_objects.append([obj, htm])

    def gen_code(self, name):

        name_link = "link_"+str(self.joint_number)+"_"+name

        string = self.model_3d.gen_code(name_link)
        string += "const "+name_link+" = {\n"
        string += "theta: "+str(self.theta)+", \n"
        string += "d: " + str(self.d) + ", \n"
        string += "a: " + str(self.a) + ", \n"
        string += "alpha: " + str(self.alpha) + ", \n"
        string += "jointType: " + str(self.joint_type) + ", \n"
        string += "model3d: object3d_" + name_link + "}; \n \n"

        return string
