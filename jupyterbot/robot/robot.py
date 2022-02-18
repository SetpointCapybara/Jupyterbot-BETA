import os, sys


currentdir = os.path.dirname(os.path.realpath(__file__))
parentdir = os.path.dirname(currentdir)
sys.path.append(parentdir)

import numpy as np
from utils import *
from links import *
from meshmaterial import *

from ._set_ani_frame import _set_ani_frame
from ._add_ani_frame import _add_ani_frame

from ._ikm import _ikm
from ._fkm import _fkm
from ._jac_geo import _jac_geo

from ._dynM import _dynM
from ._dynC import _dynC
from ._dynG import _dynG
from ._dyn_param import _dyn_param

from ._vector_field import _vector_field
from ._task_function import _task_function
from ._coop_task_function import _coop_task_function

from ._control_demo_1 import _control_demo_1
from ._control_demo_2 import _control_demo_2
from ._control_demo_3 import _control_demo_3

from ._gen_code import _gen_code
from ._update_col_object import _update_col_object
from ._add_col_object import _add_col_object
from ._attach_object import _attach_object
from ._detach_object import _detach_object

from ._compute_hg_dist import _compute_hg_dist
from .robotEnv import RobotEnv

from ._create_kukakr5 import _create_kukakr5


class Robot:
    """
  A class that contains a robot object in jupyterbot.

  Parameters
  ----------
  htm : 4x4 numpy array or 4x4 nested list
      The robot base's configuration
      (default: 4x4 identity matrix).
  name : string
      The robot name 
      (default: 'genRobot').

  """

    #######################################
    # Attributes
    #######################################

    @property
    def q(self):
        """The current joint configuration."""
        return np.array(self._q)

    @property
    def q0(self):
        """The default joint configuration."""
        return np.array(self._q0)

    @property
    def htm(self):
        """
        The current base configuration.
        A 4x4 homogeneous matrix written is scenario coordinates.
        """
        return np.array(self._htm)

    @property
    def links(self):
        """Data structures containing the links of the robot."""
        return self._links

    @property
    def attached_objects(self):
        """Data structures containing the objects attached into the robot."""
        return self._attached_objects

    @property
    def name(self):
        """Name of the object."""
        return self._str_name

    @property
    def object_3d_base(self):
        """3d object of the base"""
        return self._object_3d_base

    #######################################
    # Constructor
    #######################################

    def __init__(self, name, base_3d_obj, links, htm=np.identity(4), q0=None):
        # Error handling


        # end error handling

        self._frames = []
        self._object_3d_base = base_3d_obj
        self._htm = np.array(htm)
        self._str_name = name
        self._attached_objects = []
        self._links = links

        n = len(links)

        if not(q0 is None):
            self._q0 = q0
        else:
            self._q0 = np.zeros((n,))

        # Set initial total configuration
        self.set_ani_frame(self._q0, self._htm)

    #######################################
    # Std. Print
    #######################################

    def __repr__(self):
        n = len(self.links)

        string = "Robot with name '" + self._str_name + "': \n\n"
        string += " Number of joints: " + str(n) + "\n"
        string += " Joint types: "

        for i in range(n):
            string += "R" if self._links[i].joint_type == 0 else "P"

        string += "\n"
        string += " Current configuration: " + str([round(num, 3) for num in np.ndarray.tolist(self._q)]) + "\n"
        string += " Current base HTM: \n" + str(self._htm) + "\n"
        string += " Current end-effector HTM: \n" + str(self.fkm())
        return string

    #######################################
    # Methods for configuration changing
    #######################################

    def add_ani_frame(self, time, q=None, htm=None):
        """
    Add a single configuration to the object's animation queue.

    Parameters
    ----------
    time: positive float
        The timestamp of the animation frame, in milliseconds.
    q : nd numpy vector or array
        The manipulator's joint configuration.
    htm : 4x4 numpy array or 4x4 nested list
        The robot base's configuration
        (default: the same as the current HTM).

    Returns
    -------
    None
    """
        return _add_ani_frame(self, time, q, htm)

    def set_ani_frame(self, q=None, htm=None):
        """
    Reset object's animation queue and add a single configuration to the 
    object's animation queue.

    Parameters
    ----------
    q : nd numpy vector or array
        The manipulator's joint configuration 
        (default: the default joint configuration for the manipulator, q0).
    htm : 4x4 numpy array or 4x4 nested list
        The robot base's configuration
        (default: the same as the current HTM).

    Returns
    -------
    None
    """
        return _set_ani_frame(self, q, htm)

    #######################################
    # Methods for kinematics model
    #######################################

    def fkm(self, q=None, axis='eef', htm=None):
        """
    Compute the forward kinematics for an axis at a given joint and base
    configuration. Everything is written in the scenario coordinates.

    Parameters
    ----------
    q : nd numpy vector or array
        The manipulator's joint configuration 
        (default: the default joint configuration for the manipulator).
    axis : string
        For which axis you want to compute the FK:
        'eef': for the end-effector
        'dh': for all Denavit-Hartenberg axis
        'com': for all center-of-mass axis
        (default: 'eef').    
    htm : 4x4 numpy array or 4x4 nested list
        The robot base's configuration
        (default: the same as the current HTM).

    Returns
    -------
    FK : 4x4 or nx4x4 numpy matrix
        For axis='eef', returns a single htm. For the other cases, return
        n htms as a nx4x4 numpy matrix.
    """
        return _fkm(self, q, axis, htm)

    def ikm(self, htm_target, q0=None, p_tol=0.005, a_tol=5, no_iter_max=2000):
        """
    Try to solve the inverse kinematic problem for the end-effector, given a
    desired homogeneous transformation matrix. It returns the manipulator
    configuration.

    Important: it disregards the current htm of the base of the robot. That is,
    it assumes that robot.htm = np.identity(4). You can easily consider other
    cases by transforming htm_target as Utils.inv_htm(robot.htm) @ htm_target.

    Uses an iterative algorithm.

    The algorithm can fail, throwing an Exception when it happens.

    Parameters
    ----------
    htm_target : 4x4 numpy array or 4x4 nested list
        The target end-effector HTM, written in scenario coordinates.
    q0 : nd numpy vector or array
        Initial guess for the algorithm for the joint configuration
        (default: a random joint configuration).
    p_tol : positive float
        The accepted error for the end-effector position, in meters
        (default: 0.005 m).    
    a_tol : positive float
        The accepted error for the end-effector orientation, in degrees
        (default: 5 degrees). 
    no_iter_max : positive int
        The maximum number of iterations for the algorithm
        (default: 2000 iterations). 

    Returns
    -------
    q : nd numpy vector or array
        The configuration that solves the IK problem
    """
        return _ikm(self, htm_target, q0, p_tol, a_tol, no_iter_max)

    def jac_geo(self, q=None, axis='eef', htm=None):
        """
    Compute the geometric Jacobian for an axis at a given joint and base
    configuration. Also returns the forward kinematics as a by-product.
    Everything is written in the scenario coordinates.

    Parameters
    ----------
    q : nd numpy vector or array
        The manipulator's joint configuration 
        (default: the default joint configuration for the manipulator).
    axis : string
        For which axis you want to compute the FK:
        'eef': for the end-effector
        'dh': for all Denavit-Hartenberg axis
        'com': for all center-of-mass axis
        (default: 'eef').    
    htm : 4x4 numpy array or 4x4 nested list
        The robot base's configuration 
        (default: the same as the current htm).

    Returns
    -------
    J : 6xn or nx6xn numpy matrix
        For axis='eef', returns a single 6xn Jacobian. For the other cases, 
        return n Jacobians as a nx6xn numpy matrix.

    FK : 4x4 or nx4x4 numpy matrix
        For axis='eef', returns a single htm. For the other cases, return
        n htms as a nx4x4 numpy matrix.
    """
        return _jac_geo(self, q, axis, htm)

    #######################################
    # Methods for dynamics model
    #######################################

    def dynM(self, q=None):
        """
    Compute the generalized inertia matrix at a given joint configuration.

    Parameters
    ----------
    q : nd numpy vector or array
        The manipulator's joint configuration, a nD numpy vector or nD array 
        (default: the default joint configuration for the manipulator).

    Returns
    -------
    M : nxn numpy array
        The generalized inertia matrix at the joint configuration q.
    """
        return _dynM(self, q)

    def dynC(self, qdot, q=None):
        """
    Compute the generalized Coriolis-Centrifugal torques at a given joint 
    configuration speed and joint configuration.

    Parameters
    ----------
    qdot : nd numpy vector or array
        The manipulator's joint configuration speed. 

    q : nd numpy vector or array
        The manipulator's joint configuration 
        (default: the default joint configuration for the manipulator).

    Returns
    -------
    C : nD numpy vector
        The generalized Coriolis-Centrifugal torques at the joint 
        configuration q and joint configuration speed qdot.
    """
        return _dynC(self, qdot, q)

    def dynG(self, q=None):
        """
    Compute the generalized gravity torques at a given joint configuration.

    Parameters
    ----------
    q : nd numpy vector or array
        The manipulator's joint configuration, a nD numpy vector or nD array 
        (default: the default joint configuration for the manipulator).

    Returns
    -------
    G : nD numpy vector
        The generalized gravity torques at the joint configuration q.
    """
        return _dynG(self, q)

    def dyn_param(self, qdot, q=None):
        """
    Compute all the three dynamic parameters at a given joint configuration
    and joint configuration speed. It is more efficient than calling them
    separatedly.

    Parameters
    ----------
    qdot : nd numpy vector or array
        The manipulator's joint configuration speed .

    q : nd numpy vector or array
        The manipulator's joint configuration 
        (default: the default joint configuration for the manipulator).

    Returns
    -------
    M : nxn numpy array
        The generalized inertia matrix at the joint configuration q.

    G : nD numpy vector
        The generalized gravity torques at the joint configuration q.

    C : nD numpy vector
        The generalized Coriolis-Centrifugal torques at the joint 
        configuration q and joint configuration speed qdot.
    """
        return _dyn_param(self, qdot, q)

    #######################################
    # Methods for control
    #######################################
    @staticmethod
    def vector_field(curve, alpha=1, const_vel=1):
        """
    Computes a handle to a vector field function fun(p). Uses the vector field
    presented in 
    
    "Adriano M. C. Rezende; Vinicius M. Goncalves; Luciano C. A. Pimenta: 
    Constructive Time-Varying Vector Fields for Robot Navigation 
    IEEE Transactions on Robotics (2021)". 
    
    The vector field has constant velocity and use the function 
    G(u) = (2/pi)*atan(alpha*u).


    Parameters
    ----------
    curve : nxm numpy array or nxm nested list
        Curve, described as sampled points. Each one of the n rows should 
        contain a m-dimensional float vector that is the n-th m-dimensional
        sampled point of the curve. 
 
    alpha : positive float
        Controls the vector field behaviour. Greater alpha's imply more 
        robustness to the vector field, but increased velocity and acceleration
        behaviours. Used in G(u) = (2/pi)*atan(alpha*u)
        (default: 1).

    const_vel : positive float
        The constant velocity of the vector field. The signal of this number 
        controls the direction of rotation 
        (default: 1).

    Returns
    -------
    fun: a function handle that you can call as f(p). 'p' should be a
         m-dimensional vector.
    """

        return _vector_field(curve, alpha, const_vel)

    def task_function(self, htm_des, q=None, htm=None):
        """
    Computes the 6-dimensional task function for end-effector pose control,  
    given a joint configuration, a base configuration and the desired pose 
    'htm_des'.

    The first three entries are position error, and the three last entries are
    orientation error.

    Everything is written in scenario coordinates. 

    Also returns the Jacobian of this function.

    Parameters
    ----------
    htm_des : 4x4 numpy array or 4x4 nested list
        The desired end-effector pose. 
 
    q : nd numpy vector or array
        The manipulator's joint configuration 
        (default: the default joint configuration for the manipulator).

    htm : 4x4 numpy array or 4x4 nested list
        The robot base's configuration 
        (default: the same as the current htm).

    Returns
    -------
    r : 6-dimensional numpy vector
        The task function.

    Jr : 6xn numpy matrix
        The respective task Jacobian.
    """
        return _task_function(self, htm_des, q, htm)

    @staticmethod
    def coop_task_function(robot_a, robot_b, htm_a_des, htm_a_b_des, q_a=None, q_b=None):
        """
    Computes the 12-dimensional task function for end-effector pose control
    of two robots 'robot_a', 'robot_b' given their respectives configurations
    q_a and q_b.

    
    The first three components are relative position error.  

    The second three components are relative orientation error.  

    The third three components are position error for 'robot_a'.

    The third three components are orientation error for 'robot_a'.
    
    Everything is written in scenario coordinates. 
    
    Also returns the Jacobian of this function.

    Parameters
    ----------
    robot_a :robot object
        The first robot.

    robot_b :robot object
        The second robot.

    htm_a_des :4x4 numpy array or 4x4 nested list
        The desired pose for the end-effector of 'robot_a'.

    htm_a_b_des :4x4 numpy array or 4x4 nested list
        The desired relative pose between the end-effector of 'robot_a' and
        'robot_b'. That is, inv(htmA) * htmB.

    q_a : nd numpy vector or array
        'robot_a'' joint configuration
        (default: the default joint configuration for 'robot_a').

    q_b : md numpy vector or array
        'robot_b'' joint configuration
        (default: the default joint configuration for 'robot_b').

    Returns
    -------
    r : 12-dimensional numpy vector
        The task function.

    Jr : 12x(n+m) numpy matrix
        The respective task Jacobian.
    """
        return _coop_task_function(robot_a, robot_b, htm_a_des, htm_a_b_des, q_a, q_b)

    #######################################
    # Control demos
    #######################################

    @staticmethod
    def control_demo_1():
        """
    Show a robotic manipulator drawing a circle in a drawboard.
    The task demands that the manipulator keep a constant pose while
    drawing.
    """
        return _control_demo_1()

    @staticmethod
    def control_demo_2():
        """
    Show two robotic manipulators cooperating to transport an object between
    two plates. Note that the relative pose between the two end-effectors must
    be kept.

    The control uses task priority framework using the null space of the task
    Jacobian. The task of keeping the relative pose is prioritized over the
    task of moving the plate.
    """
        return _control_demo_2()

    @staticmethod
    def control_demo_3():
        """
    Show a robot manipulator that must achieve a pose while avoiding an
    obstacle.

    Perform second order kinematic control with constraints that enforce
    collision avoidance.
    """
        return _control_demo_3()

    #######################################
    # Methods for simulation
    #######################################

    def gen_code(self):
        """Generate code for injection."""
        return _gen_code(self)

    def update_col_object(self):
        """
        Update internally the objects that compose the collision model to the
        current configuration of the robot.
        """
        _update_col_object(self)

    def add_col_object(self, sim):
        """
        Add the objects that compose the collision model to a simulation.

        Parameters
        ----------
        sim : 'Simulation' object
            'Simulation' object.
    """
        _add_col_object(self, sim)

    def attach_object(self, obj):
        """
        Attach an object (ball, box or cylinder) to the end-effector.

        Parameters
        ----------
        obj : 'Ball', 'Box' or 'Cylinder' object
            Object to be attached.
    """
        _attach_object(self, obj)

    def detach_object(self, obj):
        """
        Detach an object (ball, box or cylinder) to the end-effector.

        Parameters
        ----------
        obj : 'Ball', 'Box' or 'Cylinder' object
            Object to be detached.
    """
        _detach_object(self, obj)

    #######################################
    # Robot constructors
    #######################################

    @staticmethod
    def create_kukakr5(htm=np.identity(4), name='kukakr5', color="#df6c25", opacity=1):
        """
    Create a Kuka KR5, a six-degree of freedom manipulator.

    Parameters
    ----------
    htm : 4x4 numpy array or 4x4 nested list
        The initial base configuration for the robot. 
 
    name : string
        The robot name
        (default: 'sdofRobot').

    htm : 4x4 numpy array or 4x4 nested list
        The robot base's configuration 
        (default: the same as the current htm).

    Returns
    -------
    robot : Robot object
        The robot.

    """
        base_3d_obj, links, q0 = _create_kukakr5(htm, name, color, opacity)
        return Robot(name, base_3d_obj, links, htm, q0)

    #######################################
    # Advanced methods
    #######################################

    def compute_hg_dist(self, obj, h=0.001, g=0.001, q=None, htm=None, old_dist_struct=None, tol=0.0005, no_iter_max=20):
        """
    Compute (h,g) distance structure from each one of the robot's link to a 
    'simple' external object (ball, box or cylinder), given a joint and base 
    configuration.

    Use an iterative algorithm, based no h-projections and a modification of 
    Von Neumann's cyclic projection algorithm.

    Parameters
    ----------
    obj : a simple object (ball, box or cylinder)
        The external object for which the distance structure is going to be 
        computed, for each robot link.
    h : positive float
        Smoothing parameter for the robot's links, in meters
        (default: 0.001 m).
    g : positive float
        Smoothing parameter for the external object
        (default: 0.001 m).
    q : nd numpy vector or array
        The manipulator's joint configuration 
        (default: the default joint configuration for the manipulator).
    htm : 4x4 numpy array or 4x4 nested list
        The robot base's configuration 
        (default: the same as the current htm).
    old_dist_struct : 'DistStruct' object
        'DistStruct' obtained previously for the same robot and external object.
        Can be used to enhance the algorith speed using the previous closest 
        point as an initial guess
        (default: None).
    tol : positive float
        Tolerance for convergence in the iterative algorithm, in meters.
        (default: 0.0005 m).        
    no_iter_max : positive int
        The maximum number of iterations for the algorithm
        (default: 20 iterations). 

    Returns
    -------
    distStruct : 'DistStruct' object
        Distance struct for each one of the m objects that compose the robot's
        collision model. Contains m dictionaries. Each one of these dictionaries
        contains the following entries:
        
        'linkNumber', containing the robot link index for this entries

        'linkColObjNumber', containing the link index for the covering object 
          at that link

        'hgDistance', containing the smoothed hg distance

        'pObj', containing the closest point at the external object

        'pObjCol', containing the closest point at the covering object in the 
          link
        
        'jacDist', containing the Jacobian (gradient) of this distance function.       
    """

        return _compute_hg_dist(self, obj, h, g, q, htm, old_dist_struct, tol, no_iter_max)
