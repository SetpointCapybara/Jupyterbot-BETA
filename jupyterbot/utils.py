from math import *
import numpy as np
from colour import Color
import plotly.express as px
import httplib2

class Utils:

    """A library that contains some utilities for jupyterbox."""

    #######################################
    # Constants
    #######################################

    _PI = 3.1415926
    _SQRTHALFPI = 1.2533141
    _SQRT2 = 1.4142135
    _CONSTJA = 2.7889
    _CONSTI0HAT1 = 0.24273
    _CONSTI0HAT2 = 0.43023

    #######################################
    # Basic functions
    #######################################

    @staticmethod
    def S(v):
        """
      Returns a 3x3 matrix that implements the cross product for a 3D vector  
      as a matricial product, that is, a matrix S(v) such that for any other 
      3D column  vector w, S(v)w = cross(v,w).
      
      Parameters
      ----------
      v : a 3D vector
          The vector for which the S matrix will be created.

      Returns
      -------
      S : 3x3 numpy matrix
          A matrix that implements the cross product with v.
      """
        vv = np.reshape(v, (3,))
        return np.array([[0, -vv[2], vv[1]],
                         [vv[2], 0, -vv[0]],
                         [-vv[1], vv[0], 0]])

    @staticmethod
    def rot(axis, angle):
        """
      Homogeneous transformation matrix that represents the rotation of an
      angle in an axis.
      
      Parameters
      ----------
      axis : a 3D vector
          The axis of rotation.
      
      angle: float
          The angle of rotation, in radians.

      Returns
      -------
      htm : 4x4 numpy matrix
          The homogeneous transformation matrix.
      """
        a = np.reshape(axis, (3,))
        a = a / np.linalg.norm(a)
        K = Utils.S(a)
        Q = np.identity(3) + sin(angle) * K + (1 - cos(angle)) * (K @ K)
        return np.hstack([np.vstack([Q, np.array([0, 0, 0])]), np.array([[0], [0], [0], [1]])])

    @staticmethod
    def trn(vector):
        """
      Homogeneous transformation matrix that represents the displacement
      of a vector
      
      Parameters
      ----------
      vector : a 3D vector
          The displacement vector.
      
      Returns
      -------
      htm : 4x4 numpy matrix
          The homogeneous transformation matrix.
      """
        v = np.reshape(vector, (3,))
        return np.array([[1, 0, 0, v[0]],
                         [0, 1, 0, v[1]],
                         [0, 0, 1, v[2]],
                         [0, 0, 0,   1]])

    @staticmethod
    def rotx(angle):
        """
      Homogeneous transformation matrix that represents the rotation of an
      angle in the 'x' axis.
      
      Parameters
      ----------
      angle: float
          The angle of rotation, in radians.

      Returns
      -------
      htm : 4x4 numpy matrix
          The homogeneous transformation matrix.
      """
        return np.array([[1, 0, 0, 0],
                         [0, cos(angle), -sin(angle), 0],
                         [0, sin(angle), cos(angle), 0],
                         [0, 0, 0, 1]])

    @staticmethod
    def roty(angle):
        """
      Homogeneous transformation matrix that represents the rotation of an
      angle in the 'y' axis.
      
      Parameters
      ----------
      angle: float
          The angle of rotation, in radians.

      Returns
      -------
      htm : 4x4 numpy matrix
          The homogeneous transformation matrix.
      """
        return np.array([[cos(angle), 0, sin(angle), 0],
                         [0, 1, 0, 0],
                         [-sin(angle), 0, cos(angle), 0],
                         [0, 0, 0, 1]])

    @staticmethod
    def rotz(angle):
        """
      Homogeneous transformation matrix that represents the rotation of an
      angle in the 'z' axis.
      
      Parameters
      ----------
      angle: float
          The angle of rotation, in radians.

      Returns
      -------
      htm : 4x4 numpy matrix
          The homogeneous transformation matrix.
      """
        return np.array([[cos(angle), -sin(angle), 0, 0],
                         [sin(angle), cos(angle), 0, 0],
                         [0, 0, 1, 0],
                         [0, 0, 0, 1]])

    @staticmethod
    def inv_htm(htm):
        """
      Given an homogeneous transformation matrix, compute its inverse.
      
      Parameters
      ----------
      htm: 4X4 numpy array or nested list 
          Homogeneous transformation matrix of the rotation.

      Returns
      -------
      inv_htm: 4X4 numpy array
          The inverse of the transformation matrix.       
      """
        Q = htm[0:3, 0:3]
        p = htm[0:3, 3]

        inv_htm = np.zeros((4, 4))
        inv_htm[0:3, 0:3] = np.transpose(Q)
        inv_htm[0:3, 3] = - np.transpose(Q) @ p
        inv_htm[3, 3] = 1

        return inv_htm

    @staticmethod
    def axis_angle(htm):
        """
      Given an homogeneous transformation matrix representing a rotation, 
      return the rotation axis angle.
      
      Parameters
      ----------
      htm: 4X4 numpy array or nested list 
          Homogeneous transformation matrix of the rotation.

      Returns
      -------
      axis : 3D numpy vector
          The rotation axis.

      angle : float
          The rotation angle, in radians.        
      """
        Q = htm[0:3, 0:3]
        trace = Q[0, 0] + Q[1, 1] + Q[2, 2]
        angle = acos((trace - 1) / 2)
        G = Q @ Q - 2 * cos(angle) * Q + np.identity(3)
        ok = False
        while not ok:
            v = np.random.uniform(-100, 100, size=(3,))
            w = np.random.uniform(-100, 100, size=(3,))
            r = G @ v
            nr = np.linalg.norm(r)
            prod = np.transpose(w) @ r
            if nr > 0.01:
                ortr = w - prod * r / (nr * nr)
                axis = Utils.S(ortr) @ (Q @ ortr)
                naxis = np.linalg.norm(axis)
                ok = naxis > 0.01

        axis = axis / naxis
        return axis, angle

    @staticmethod
    def dp_inv(A, eps):
        """
      Compute the damped pseudoinverse of the matrix A.
      
      Parameters
      ----------
      A: nxm numpy array
          The matrix to compute the damped pseudoinverse.
      
      eps: positive float
          The damping factor.

      Returns
      -------
      pinvA: mxn numpy array
          The damped pseudoinverse of A.    
      """
        n = len(A[0, :])
        return np.linalg.inv(np.transpose(A) @ A + eps * np.identity(n)) @ np.transpose(A)

    @staticmethod
    def jac(f, x, delta=0.0001):
        """
      Compute the numerical Jacobian of a function f at the point x.
      Uses centralized finite difference to compute the derivatives.
      
      Parameters
      ----------
      f: function handle
          The function handle. It should accepts a single 'm' dimensional vector
          and returns a single 'n' dimensional vector. If 'f' outputs more than
          one parameter, all others are ignored.
      
      x: m dimensional numpy vector or list
          Point in which the Jacobian will be computed.

      delta: float
          Variation used in the numerical differentiation
          (default: 0.0001)

      Returns
      -------
      jac: nxm numpy array
          The numerical Jacobian of f at point x.   
      """

        if not (str(type(f)) == "<class 'function'>"):
            raise Exception("The parameter 'f' must be a function that maps 'm' dimensional vectors" \
                            ", in which 'm' is the size of the vector 'x', into 'n' dimensional vectors")

        if not Utils.is_a_vector(x):
            raise Exception("Parameter 'x' should be a vector")

        if str(type(x)) == "<class 'numpy.ndarray'>":
            m = max(np.shape(x))

        if str(type(x)) == "<class 'list'>":
            m = len(x)

        try:
            y = f(x)[0]
            if str(type(y)) == "<class 'numpy.ndarray'>":
                n = max(np.shape(y))

            if str(type(y)) == "<class 'list'>":
                n = len(y)
        except:
            raise Exception("The parameter 'f' must be a function that maps 'm' dimensional vectors" \
                            ", in which 'm' is the size of the vector 'x', into 'n' dimensional vectors")

        jac = np.zeros((n, m))
        idm = np.identity(m)

        for j in range(m):
            yp = np.array(f(x + idm[:, j] * delta)[0]).reshape((n,))
            yn = np.array(f(x - idm[:, j] * delta)[0]).reshape((n,))
            jac[:, j] = (yp - yn) / (2 * delta)

        return jac

    @staticmethod
    def signed_shape_pow(x, x_bar, n):
        """
      Compute the following function, for x>=0

      x^n/(x_bar^(n-1)*n) for x <= x_bar
      x + (1/n-1)*x_bar   for x >  x_bar

      The function is continuous. It is also differentiable once, except for
      n<1 at x=0.

      For x<=0, returns -signed_shape_pow(abs(x), x_bar, n), so the function
      is odd.

      This function is useful in control.

      Parameters
      ----------
      x: float
          Argument of the function.

      x_bar: positive float
          Parameter of the function.

      n: positive float
          Parameter of the function.

      Returns
      -------
      y: float
          The return of the function.
      """

        if x >= 0:
            if x <= x_bar:
                return pow(x, n)/(pow(x_bar, n - 1) * n)
            else:
                return x+(1/n - 1) * x_bar
        else:
            return -Utils.signed_shape_pow(abs(x), x_bar, n)

    #######################################
    # Type check functions
    #######################################

    @staticmethod
    def is_a_number(obj):
        """
      Check if the argument is a float or int number
      
      Parameters
      ----------
      obj: object
          Object to be verified.
      
      Returns
      -------
      is_type: boolean
          If the object is of the type.   
      """

        return str(type(obj)) == "<class 'int'>" or str(type(obj)) == "<class 'float'>"

    @staticmethod
    def is_a_natural_number(obj):
        """
      Check if the argument is a natural number (integer and >=0)
      
      Parameters
      ----------
      obj: object
          Object to be verified.
      
      Returns
      -------
      is_type: boolean
          If the object is of the type.   
      """

        return str(type(obj)) == "<class 'int'>" and obj >= 0

    @staticmethod
    def is_a_matrix(obj, n=None, m=None):
        """
      Check if the argument is a nxm matrix of floats.
      
      Parameters
      ----------
      obj: object
          Object to be verified.

      n: positive int
          Number of rows
          (default: it does not matter).

      m: positive int
          Number of columns
          (default: it does not matter).

      Returns
      -------
      is_type: boolean
          If the object is of the type.   
      """

        value = str(type(obj)) == "<class 'numpy.ndarray'>" or str(type(obj)) == "<class 'list'>"

        if str(type(obj)) == "<class 'numpy.ndarray'>":
            shp = np.shape(obj)
            if len(shp) == 1:
                if (not n is None) and (not m is None):
                    value = (shp[0] == n and m == 1) or (shp[0] == m and n == 1)

                if (not n is None) and (m is None):
                    value = (shp[0] == n) or (n == 1)

                if (n is None) and (not m is None):
                    value = (m == 1) or (shp[0] == m)

                if (n is None) and (m is None):
                    value = True

            if len(shp) == 2:
                if (not n is None) and (not m is None):
                    value = shp[0] == n and shp[1] == m

                if (not n is None) and (m is None):
                    value = shp[0] == n

                if (n is None) and (not m is None):
                    value = shp[1] == m

                if (n is None) and (m is None):
                    value = True

            if len(shp) > 2:
                value = False

        if str(type(obj)) == "<class 'list'>":
            return Utils.is_a_matrix(np.array(obj), n, m)
        else:
            return value

    @staticmethod
    def is_a_vector(obj, n=None):
        """
      Check if the argument is a n vector of floats.
      
      Parameters
      ----------
      obj: object
          Object to be verified.

      n: positive int
          Number of elements
          (default: it does not matter).

      Returns
      -------
      is_type: boolean
          If the object is of the type.   
      """
        return Utils.is_a_matrix(obj, n, 1) or Utils.is_a_matrix(obj, 1, n)

    @staticmethod
    def is_a_pd_matrix(obj, n=None):
        """
      Check if the argument is a symmetric nxn positive (semi)-definite matrix.
      
      Parameters
      ----------
      obj: object
          Object to be verified.

      n: positive int
          Dimension of the square matrix
          (default: it does not matter).
    
      Returns
      -------
      is_type: boolean
          If the object is of the type.   
      """
        value = Utils.is_a_matrix(obj, n, n)

        if value:
            value = np.allclose(obj, np.transpose(obj), rtol=1e-05, atol=1e-08)

        if value:
            try:
                np.linalg.cholesky(obj)
            except:
                value = False

        return value

    @staticmethod
    def is_a_color(obj):
        """
      Check if the argument is a HTML-compatible string that represents a color.
      
      Parameters
      ----------
      obj: object
          Object to be verified.
      
      Returns
      -------
      is_type: boolean
          If the object is of the type.   
      """

        try:
            obj = obj.replace(" ", "")
            Color(obj)
            return True
        except ValueError:
            return False

    @staticmethod
    def get_jupyterbot_type(obj):
        """
      Return the jupyterbot type of the object. 
      Return the empty string if it is not a jupyterbot object.
      
      Parameters
      ----------
      obj: object
          Object to be verified.
      
      Returns
      -------
      obj_type: string
          jupyterbot type.   
      """
        type_str = str(type(obj))
        ind = type_str.find("jupyterbot.")
        if ind == -1:
            ind = type_str.find("cylinder.")
        if ind == -1:
            ind = type_str.find("box.")
        if ind == -1:
            ind = type_str.find("sphere.")
        if ind == -1:
            ind = type_str.find("robot.")
        if ind == -1:
            ind = type_str.find("simulation.")
        if ind == -1:
            ind = type_str.find("meshmaterial.")
        if ind == -1:
            ind = type_str.find("texture.")
        if ind == -1:
            ind = type_str.find("pointlight.")
        if ind == -1:
            ind = type_str.find("frame.")
        if ind == -1:
            ind = type_str.find("model3d.")

        if ind == -1:
            return ""
        else:
            ind1 = type_str.rfind('.')
            ind2 = type_str.rfind('>')
            return "jupyterbot." + type_str[ind1 + 1:ind2 - 1]

    @staticmethod
    def is_a_simple_object(obj):
        """
      Check if the argument is a ball, box or cylinder.
      
      Parameters
      ----------
      obj: object
          Object to be verified.
      
      Returns
      -------
      is_type: boolean
          If the object is of the type.   
      """
        return (Utils.get_jupyterbot_type(obj) == "jupyterbot.Ball") or \
               (Utils.get_jupyterbot_type(obj) == "jupyterbot.Box") or \
               (Utils.get_jupyterbot_type(obj) == "jupyterbot.Cylinder")

    @staticmethod
    def is_a_obj_sim(obj):
        """
      Check if the argument is an object that can be put into the simulator.
      
      Parameters
      ----------
      obj: object
          Object to be verified.
      
      Returns
      -------
      is_type: boolean
          If the object is of the type.   
      """
        return (Utils.get_jupyterbot_type(obj) == "jupyterbot.Ball") or \
               (Utils.get_jupyterbot_type(obj) == "jupyterbot.Box") or \
               (Utils.get_jupyterbot_type(obj) == "jupyterbot.Cylinder") or \
               (Utils.get_jupyterbot_type(obj) == "jupyterbot.Robot") or \
               (Utils.get_jupyterbot_type(obj) == "jupyterbot.PointLight") or \
               (Utils.get_jupyterbot_type(obj) == "jupyterbot.Frame")

    @staticmethod
    def is_url_available(url, types):

        if not (str(type(url)) == "<class 'str'>"):
            return " is not a valid url"

        ind = url.rfind('.')
        filetype = url[ind+1:]

        if not(filetype in types):
            return " must be an file of the types: '"+str(types)+"'"

        try:
            h = httplib2.Http()
            resp = h.request(url, 'HEAD')
            if int(resp[0]['status']) < 400:
                return "ok!"
            else:
                return " : not able to access '" + url + "'"
        except:
            return " : not able to access '" + url + "'"



    #######################################
    # Plotting functions
    #######################################

    @staticmethod
    def plot(xv, yv, title="", xname="x", yname="y", labels=""):

        fig = px.line(width=800, height=400)

        # Error handling

        if not Utils.is_a_vector(xv):
            raise Exception("The parameter 'xv' should be a vector")

        m = max(np.shape(xv))

        if not Utils.is_a_matrix(yv, None, m):
            raise Exception("The parameter 'yv' should be a matrix with " + str(m) + " columns.")

        n = 1 if len(np.shape(yv)) == 1 else np.shape(yv)[0]

        list_names = []

        if str(type(labels)) == "<class 'str'>":
            for i in range(n):
                list_names.append(labels + "_" + str(i + 1))
        else:
            if str(type(labels)) == "<class 'list'>" and len(labels) == n:
                for i in range(n):
                    if str(type(labels[i])) == "<class 'str'>":
                        list_names.append(labels[i])
                    else:
                        raise Exception(
                            "Optional parameter 'labels' must be either a string or a list of " + str(n) + " strings")
            else:
                raise Exception(
                    "Optional parameter 'labels' must be either a string or a list of " + str(n) + " strings")

        # end error handling
        if n > 1:
            for i in range(n):
                fig.add_scatter(x=xv, y=yv[i], mode="lines", name=list_names[i])
        else:
            fig.add_scatter(x=xv, y=yv, mode="lines", name=list_names[0])

        fig.update_xaxes(title_text=xname)
        fig.update_yaxes(title_text=yname)
        fig.show()

        return fig

    #######################################
    # Distance computation functions
    #######################################

    @staticmethod
    def _fun_J(u):
        return Utils._CONSTJA / ((Utils._CONSTJA - 1) * sqrt(Utils._PI * u * u) + sqrt(
            Utils._PI * u * u + Utils._CONSTJA * Utils._CONSTJA))

    @staticmethod
    def fun_Int(v, h, L):
        v = abs(v)
        if v <= L:
            a1 = exp(-(L - v) * (L - v) / (2 * h * h)) * Utils._fun_J((v - L) / (Utils._SQRT2 * h))
            a2 = exp(-(L + v) * (L + v) / (2 * h * h)) * Utils._fun_J((v + L) / (Utils._SQRT2 * h))
            return -h * h * log(Utils._SQRTHALFPI * (h / (2 * L)) * (2 - a1 - a2))
        else:
            a1 = Utils._fun_J((v - L) / (Utils._SQRT2 * h))
            a2 = exp(-2 * L * v / (h * h)) * Utils._fun_J((v + L) / (Utils._SQRT2 * h))
            return 0.5 * (v - L) * (v - L) - h * h * log(Utils._SQRTHALFPI * (h / (2 * L)) * (a1 - a2))

    @staticmethod
    def _fun_I0hat(u):
        return pow(1 + 0.25 * u * u, -0.25) * (1 + Utils._CONSTI0HAT1 * u * u) / (1 + Utils._CONSTI0HAT2 * u * u)

    @staticmethod
    def _fun_f(nu, rho):
        a1 = exp(-0.5 * (rho - nu) * (rho - nu))
        a2 = exp(-0.5 * (rho + nu) * (rho + nu))
        return rho * (a1 + a2) * Utils._fun_I0hat(rho * nu)

    @staticmethod
    def _fun_fhat(nu, rho, rhobar):
        a1 = exp(-0.5 * (rho - nu) * (rho - nu) + 0.5 * (rhobar - nu) * (rhobar - nu))
        a2 = exp(-0.5 * (rho + nu) * (rho + nu) + 0.5 * (rhobar - nu) * (rhobar - nu))
        return rho * (a1 + a2) * Utils._fun_I0hat(rho * nu)

    @staticmethod
    def fun_Cir(v, h, r):

        v = abs(v)
        N = 7
        node = [0.94910, -0.74153, -0.40584, 0, 0.40584, 0.74153, 0.94910]
        weight = [0.12948, 0.27970, 0.38183, 0.4179, 0.38183, 0.27970, 0.12948]

        if v <= r:
            f_low = max(0, sqrt((v / h) * (v / h) + 1) - 3)
            f_up = min(r / h, sqrt((v / h) * (v / h) + 1) + 3)
            delta = 0.5 * (f_up - f_low)
            y = 0
            for i in range(N):
                y = y + weight[i] * Utils._fun_f(v / h, f_low + delta * (node[i] + 1))

            y = delta * y
            return -h * h * log(y * (h / r) * (h / r))
        else:
            f_low = 0
            f_up = r / h
            delta = 0.5 * (f_up - f_low)
            rhobar = f_low + delta * (node[N - 1] + 1)
            y = 0
            for i in range(N):
                y = y + weight[i] * Utils._fun_fhat(v / h, f_low + delta * (node[i] + 1), rhobar)

            y = delta * y
            return 0.5 * (v - h * rhobar) * (v - h * rhobar) - h * h * log(y * (h / r) * (h / r))

    @staticmethod
    def fun_Sph(v, h, r):

        v = abs(v)
        c = 3 * (h * h) / (2 * r * r * r)

        if v <= r:
            if v == 0:
                return -h * h * log(
                    c * (-2 * r * exp(-(r * r) / (2 * h * h)) + 2 * r * exp(-Utils.fun_Int(0, h, r) / (h * h))))
            else:
                a1 = exp(-((r + v) * (r + v) / (2 * h * h)))
                a2 = exp(-((r - v) * (r - v) / (2 * h * h)))
                return -h * h * log(c * (h * h * (a1 - a2) / v + 2 * r * exp(-Utils.fun_Int(v, h, r) / (h * h))))

        else:
            a1 = exp(-(2 * r * v / (h * h)))
            a2 = 1
            return 0.5 * (v - r) * (v - r) - h * h * log(
                c * (h * h * (a1 - a2) / v + 2 * r * exp((0.5 * (v - r) * (v - r) - Utils.fun_Int(v, h, r)) / (h * h))))

    @staticmethod
    def compute_hg_dist(obj_a, obj_b, h=0.001, g=0.001, p_a_init=None, tol=0.001, no_iter_max=20):

        p_a = np.array(p_a_init)
        if p_a is None:
            p_a = np.random.uniform(-3, 3, size=(3,))

        # Error handling

        if not Utils.is_a_simple_object(obj_a):
            raise Exception("The parameter 'obj_a' must be either a box, a ball or a cylinder")

        if not Utils.is_a_simple_object(obj_b):
            raise Exception("The parameter 'obj_b' must be either a box, a ball or a cylinder")

        if not Utils.is_a_number(h) or h <= 0:
            raise Exception("The optional parameter 'h' must be a nonnegative number")

        if not Utils.is_a_number(g) or g <= 0:
            raise Exception("The optional parameter 'g' must be a nonnegative number")

        if not Utils.is_a_vector(p_a_init, 3):
            raise Exception("The optional parameter 'p_a_init' must be a 3D vector")

        if not Utils.is_a_number(tol) or tol <= 0:
            raise Exception("The optional parameter 'tol' must be a nonnegative number")

        if not Utils.is_a_natural_number(no_iter_max):
            raise Exception("The optional parameter 'no_iter_max' must be a nonnegative integer")

        # end error handling

        converged = False
        i = 0

        while (not converged) and i < no_iter_max:
            p_a_ant = p_a
            p_b, dp_a = obj_b.h_projection(p_a, g)
            p_a, dp_b = obj_a.h_projection(p_b, h)
            converged = np.linalg.norm(p_a - p_a_ant) < tol
            i += 1

        dist = np.linalg.norm(p_a - p_b)
        hg_dist = sqrt(2 * (dp_a + dp_b - 0.5 * dist * dist))

        return p_a, p_b, hg_dist
