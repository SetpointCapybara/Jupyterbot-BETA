from utils import *
from meshmaterial import *
import numpy as np


class Model3D:

    #######################################
    # Attributes
    #######################################

    @property
    def url(self):
        """The obj url."""
        return self._url

    @property
    def scale(self):
        """The object scale."""
        return self._scale

    @property
    def htm(self):
        """Model pose written in the Denavit-Hartenberg frame of the link.
        A 4x4 homogeneous transformation matrix."""
        return self._htm

    @property
    def mesh_material(self):
        """The model mesh material"""
        return self._mesh_material

    #######################################
    # Constructor
    #######################################

    def __init__(self, url="", scale=1, htm = np.identity(4), mesh_material=None):

        # Error handling

        error = Utils.is_url_available(url, ['obj'])
        if not (error == "ok!"):
            raise Exception("The parameter 'url' " + error)

        if not Utils.is_a_matrix(htm, 4, 4):
            raise Exception("The parameter 'htm' should be a 4x4 homogeneous transformation matrix")

        if not Utils.is_a_number(scale) or scale < 0:
            raise Exception("The parameter 'scale' should be a float")

        if not (Utils.get_jupyterbot_type(mesh_material) == "jupyterbot.MeshMaterial" or (mesh_material is None)):
            raise Exception("The parameter 'mesh_material' should be a 'MeshMaterial' object or none")

        # end error handling

        self._url = url
        self._scale = scale
        self._htm = htm

        if mesh_material is None:
            self._mesh_material = MeshMaterial()
        else:
            self._mesh_material = mesh_material

    #######################################
    # Std. Print
    #######################################

    def __repr__(self):

        string = "3D object model with url: '" + self.url + "': \n\n"

        return string

    #######################################
    # Methods
    #######################################

    def gen_code(self, name):
        """Generate code for injection."""

        string = self.mesh_material.gen_code(name)
        string += "const htm_" + name + " = new Matrix4(); \n"
        string += "htm_" + name + ".set(" \
        + str(self.htm[0][0]) + "," + str(self.htm[0][1]) + "," + str(self.htm[0][2]) + "," + str(self.htm[0][3]) + "," \
        + str(self.htm[1][0]) + "," + str(self.htm[1][1]) + "," + str(self.htm[1][2]) + "," + str(self.htm[1][3]) + "," \
        + str(self.htm[2][0]) + "," + str(self.htm[2][1]) + "," + str(self.htm[2][2]) + "," + str(self.htm[2][3]) + "," \
        + str(self.htm[3][0]) + "," + str(self.htm[3][1]) + "," + str(self.htm[3][2]) + "," + str(self.htm[3][3]) + ");\n"

        string += "const object3d_" + name + "={\n"
        string += "url: '" + self.url + "',\n"
        string += "scale: " + str(self.scale) + ",\n"
        string += "matrix: htm_" + name + ",\n"
        string += "mesh_material: material_" + name + "};\n\n"

        return string
