from IPython.core.display import display, HTML
import re
from utils import *
import sys


class Simulation:
    """A simulation variable."""

    # Import the javascript code as a string
    _STRJAVASCRIPT = "<canvas id=\"scene\" width=\"800\" height=\"600\"></canvas>"
    _STRJAVASCRIPT += "\n <script type=\"module\">"
    for line in open("D:\\PycharmProjects\\pyProj\\jupyterbot\\threejs_sim.js").readlines():
        _STRJAVASCRIPT += line

    _STRJAVASCRIPT += "\n </script>"

    #######################################
    # Attributes
    #######################################

    @property
    def list_of_objects(self):
        """A list of all objects."""
        return self._list_of_objects

    @property
    def list_of_names(self):
        """A list of all object names."""
        return self._list_of_names

    @property
    def ambient_light_intensity(self):
        """A list of all object names."""
        return self._ambient_light_intensity

    #######################################
    # Constructor
    #######################################

    def __init__(self, obj_list=[], ambient_light_intensity=12):

        if not Utils.is_a_number(ambient_light_intensity) or ambient_light_intensity < 0:
            raise Exception("The parameter 'ambient_light_intensity' should be a nonnegative float")

        self._list_of_objects = []
        self._list_of_names = []
        self._ambient_light_intensity = ambient_light_intensity

        if str(type(obj_list)) == "<class 'list'>":
            for obj in obj_list:
                self.add(obj)
        else:
            raise Exception("The parameter 'obj_list' should be a list of objects")




    #######################################
    # Std. Print
    #######################################

    def __repr__(self):

        string = "Simulation: \n\n"
        string += " Variables: \n"
        string += str(self.list_of_names)

        return string

    #######################################
    # Methods
    #######################################

    def gen_code(self, speed_multiplier=1):

        if not Utils.is_a_number(speed_multiplier) or speed_multiplier <= 0:
            raise Exception("The parameter 'speed_multiplier' should be a positive number")

        string = Simulation._STRJAVASCRIPT

        string = re.sub("//USER INPUT GOES HERE", "ambientLight.intensity = "+str(self.ambient_light_intensity)+";\n //USER INPUT GOES HERE",
                            string)

        string = re.sub("//SIMULATION PARAMETERS GO HERE", "const var_multiplier = " + str(speed_multiplier) + "; \n //SIMULATION PARAMETERS GO HERE",
                            string)

        string = re.sub("//SIMULATION PARAMETERS GO HERE", "const delay = 2000; \n //SIMULATION PARAMETERS GO HERE",
                            string)


        for obj in self.list_of_objects:
            string = re.sub("//USER INPUT GOES HERE", obj.gen_code(), string)

        return string

    def run(self, speed_multiplier=1):
        """Run simulation."""

        if not Utils.is_a_number(speed_multiplier) or speed_multiplier <= 0:
            raise Exception("The parameter 'speed_multiplier' should be a positive number")

        display(HTML(self.gen_code()))

    def save(self, address, file_name, speed_multiplier=1):

        if not Utils.is_a_number(speed_multiplier) or speed_multiplier <= 0:
            raise Exception("The parameter 'speed_multiplier' should be a positive number")

        file = open(address+"/"+file_name+".html", "w+")
        file.write(self.gen_code())
        file.close()

    def add(self, obj_sim):
        """
    Add an object to the simulation. It should be an object that
    can be simulated (Utils.is_a_obj_sim(obj) is true).

    Parameters
    ----------
    obj_sim : obj
        The object to be added.
    """

        # Error handling
        if not Utils.is_a_obj_sim(obj_sim):
            raise Exception("The parameter 'obj' should be able to be simulated (Utils.is_a_obj_sim(obj) is true)")

        if obj_sim.name in self.list_of_names:
            raise Exception("The name '" + obj_sim.name + "' is already in the list of symbols")

        # end error handling

        self._list_of_names.append(obj_sim.name)
        self._list_of_objects.append(obj_sim)
