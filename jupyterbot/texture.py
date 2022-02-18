from utils import *


class Texture:

    #######################################
    # Attributes
    #######################################

    @property
    def url(self):
        """The address of the texture."""
        return self._url

    @property
    def mapping(self):
        """A list of all object names."""
        return self._mapping

    @property
    def wrap_s(self):
        """A list of all object names."""
        return self._wrap_s

    @property
    def wrap_t(self):
        """A list of all object names."""
        return self._wrap_t

    @property
    def mag_filter(self):
        """A list of all object names."""
        return self._mag_filter

    @property
    def min_filter(self):
        """A list of all object names."""
        return self._min_filter

    @property
    def offset(self):
        """A list of all object names."""
        return self._offset

    @property
    def repeat(self):
        """A list of all object names."""
        return self._repeat

    @property
    def rotation(self):
        """A list of all object names."""
        return self._rotation

    @property
    def center(self):
        """A list of all object names."""
        return self._center

    #######################################
    # Constructor
    #######################################

    def __init__(self, url, mapping="UVMapping", wrap_s="ClampToEdgeWrapping", wrap_t="ClampToEdgeWrapping",
                 mag_filter="LinearFilter", min_filter="LinearMipmapLinearFilter", offset=[0,0], repeat=[1,1],
                 rotation=0, center=[0,0]):

        image_types = ["png", "bmp", "jpg", "jpeg"]
        mapping_list = ['UVMapping', 'CubeReflectionMapping', 'CubeRefractionMapping',
                        'EquirectangularReflectionMapping', 'EquirectangularRefractionMapping',
                        'CubeUVReflectionMapping', 'CubeUVRefractionMapping']
        wrapping_mode_list = ['RepeatWrapping', 'ClampToEdgeWrapping', 'MirroredRepeatWrapping']
        mag_filter_list = ['NearestFilter', 'LinearFilter']
        min_filter_list = ['NearestFilter', 'NearestMipmapNearestFilter', 'NearestMipmapLinearFilter', 'LinearFilter',
                           'LinearMipmapNearestFilter', 'LinearMipmapLinearFilter']

        error = Utils.is_url_available(url, image_types)
        if not(error == "ok!"):
            raise Exception("The parameter 'url' "+error)

        if not (str(type(mapping)) == "<class 'str'>" and (mapping in mapping_list)):
            raise Exception("The parameter 'mapping' should be one of the following strings: "+str(mapping_list))

        if not (str(type(wrap_s)) == "<class 'str'>" and (wrap_s in wrapping_mode_list)):
            raise Exception("The parameter 'wrap_s' should be one of the following strings: "+str(wrapping_mode_list))

        if not (str(type(wrap_t)) == "<class 'str'>" and (wrap_t in wrapping_mode_list)):
            raise Exception("The parameter 'wrap_t' should be one of the following strings: "+str(wrapping_mode_list))

        if not (str(type(mag_filter)) == "<class 'str'>" and (mag_filter in mag_filter_list)):
            raise Exception("The parameter 'mag_filter' should be one of the following strings: "+str(mag_filter_list))

        if not (str(type(min_filter)) == "<class 'str'>" and (min_filter in min_filter_list)):
            raise Exception("The parameter 'min_filter_list' should be one of the following strings: "+str(min_filter_list))

        if not Utils.is_a_vector(offset, 2):
            raise Exception("The parameter 'offset' should be 2D vector")

        if not Utils.is_a_vector(repeat, 2):
            raise Exception("The parameter 'repeat' should be 2D vector")

        if not Utils.is_a_number(rotation):
            raise Exception("The parameter 'rotation' should be a float")

        if not Utils.is_a_vector(center, 2):
            raise Exception("The parameter 'center' should be 2D vector")

        self._url = url
        self._mapping = mapping
        self._wrap_s = wrap_s
        self._wrap_t = wrap_t
        self._mag_filter = mag_filter
        self._min_filter = min_filter
        self._offset = offset
        self._repeat = repeat
        self._rotation = rotation
        self._center = center

    #######################################
    # Std. Print
    #######################################

    def __repr__(self):

        string = "Texture: '"+self.url+"'"

        return string

    #######################################
    # Methods
    #######################################

    def gen_code(self, name):

        string = "const texture_"+name+" = new TextureLoader().load('"+self.url+"');\n"
        string += "texture_" + name + ".wrapS = "+str(self.wrap_s)+";\n"
        string += "texture_" + name + ".wrapT = " + str(self.wrap_t) + ";\n"
        string += "texture_" + name + ".magFilter = " + str(self.mag_filter) + ";\n"
        string += "texture_" + name + ".minFilter = " + str(self.min_filter) + ";\n"
        string += "texture_" + name + ".rotation = " + str(self.rotation) + ";\n"
        string += "texture_" + name + ".offset.set("+str(self.offset[0])+","+str(self.offset[1])+");\n"
        string += "texture_" + name + ".repeat.set("+str(self.repeat[0])+","+str(self.repeat[1])+");\n"
        string += "texture_" + name + ".center.set("+str(self.center[0])+","+str(self.center[1])+");\n\n"

        return string