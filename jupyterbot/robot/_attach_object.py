from utils import *


def _attach_object(self, obj):
    # Error handling
    if not (Utils.is_a_simple_object(obj)):
        raise Exception("The parameter 'obj' must be either a box, a ball or a cylinder")
    # end error handling

    self._attached_objects.append([obj, Utils.inv_htm(self.fkm()) @ obj.htm])
