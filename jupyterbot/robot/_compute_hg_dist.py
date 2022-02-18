from utils import *
import numpy as np
from ._diststruct import _DistStruct


# Compute the distance from each link to an object, for the current configuration
# of the robot
def _compute_hg_dist(self, obj, h, g, q=None, htm=None, old_dist_struct=None, tol=0.0005, no_iter_max=20):
    n = len(self.links)

    if q is None:
        q = self.q

    if htm is None:
        htm = self.htm

    # Error handling
    if not Utils.is_a_vector(q, n):
        raise Exception("The optional parameter 'q' should be a " + str(n) + " dimensional vector")

    if not Utils.is_a_matrix(htm, 4, 4):
        raise Exception("The optional parameter 'htm' should be a 4x4 homogeneous transformation matrix")

    if not (old_dist_struct is None):
        try:
            id_robot = old_dist_struct.id_robot
            id_obj = old_dist_struct.id_obj
            if not (id_obj == id(obj) and id_robot == id(self)):
                Exception("The optional parameter 'old_dist_struct' is a 'DistStruct' object, but it " \
                          "must have to be relative to the SAME robot object and SAME external object, and " \
                          "this is not the case")
        except:
            raise Exception("The optional parameter 'old_dist_struct' must be a 'DistStruct' object")

    # end error handling

    dist_struct = _DistStruct(obj, self)

    jac_dh, mth_dh = self.jac_geo(q, "dh", htm)

    col_object_copy = []

    # Update all collision objects of all links
    for i in range(n):
        col_object_copy.append([])
        for j in range(len(self.links[i].col_objects)):
            temp_copy = self.links[i].col_objects[j][0].copy()
            htmd = self.links[i].col_objects[j][1]
            temp_copy.set_ani_frame(mth_dh[i, :, :] @ htmd)
            col_object_copy[i].append(temp_copy)

    # Compute the distance structure
    for i in range(n):
        for j in range(len(self.links[i].col_objects)):

            if old_dist_struct is None:
                p_obj_0 = np.random.uniform(-100, 100, size=(3,))
            else:
                p_obj_0 = old_dist_struct.get_item(i, j)["pObj"]

            p_obj, p_obj_col, d = Utils.compute_hg_dist(obj, col_object_copy[i][j], h, g, p_obj_0, tol, no_iter_max)

            jac_obj_col = jac_dh[i, 0:3, :] - Utils.S(p_obj_col - mth_dh[i, 0:3, 3]) @ jac_dh[i, 3:6, :]
            jac_dist = (np.transpose(p_obj_col - p_obj) @ jac_obj_col) / d

            dist_struct._append(i, j, d, p_obj, p_obj_col, jac_dist)

    return dist_struct
