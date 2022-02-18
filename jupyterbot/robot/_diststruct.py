import numpy as np


class _DistStruct:

    #######################################
    # Attributes
    #######################################

    @property
    def id_obj(self):
        """Return the memory address for the associated object."""
        return self._id_obj

    @property
    def id_robot(self):
        """Return the memory address for the associated robot."""
        return self._id_robot

    @property
    def jac_dist_mat(self):
        """
		Return the matrix in which each row we have the distance Jacobian (gradient) for each robot link.
		"""
        return self._jac_dist_mat

    @property
    def dist_vect(self):
        """
		Return the vector in which each row we have the distance for each robot link.
		"""
        return np.array(self._dist_vect).reshape((self.no_items,))

    @property
    def no_items(self):
        """Return the number of items."""
        return self._no_items

    def __getitem__(self, key):
        return self._list_dict[key]

    #######################################
    # Constructor
    #######################################

    def __init__(self, obj, robot):

        self._id_obj = id(obj)
        self._id_robot = id(robot)
        self._obj_name = obj.name
        self._robot_name = robot.name
        self._no_items = 0

        n = len(robot.links)
        self._list_dict = [];
        self._jac_dist_mat = np.zeros((0, n))
        self._dist_vect = []

    #######################################
    # Std. Print
    #######################################

    def __repr__(self):

        return "Distance struct between robot '" + self._robot_name + "' and object '" \
               + self._obj_name + "', with " + str(self.no_items) + " items"

    #######################################
    # Methods
    #######################################

    def _append(self, i, j, d, pObj, pObjCol, Jd):
        self._list_dict.append({
            "linkNumber": i,
            "linkColObjNumber": j,
            "hgDistance": d,
            "pObj": pObj,
            "pObjCol": pObjCol,
            "jacDist": Jd
        })

        self._jac_dist_mat = np.vstack((self._jac_dist_mat, Jd))
        self._dist_vect.append(d)
        self._no_items += 1

    def get_item(self, i, j):
        for d in self._list_dict:
            if i == d["linkNumber"] and j == d["linkColObjNumber"]:
                return d

        raise Exception("Item not found!")

    def get_closest_item(self):
        dmin = float('inf')
        imin = -1
        for i in range(self._no_items):
            if self[i]["hgDistance"] < dmin:
                dmin = self[i]["hgDistance"]
                imin = i

        return self[imin]
