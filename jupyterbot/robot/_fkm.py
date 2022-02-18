from utils import *
import numpy as np


def _fkm(self, q=None, axis='eef', htm=None):
    if q is None:
        q = self.q
    if htm is None:
        htm = np.array(self.htm)

    n = len(self._links)

    # Error handling
    if not Utils.is_a_vector(q, n):
        raise Exception("The optional parameter 'q' should be a " + str(n) + " dimensional vector")

    if not (axis == "eef" or axis == "dh" or axis == "com"):
        raise Exception("The optional parameter 'axis' should be one of the following strings:\n" \
                        "'eef': End-effector \n" \
                        "'dh': All " + str(n) + " axis of Denavit-Hartenberg\n" \
                                                "'com': All " + str(
            n) + " axis centered at the center of mass of the objects")

    if not Utils.is_a_matrix(htm, 4, 4):
        raise Exception("The optional parameter 'htm' should be a 4x4 homogeneous transformation matrix")
    # end error handling

    #q = q + self._dq0

    htm_dh = np.zeros((n, 4, 4))

    for i in range(n):
        if i == 0:
            htm_dh[i, :, :] = htm
        else:
            htm_dh[i, :, :] = htm_dh[i - 1, :, :]

        if self._links[i].joint_type == 0:
            htm_dh[i, :, :] = htm_dh[i, :, :] @ Utils.rotz(q[i])
        else:
            htm_dh[i, :, :] = htm_dh[i, :, :] @ Utils.rotz(self._links[i].theta)

        if self._links[i].joint_type == 1:
            htm_dh[i, :, :] = htm_dh[i, :, :] @ Utils.trn([0, 0, q[i]])
        else:
            htm_dh[i, :, :] = htm_dh[i, :, :] @ Utils.trn([0, 0, self._links[i].d])

        htm_dh[i, :, :] = htm_dh[i, :, :] @ Utils.rotx(self._links[i].alpha)
        htm_dh[i, :, :] = htm_dh[i, :, :] @ Utils.trn([self._links[i].a, 0, 0])

    if axis == 'com':
        for i in range(n):
            htm_dh[i, 3, 0:3] = htm_dh[i, 3, 0:3] + htm_dh[i, 0:3, 0:3] @ self._links[i].center_shift

    if axis == 'eef':
        htm_dh = htm_dh[-1, :, :]

    return htm_dh
