from utils import *
import numpy as np


# (Private function) Used in IK
def _evolve_config(self, q, p_tol, a_tol, htm_target, iter_remain):
    n = len(self.links)
    found = False
    zero_u = False
    iter_end = False

    dt = min(p_tol, 0.005 * a_tol, 0.01)
    K = 2
    eps = 0.001
    tol_u = 0.0001 / dt
    i = 0

    while (not found) and (not zero_u) and (not iter_end):
        r, Jr = self.task_function(htm_target, np.array(q))
        r[3:6] = np.sqrt(r[3:6])

        u = Utils.dp_inv(Jr, eps) @ (-K * r)
        q = q + u * dt

        e_pos = max(abs(r[0:3]))
        e_ori = max([(180 / np.pi) * acos(min(max(1 - num * num, -1), 1)) for num in r[3:6]])
        i += 1

        found = (e_pos < p_tol) and (e_ori < a_tol)
        zero_u = max(abs(u)) < tol_u
        iter_end = i > iter_remain

    return found, i, q


# Inverse kinematics for the end-effector
def _ikm(self, htm_target, q0=None, p_tol=0.005, a_tol=5, no_iter_max=2000):
    n = len(self._links)
    if q0 is None:
        q0 = (2 * np.pi) * np.random.rand(n) - np.pi

        # Error handling
    if not Utils.is_a_matrix(htm_target, 4, 4):
        raise Exception("The parameter 'htm' should be a 4x4 homogeneous transformation matrix")

    if not Utils.is_a_natural_number(no_iter_max):
        raise Exception("The optional parameter 'no_iter_max' should be a nonnegative integer number")

    if not Utils.is_a_vector(q0, n):
        raise Exception("The optional parameter 'q0' should be a " + str(n) + " dimensional vector")

    if (not Utils.is_a_number(p_tol)) or p_tol <= 0:
        raise Exception("The optional parameter 'pol' should be a nonnegative number")

    if (not Utils.is_a_number(a_tol)) or a_tol <= 0:
        raise Exception("The optional parameter 'a_tol' should be a nonnegative number")

    # end error handling

    j = 0
    found = False
    q = q0
    no_iter_remain = no_iter_max

    while not found and no_iter_remain >= 0:
        found, i, q = _evolve_config(self, q, p_tol, a_tol, htm_target, no_iter_remain)
        no_iter_remain -= i
        if not found:
            q = (2 * np.pi) * np.random.rand(n)

    if not found:
        raise Exception("Solution for IK not found. You can try the following: \n" \
                        " Increasing the maximum number of iterations, 'no_iter_max' (currently " + str(no_iter_max) + ")\n" \
                                                                                                                   " Increasing the tolerance for the position, 'p_tol' (currently " + str(
            p_tol) + " meters)\n" \
                    " Increasing the tolerance for the orientation, 'a_tol' (currently " + str(a_tol) + " degrees)")
    else:
        return q
