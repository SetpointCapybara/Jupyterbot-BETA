from utils import *
import numpy as np


# Dynamic parameter Inertia Matrix
def _dynM(self, q=None):
    if q is None:
        q = self.q

    n = len(self._links)

    # Error handling
    if not Utils.is_a_vector(q, n):
        raise Exception("The optional parameter, 'q', should be a " + str(n) + " dimensional vector")
    # end error handling

    jac_com, htm_com = self.jac_geo(q, 'com')

    return self._dynM_aux(q, jac_com, htm_com)


def _dynM_aux(self, q, jac_com, htm_com):
    n = len(self._links)
    M = np.zeros((n, n))

    for i in range(n):
        M = M + self._links[i].mass * np.transpose(jac_com[i, 0:3, :]) @ jac_com[i, 0:3, :]
        M = M + np.transpose(jac_com[i, 3:6, :]) @ htm_com[i, 0:3, 0:3] @ \
            self._links[i].inertia_matrix @ \
            np.transpose(htm_com[i, 0:3, 0:3]) @ jac_com[i, 3:6, :]

    return M
