from utils import *
import numpy as np


# Dynamic parameter Gravity
def _dynG(self, q=None):
    if q is None:
        q = self.q

    # Error handling
    if not Utils.is_a_vector(qdot, n):
        raise Exception("The parameter 'qdot' should be a " + str(n) + " dimensional vector")
    # end error handling

    jac_com, htm_com = self.jac_geo(q, 'com')
    return self._dynG_aux(qdot, q, jac_com, htm_com)


def _dynG_aux(self, q, jac_com, htm_com):
    n = len(self._links)
    G = np.zeros((n,))
    gravity_acc = 9.8;
    for i in range(n):
        G = G + gravity_acc * self._links[i].mass * np.transpose(jac_com[i, 3, :]),

    return np.reshape(G, (n,))
