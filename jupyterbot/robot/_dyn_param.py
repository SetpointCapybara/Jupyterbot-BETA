from utils import *



# Dynamic parameter All
def _dyn_param(self, qdot, q=None):
    if q is None:
        q = self.q

    # Error handling
    if not Utils.is_a_vector(qdot, n):
        raise Exception("The parameter 'qdot' should be a " + str(n) + " dimensional vector")

    if not Utils.is_a_vector(q, n):
        raise Exception("The optional parameter 'q' should be a " + str(n) + " dimensional vector")
    # end error handling

    jac_com, htm_com = self.jac_geo(q, 'com')

    M = self._dynM(q, jac_com, htm_com)
    C = self._dynC(qdot, q, jac_com, FK)
    G = self._dynG(q, jac_com, FK)

    return M, C, G
