import numpy as np
from utils import *

_INVHALFPI = 0.63660


def _vector_field(curve, alpha, const_vel):
    # Error handling
    if not Utils.is_a_matrix(curve):
        raise Exception("The optional parameter 'curve' should be a matrix of float numbers")

    vector_size = len(curve[0])

    if not Utils.is_a_number(alpha) or alpha <= 0:
        raise Exception("The optional parameter 'alpha' should be a positive float")

    if not Utils.is_a_number(const_vel):
        raise Exception("The optional parameter 'const_vel' should be a float")

    # end error handling

    return lambda p: _vector_field_vel(p, curve, alpha, const_vel, vector_size)


def _vector_field_vel(p, curve, alpha, const_vel, vector_size):
    # Error handling
    if not Utils.is_a_vector(p, vector_size):
        raise Exception("The vector field argument should be a " + str(vector_size) + " dimensional vector")
    # end error handling

    N, T, D = _compute_NTD(curve, p)
    G = _INVHALFPI * atan(alpha * D)
    H = sqrt(1 - G * G)
    abs_const_vel = abs(const_vel)
    sgn = const_vel / (abs_const_vel + 0.00001)

    return (abs_const_vel * (G * N + sgn * H * T)).reshape((vector_size,))


def _compute_NTD(curve, p):
    min_dist = float('inf')
    ind_min = -1
    for i in range(len(curve)):
        dist_temp = np.linalg.norm(np.array(p) - np.array(curve[i]))
        if dist_temp < min_dist:
            min_dist = dist_temp
            ind_min = i

    N = np.array(curve[ind_min]) - np.array(p)
    N = N / (np.linalg.norm(N) + 0.0001)

    if ind_min == len(curve) - 1:
        T = np.array(curve[1]) - np.array(curve[ind_min])
    else:
        T = np.array(curve[ind_min + 1]) - np.array(curve[ind_min])

    T = T / (np.linalg.norm(T) + 0.0001)

    return N, T, min_dist
