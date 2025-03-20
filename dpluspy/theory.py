
import numpy as np


def U_operator(u):

    U = u * np.array([
        [0, 0, 0],
        [2, 0, 0],
        [0, 4, 0]
    ])
    return U


def D_operator(Ne):

    D = 1 / (2 * Ne) * np.array([
        [0, 0, 0],
        [0, -1, 0],
        [0, 0, -1]
    ])
    return D


def R_operator(r):

    R = r * np.array([
        [0, 0, 0],
        [0, 0, 0],
        [0, 0, -1]
    ])
    return R


def get_operator(Ne, u, r):

    U = U_operator(u)
    D = D_operator(Ne)
    R = R_operator(r)
    return U + D + R




