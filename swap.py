"""
Created on Mon Mar 31 01:00:02 2014

Author: Divyanshu Vats (vats.div@gmail.com)
"""

import numpy as np
from numpy.linalg import inv


def SWAP(X, y, s_ini):
    """
    Function SWAP:
    --------------
    Returns an estimate of the locations of the non-zero
    entries in the vector b s.t. y = X * b + noise

    Parameters
    ----------
    X: design matrix
    y: observations
    s_ini: initial estimate of the support

    Returns
    -------
    estimate of the support computing using the SWAP algorithm
    """

    k = np.size(s_ini)
    n, p = X.shape

    while 1:

        Delta = np.zeros((k, p-k))

        s_temp_c = np.delete(np.arange(0, p), s_ini)

        for i in np.arange(0, k):
            # remove the ith index from s_temp
            s_temp = np.delete(s_ini, i)
            PiPerp = np.eye(n) - np.dot(np.dot(X[:, s_temp],
                                               inv(np.dot(X[:, s_temp].T,
                                                   X[:, s_temp]))),
                                        X[:, s_temp].T)
            temp = np.dot(PiPerp, X[:, s_ini[i]])
            ytemp = np.dot(y.T, temp)
            yPi = np.dot(ytemp, ytemp.T) / np.dot(temp.T, temp)

            temp = np.dot(PiPerp, X[:, s_temp_c])
            ddt = np.dot(y.T, temp)
            ddtt = ddt * ddt
            temp_new = np.sum(temp.T * temp.T, axis=1)
            ddtt = ddtt / temp_new

            Delta[i, :] = yPi - ddtt

        i1, i2 = np.unravel_index(np.argmin(Delta), Delta.shape)

        # swap i1 with i2
        print i1, i2
        s_ini[i1] = s_temp_c[i2]

        if np.amin(Delta) > 0:
            break

    return s_ini
