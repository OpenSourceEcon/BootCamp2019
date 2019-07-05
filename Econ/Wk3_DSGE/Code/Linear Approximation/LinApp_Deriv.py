"""
Created November, 2015
Revised January, 2017

Author:  Kerk Phillips

Based on code by Spencer Lyon
MatLab code by Kerk P. (2013) was referenced in creating this file.
"""
from __future__ import division
import numpy as np

def LinApp_Deriv(funcname, param, theta0, nx, ny, nz, logX):
    """
    This function computes the matricies AA-MM in the log-linearization of
    the equations in the function 'func'.

    Parameters
    ----------
    func: function
        The function that generates a vector from the dynamic equations that are
        to be linearized. This must be written to evaluate to zero in the
        steady state. Often these equations are the Euler equations defining
        the model

    theta0: array, dtype=float
        A vector of steady state values for state parameters. Place the values
        of X in the front, then the Y values, followed by the Z's.

    nx: number, dtype=int
        The number of elements in X

    ny: number, dtype=nt
        The number of elements in Y

    nz: number, dtype=int
        The number of elements in Z

    logX: binary, dtype=int
        true if log-linearizing the X & Y variables, false for simple linearization

    Returns
    -------
    AA - MM : 2D-array, dtype=float:
        The equaitons described by Uhlig in the log-linearization.
    """
    # calculate of derivative matrix
    length = 3 * nx + 2 * (ny + nz)
    height = nx + ny
    
    # set vaue for epsilon
    # eps = 2.2E-16  # machine epsilon for double-precision
    eps = 10E-6

    # Constant term, T0, should be very close to zero if linearizing about SS
    T0 = funcname(theta0, param)  
    
    # set up plus and minus deviations matrices, disturb each input one-by-one
    devplus = np.tile(theta0.reshape(1, theta0.size), (length, 1))
    devminus = np.tile(theta0.reshape(1, theta0.size), (length, 1))
    for i in range(length):
        devplus[i, i] += eps
        devminus[i, i] -= eps
        
    # initialize first derivative  matrix
    dGammaMat = np.zeros((height,length))
    
    # find first derivatives
    for i in range(0,length):
        if i < 3 * nx + 2 * ny:
            if logX:
                dGammaMat[:,i] = \
                (theta0[i]*(funcname(devplus[i, :], param)-T0)/(1.0+T0) \
                - theta0[i]*(funcname(devminus[i, :], param)-T0)/(1.0+T0)) \
                / (2.0 * eps)
            else:
                dGammaMat[:,i] = \
                (funcname(devplus[i, :], param) \
                -funcname(devminus[i, :], param)) / (2.0 * eps)
        else:
            dGammaMat[:,i] = \
            (funcname(devplus[i, :], param) \
            -funcname(devminus[i, :], param)) / (2.0 * eps)

    # partition into parts as labeled by Uhlig
    AA = np.array(dGammaMat[0:ny, nx:2 * nx])
    BB = np.array(dGammaMat[0:ny, 2 * nx:3 * nx])
    CC = np.array(dGammaMat[0:ny, 3 * nx + ny:3 * nx + 2 * ny])
    DD = np.array(dGammaMat[0:ny, 3 * nx + 2 * ny + nz:length])
    FF = np.array(dGammaMat[ny:ny + nx, 0:nx])
    GG = np.array(dGammaMat[ny:ny + nx, nx:2 * nx])
    HH = np.array(dGammaMat[ny:ny + nx, 2 * nx:3 * nx])
    JJ = np.array(dGammaMat[ny:ny + nx, 3 * nx:3 * nx + ny])
    KK = np.array(dGammaMat[ny:ny + nx, 3 * nx + ny:3 * nx + 2 * ny])
    LL = np.array(dGammaMat[ny:ny + nx, 3 * nx + 2 * ny:3 * nx + 2 * ny + nz])
    MM = np.array(dGammaMat[ny:ny + nx, 3 * nx + 2 * ny + nz:length])

    return [AA, BB, CC, DD, FF, GG, HH, JJ, KK, LL, MM]
