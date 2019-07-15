'''
Version 1.0, written by Kerk Phillips, April 2014

Adapted by Yulong Li, November 2015 
'''
from __future__ import division
import numpy as np

def LinApp_Sim(Xm, Z, PP, QQ, RR, SS):
    '''
    Uses the coefficients from a linear approximation to generate data for
    next period given today's state. The set of endogenous state variables 
    known today is Xm and the set of exogenous state variables is Z.
    This program generates X.  

    The input and output values are in deviation from the linearization point 
    (almost always the steady  state, but not necessarily so).  This means 
    you will need to add back the steady state or other values after you have 
    called this function.  How you do this depends on whether you used 
    log-linearization or simple linearization in deriving the values of the 
    input coefficients.

    Parameters
    -----------
    Xm: array, dtype=float
        nx vector of X(t-1) values

    Z: array, dtype=float
        nz vector of Z(t) values

    PP: 2D-array, dtype=float
        nx-by-nx  matrix of X(t-1) on X(t) coefficients

    QQ: 2D-array, dtype=float
        nx-by-nz  matrix of Z(t) on X(t) coefficients

    RR: 2D-array, dtype=float
        ny-by-nx  matrix of X(t-1) on Y(t) coefficients

    SS: 2D-array, dtype=float
        ny-by-nz  matrix of Z(t) on Y(t) coefficients

    Returns
    --------
    X: array, dtype=float
        nx vector containing the value of the endogenous
        state variables for next period
    
    Y: array, dtype=float
        ny vector containing the value of the endogenous
        non-state variables for the current period
    '''
    # Find the number of each kind of state variable
    # Using Xm find nx
    if len(Xm.shape)!=1:
        print('Xm must be a one-dimensional array')
    else:
        nx = Xm.shape[0]
        # Using RR find ny
        ny = RR.shape[0]

        # Using Z find nz
        if len(Z.shape)!=1:
            print('Z must be a one-dimensional array')
        else:
            nz = Z.shape[0]

            # Check conformity of input coefficient matrices
            d1,d2 = PP.shape
            if d1 != nx or d2 != nx:
                print('Dimensions of PP incorrect')

            d1,d2 = QQ.shape
            if d1 != nx or d2 != nz:
                print('dimensions of QQ incorrect')

            # Generate data for next period, one equation at a time
            X = PP.dot(Xm) + QQ.dot(Z)

            if ny>0:
                Y = RR.dot(Xm) + SS.dot(Z)
            else:
                Y = []
                
    return np.array(X), np.array(Y)
