'''
MATLAB version 1.0 written by Kerk Phillips, April 2014

PYTHON version adapted by Yulong Li, November 2015 
'''
from __future__ import division
from numpy import tile, array, zeros, log, exp
from LinApp_Sim import LinApp_Sim

def LinApp_SSL(X0, Z, XYbar, logX, PP, QQ, RR, SS):
    '''
    Generates a history of X & Y variables by linearizing the policy function
    about the steady state as in Uhlig's toolkit.
    
    Parameters
    -----------    
    X0: array, dtype=float
        nx vector of X(1) starting values values
    
    Z: 2D-array, dtype=float
        nobs-by-nz matrix of Z values
    
    XYbar: array, dtype=float
        nx+ny vector of X and Y steady state values
    
    logX: binary, dtype=int
        an indicator that determines if the X & Y variables are
        log-linearized (true) or simply linearized (false).  Z variables
        are always simply linearized.
    
    PP: 2D-array, dtype=float
        nx-by-nx matrix of X(t-1) on X(t) coefficients
    
    QQ: 2D-array, dtype=float
        nx-by-nz  matrix of Z(t) on X(t) coefficients
    
    Y0: array, dtype=float
        ny vector of Y(1) starting values values.
    
    RR: 2D-array, dtype=float
        ny-by-nx  matrix of X(t-1) on Y(t) coefficients
    
    SS: 2D-array, dtype=float
        ny-by-nz  matrix of Z(t) on Y(t) coefficients
    
    Returns
    --------
    X: 2D-array, dtype=float
        nobs-by-nx matrix containing the value of the endogenous
        state variables
    
    Y: 2D-array, dtype=float
        nobs-by-ny matrix vector containing the value of the endogenous
        non-state variables
    '''
    # Formating
    X0 = array(X0)

    # get values for nx, ny, nz and nobs
    nobs,nz = Z.shape
    nx = X0.shape[0]
    nxy = XYbar.shape[0]
    ny = nxy - nx

    # get Xbar and Ybar
    Xbar = XYbar[0:nx]
    Ybar = XYbar[nx:nx+ny]

    # Generate a history of X's and Y's
    X = zeros((nobs,nx))
    Y = zeros((nobs,ny))
    Xtil = zeros((nobs,nx))
    Ytil = zeros((nobs,ny))
    
    # set starting values
    if logX:
        Xtil[0,:] = log(X0/Xbar)
    else:
        Xtil[0,:] = X0 - Xbar

    # simulate
    for t in range(0, nobs-1):
        Xtemp, Ytemp = LinApp_Sim(Xtil[t,:], Z[t,:], PP, QQ, RR, SS)
        Xtil[t+1,:] = Xtemp
        if ny>0:
            Ytil[t,:] = Ytemp
        
    # Convert to levels
    if logX:
        X = tile(Xbar,(nobs,1))*exp(Xtil)
        if ny> 0:
            Y = tile(Ybar,(nobs,1))*exp(Ytil)
    else:
        X = tile(Xbar,(nobs,1))+Xtil
        if ny>0:
            Y = tile(Ybar,(nobs,1))+Ytil

    return X, Y #Note if ny=0, Y is a nobs by 0 empty matrix 
    