# -*- coding: utf-8 -*-
"""
Version 1.0
Fri Nov 25 2016
@author: Kerk L. Phillips
This program calculates Euler Errors from a DSGE model.  

"""
import numpy as np
from scipy.stats import norm

def EEcalc(Xdata, Zdata, efunc, epars, tfunc, tpars, lfunc, lpars):
    '''
    This function calculates Euler Errors for a DSGE model.
    It takes the following as inputs:
    1)  Xdata: an nobs-by-nx matrix of endogenous state variable data with with 
        nobs indexing the observations in rows and nx indexing the variables in 
        columns.
    2)  Zdata: an nobs-by-nz matrix of exogenous state variable data with with 
        nobs indexing the observations in rows and nz indexing the variables in 
        columns.
    3)  efunc: the name of a function which takes a single observation of data 
        and returns the Euler error for a given realization of future exogenous 
        variables.  This function must take 3 nx-element vectors of endogenous
        state variables and and 2 nz-element vector of exogenous state 
        variables as inputs and output an neq-element vector of Euler equation 
        errors.  The order of input is X(t+2), X(t+1), X(t), Z(t+1), Z(t),
        epars
    4)  epars: a list of parameters passed to efunc.
    5)  tfunc: the name of a transition function which takes a single 
        observation of the state variables and returns next period's value of 
        the endogenous state variables.  This function must take an nx-element 
        vector of endogenous state variables and and nz-element vector of 
        exogenous state variables as inputs and output an nz-element vector of 
        next-period endogenous state variables.  The order of inputs is X(t), 
        Z(t), tpars
    6)  tpars: a list of parameters passed to tfunc.
    7)  lfunc: the name of a law-of-motion function which takes a single 
        observation of the exogenous state variables and returns next period's 
        value of the exogenous state variables.  This function must take an 
        nz-element of exogenous state and an scalar iid shock as inputs 
        and output an nz-element vector of next-period endogenous state 
        variables. The order of inputs is Z(t), Eps(t+1), epars
    8)  lpars: a list of parameters passed to lfunc.
    
    The function returns the following outputs:
    1)  Eerr: an nobs-by-neq matrix of Euler errors with nobs indexing the 
        observations in rows and neq indexing the elements from the function 
        efunc in columns.
        
    Notes:
    Xdata and Zdata must have the same number of rows.
    Neither Xdata nor Zdata may have missing, nan, or complex values. 
    Innovations to the law of motion are drawn from a standard normal 
    distribution.
    Currently this function only works with one innovation shock, i.e. ne=1

    To Do:
    1) Allow for more than one shock process.  May require the use of sparse 
    grids for quadrature.
    2) Use a more efficient quadrature method.  Gaussian?
    '''
    
    # set parameter values
    npts = 10 # number of point for rectangular quadrature
    
    # check sizes of data matrices
    (Xnobs, nx) = Xdata.shape
    (Znobs, nz) = Zdata.shape
    if Xnobs == Znobs:
        nobs = Xnobs
    else:
        print('Data matrices have different numbers of observations')
        nobs = min(Xnobs, Znobs)
        Xdata = Xdata[0:nobs]
        Zdata = Zdata[0:nobs]
    
    # generate discret support for epsilon to be used in Euler error
    # Eps are the central values
    # Phi are the associated probabilities
    Eps = np.zeros(npts);
    Cum = np.linspace(0.0, 1.0, num=npts+1)+.5/npts
    Cum = Cum[0:npts]
    Phi = np.ones(npts)/npts
    Eps = norm.ppf(Cum)

    neq = nx
    
    # initialize matrix of Euler errors
    Eerr = np.zeros((nobs,neq))
    # begin loop over time periods
    for t in range(0, nobs):
        # begin loop over possible va,lues of shock next period
        for i in range(0, npts):
            # find value of next period Z
            Zp = lfunc(Zdata[t,:],Eps[i],lpars)
            # find the value of X next period
            Xp = tfunc(Xdata[t,:],Zdata[t,:],tpars)
            # find the value of X in two periods
            Xpp = tfunc(Xp,Zp,tpars)
            # find the Euler errors
            Eerr[t,:] = Eerr[t,:] + Phi[i]*efunc(Xpp,Xp,Xdata[t,:], \
                Zp,Zdata[t,:],epars)
    return Eerr
