# -*- coding: utf-8 -*-
"""
Version 1.0
Mon Nov 28 2016
@author: Kerk L. Phillips
The code below is for demonstrating how to implement the EEcalc function for a 
simple RBC model that is solved and simulated using value-function iteration
on a grid.


"""
import numpy as np
from EulerErrors import EEcalc
    
def example_def(kp, k, z, param):
    # calculate definitions for GDP, wages, rental rates, consumption
    np.exp((1-alpha)*z)
    y = k**alpha*np.exp((1-alpha)*z)
    w = (1-alpha)*y
    r = alpha*y/k
    c = w + (1+r-delta)*k - kp*(1+g)*(1+n)
    i = y - c
    return y, c, i, r, w

def example_dyn(invec, param):
    # unpack in
    kplus = invec[0]
    know = invec[1]
    kminus = invec[2]
    zplus = invec[3]
    znow = invec[4]

    # get definitions each period
    ynow, cnow, inow, rnow, wnow = example_def(know, kminus, znow, param)
    yplus, cplus, iplus, rplus, wplus = example_def(kplus, know, zplus, param)

    # calculate Gamma function
    Gamma = cnow**(-theta) / ((cplus**(-theta)*(1+rplus-delta)/((1+g)**theta* \
                            (1+rho)))) - 1.0
    return Gamma
    
def example_efunc(kpp, kp, k, zp, z, epars):
    invec = np.array([kpp, kp, k[0], zp[0], z[0]])
    outvec = example_dyn(invec, epars)
    return outvec
 
def example_tfunc(k, z, tpars):
    k = np.reshape(k,(1))
    z = np.reshape(z,(1))
    if gridsim:
        # perform simulation with grid
        kidx = np.argmin(np.abs(kvec-k)) # index of closest value for k
        zidx = np.argmin(np.abs(zvec-z)) # index of closest value for z
        kp = trans[zidx,kidx];
    else:
        # perform simulation with polynomial fit
        temp =  np.stack(([1.0], k, k**2, k**3, z, z**2, z**3, k*z, k**2*z, \
                              k*z**2))
        kp = np.dot(coeffs,temp)
    return kp
     
def example_lfunc(z, eps, lpars):
    zp = phi*z + sigma*eps
    return zp

# demonstrate the EEcalc function for a simple RBC model
# set parameter values
#  model
g = .025
n = .01
delta = .08
alpha = .33
theta = 2.5
rho = .05
phi = .9
sigma = .0075
beta = (1+g)**(1-theta)*(1+n)/(1+rho)
param = [g, n, delta, alpha, theta, rho, phi, sigma, beta]

# program
nobs = 300         # number of periods in simulation
kstart = 1         # starting value for simulation (proportional to kbar)
low = .2           # low end of k grid (proportional to kbar)
high = 5           # high end of k grid (proportional to kbar)
spread = 3         # spread of the z grid (proportional to sigma)
nptsk = 101        # number of points in the k grid
nptsz = 31         # number of points in the z grid
converge = .0001;  # convergence criterion for value-function iteration
gridsim = 0        # set to 1 to use grid, 0 to use polynomial fit
solve = 1          # set to 1 to compute grid, 0 if grid already in memory

if solve == 1:
    # calculate steady state values
    kbar = (((1+rho)*(1+g)**theta-1+delta)/alpha)**(1/(alpha-1))
    ybar = kbar**alpha
    rbar = alpha*ybar/kbar
    wbar = (1-alpha)*ybar
    cbar = wbar + (1+rbar-delta)*kbar - (1+g)*(1+n)*kbar
    ibar = ybar - cbar
    reportbar = np.array([[ybar],[cbar],[ibar],[kbar],[wbar],[rbar]])
    
    # set up grid for k
    klow = low*kbar   # low end of grid
    khigh = high*kbar # high end of grid
    kincr = np.log(khigh/klow)/(nptsk-1);
    kvec = np.zeros(nptsk)
    kvec[0] = klow
    for i in range (1, nptsk):
        kvec[i] = kvec[i-1]*(1+kincr)
        
    # set up grid for z
    zlow = -spread*max(sigma, .001)
    zhigh = spread*max(sigma, .001)
    zincr = (zhigh-zlow)/(nptsz-1);
    zvec = np.zeros(nptsz)
    zvec[0] = zlow
    for i in range (1, nptsz):
        zvec[i] = zvec[i-1] + zincr
        
    # create meshgrid
    kmesh, zmesh = np.meshgrid(kvec, zvec)
    
    # find value function and transition function
    distance = 1.0;
    maxwhile = 100;
    count = 0;
    value = np.zeros((nptsz,nptsk));
    newval = np.zeros((nptsz,nptsk));
    trans = np.zeros((nptsz,nptsk));
    print(count, distance)
    while distance > converge:
        count = count + 1
        if count > maxwhile:
            break
        for i in range(0, nptsz): # for all values of z(t)
            # find index of closest value for E{z(t+1)}
            m = np.argmin(np.abs(zvec - zvec[i]*phi))
            for j in range(0, nptsk): # for all values of k(t)
                maxval = -10000000000;
                for l in range(0, nptsk): # search over all values of k(t+1)
                    if theta == 1:
                        temp = np.log(kvec[j]**alpha*np.exp((1-alpha)*zvec[i]) 
                        +(1-delta)*kvec[j]-kvec[l]*(1+g)*(1+n)) + beta*value[m,l]
                    else:
                        temp = (((kvec[j]**alpha*np.exp((1-alpha)*zvec[i]) 
                        +(1-delta)*kvec[j]-kvec[l]*(1+g)*(1+n))**(1-theta)-1)) \
                        / (1-theta) + beta*value[m,l]
                    # print i, j, temp
                    if np.iscomplex(temp):
                       temp = -100000000
                    if np.isnan(temp):
                        temp = -100000000
                    if temp > maxval:
                        maxval = temp
                        newval[i,j] = temp
                        trans[i,j] = kvec[l]
        # print newval
        distance = np.mean(np.abs(value/newval - 1.0))
        print count, distance
        for i in range(0, nptsz):
            for j in range(0, nptsk):
                value[i,j] = newval[i,j]
    
    # fit a polynomial
    Y = trans.flatten()
    
    X = np.ones(nptsk*nptsz)
    
    temp = kmesh.flatten()
    X = np.vstack((X,temp))
    
    temp = kmesh**2
    temp = temp.flatten()
    X = np.vstack((X,temp))
    
    temp = kmesh**3
    temp = temp.flatten()
    X = np.vstack((X,temp))
    
    temp = zmesh.flatten()
    X = np.vstack((X,temp))
    
    temp = zmesh**2
    temp = temp.flatten()
    X = np.vstack((X,temp))
    
    temp = zmesh**3
    temp = temp.flatten()
    X = np.vstack((X,temp))
    
    temp = kmesh*zmesh
    temp = temp.flatten()
    X = np.vstack((X,temp))
    
    temp = kmesh**2*zmesh
    temp = temp.flatten()
    X = np.vstack((X,temp))
    
    temp = kmesh*zmesh**2
    temp = temp.flatten()
    X = np.vstack((X,temp))
    
    XtX = np.dot(X,np.transpose(X))
    XtY = np.dot(X,Y)
    coeffs = np.dot(np.linalg.inv(XtX),XtY)
    tpoly = np.zeros((nptsz,nptsk))
    for i in range(0, nptsz):
        for j in range(0, nptsk):
            tpoly[i,j] = np.dot(np.stack((1, kvec[j], kvec[j]**2, kvec[j]**3, \
            zvec[i], zvec[i]**2, zvec[i]**3, \
            kvec[j]*zvec[i], kvec[j]**2*zvec[i], kvec[j]*zvec[i]**2)),coeffs)
            
    # calcuate R-squared
    Rsq = 1 - np.sum((trans-tpoly)**2)/np.sum(trans**2)
    print 'R-squared', Rsq

# perform simulation
eps = np.random.randn(nobs)*sigma
z = np.zeros(nobs+1)
k = np.zeros(nobs+1)
k[0] = kstart*kbar
z[0] = eps[0]
tpars = [kvec, zvec, trans, coeffs]
lpars = [phi, sigma]
for t in range(0, nobs):
    k[t+1] = example_tfunc(k[t], z[t], tpars)
    z[t+1] = example_lfunc(z[t], eps[t], lpars)    

    
# remove final k & z
k = k[0:nobs]
z = z[0:nobs]



# find Euler Errors
Xdata = k.reshape(len(k), 1)
Zdata = z.reshape(len(z), 1)
efunc = example_efunc
epars = param
tfunc = example_tfunc
lfunc = example_lfunc
EErrs = EEcalc(Xdata, Zdata, efunc, epars, tfunc, tpars, lfunc, lpars)

MaxAbsEE = np.max(np.abs(EErrs))
MeanAbsEE = np.mean(np.abs(EErrs))
RootMeanSqEE = (np.mean(EErrs**2))**.5

print('Euler Error Summary Statistics')
print('Maximum Absolute Euler Error:', MaxAbsEE)
print('Mean Absolute Euler Error:', MeanAbsEE)
print('Root Mean Squared Euler Error:', RootMeanSqEE)