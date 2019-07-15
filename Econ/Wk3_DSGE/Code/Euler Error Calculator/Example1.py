# -*- coding: utf-8 -*-
"""
Version 1.0
Mon Nov 28 2016
@author: Kerk L. Phillips
The code below is for demonstrating how to implement the EEcalc function for a 
simple RBC model that is solved and simulated using linear approximation about 
the steady state.


"""
import numpy as np
from LinApp_Deriv import LinApp_Deriv
from LinApp_Solve import LinApp_Solve
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
    invec = np.array([kpp[0,0], kp[0,0], k[0], zp[0], z[0]])
    outvec = example_dyn(invec, epars)
    return outvec
 
def example_tfunc(k, z, tpars):
    kp = PP*(k-kbar) + QQ*z + kbar
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
nobs = 300      # number of periods in simulatipon
kstart = 1      # starting value for simulation (proportional to kbar)
solve = 1       # set to 1 to compute coeffsd, 0 if coeffs already in memory


# calculate steady state values
kbar = (((1+rho)*(1+g)**theta-1+delta)/alpha)**(1/(alpha-1))
ybar = kbar**alpha
rbar = alpha*ybar/kbar
wbar = (1-alpha)*ybar
cbar = wbar + (1+rbar-delta)*kbar - (1+g)*(1+n)*kbar
ibar = ybar - cbar
reportbar = np.array([[ybar],[cbar],[ibar],[kbar],[wbar],[rbar]])

# check SS values
invec = np.array([kbar, kbar, kbar, 0, 0])
check = example_dyn(invec, param)
print('SS check', check)

# find derivatives
AA, BB, CC, DD, FF, GG, HH, JJ, KK, LL, MM = \
    LinApp_Deriv(example_dyn,param,invec,1,0,1,0);

# find policy function coefficients
PP, QQ, UU, RR, SS, VV = \
    LinApp_Solve(AA,BB,CC,DD,FF,GG,HH,JJ,KK,LL,MM,phi,0,1);

# perform simulation
tpars = (PP, QQ, kbar)
lpars = (phi)
eps = np.random.randn(nobs)*sigma
z = np.zeros((nobs+1))
k = np.zeros((nobs+1))
y = np.zeros(nobs)
r = np.zeros(nobs)
w = np.zeros(nobs)
i = np.zeros(nobs)
c = np.zeros(nobs)
k[0] = kbar*kstart
z[0] = eps[0]
for t in range(0, nobs):
    z[t+1] = example_lfunc(z[t], eps[t], lpars)
    k[t+1] = example_tfunc(k[t], z[t], tpars)
    y[t], c[t], i[t], r[t], w[t] = example_def(k[t+1], k[t], z[t], param)
k = k[0:nobs];
z = z[0:nobs];

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