'''
Example Using the LinApp Bundle
Linear Approximation of a simple RBC model.  This is the neoclassical
growth model of Ramesey, Cass and Koopmans with a stochastic productivity
shock added.
X(t-1) = k(t) = k
X(t) = k(t+1) = kp
X(t+1) = k(t+2) = kpp
Y(t) is empty
Z(t) = z(t) = z
Z(t+1) = z(t+1) = zp
Definitions for y(t), w(t), r(t), c(t) and i(t) are given in the 
eample_def function.
The Euler equation is given in the example_dyn function
'''

import numpy as np
import matplotlib.pyplot as plt
from LinApp import LinApp_FindSS, LinApp_Deriv, LinApp_Solve, LinApp_SSL 

def example_def(kp, k, z, param):
	# calculate definitions
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
	Gamma = cnow**(-theta) - beta*cplus**(-theta)*(1+rplus-delta)/((1+g)**theta*(1+rho))
	return Gamma

# set parameter values
#  model parameters
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
# soution algorithm paramteres
nx = 1
ny = 0
nz = 1
Zbar = np.array([0.])
Z0 = Zbar
logX = 0   # do NOT log-linearize
Sylv = 1   # use the built-in sylvester solver

# program
nobs = 300      # number of periods in simulation
kstart = 1      # starting value for simulation (proportional to kbar)
solve = 1       # set to 1 to compute coeffs, 0 if coeffs already in memory

if solve == 1:
    # calculate steady state values
    guessXY = np.array([1.])
    kbar = LinApp_FindSS(example_dyn, param, guessXY, Zbar, nx, ny)
    print ('kbar value is ', kbar)
    zbar = Zbar
    ybar, cbar, ibar, rbar, wbar = example_def(kbar, kbar, zbar, param)
    reportbar = np.array([[ybar],[cbar],[ibar],[kbar],[wbar],[rbar]])
    
    # check SS values
    invec = np.concatenate([kbar, kbar, kbar, zbar, zbar])
    check = example_dyn(invec, param)
    print ('SS check value is ', check)
    
    # find derivatives
    [AA, BB, CC, DD, FF, GG, HH, JJ, KK, LL, MM] = \
        LinApp_Deriv(example_dyn, param, invec, nx, ny, nz, logX);
    print ('FF value is ', FF)
    print ('GG value is ', GG)
    print ('HH value is ', HH)
    print ('LL value is ', LL)
    print ('MM value is ', MM)
    
    # find policy function coefficients
    PP, QQ, RR, SS = LinApp_Solve(AA, BB, CC, DD, FF, GG, HH, JJ, KK, LL, MM, \
                                  phi, Zbar, Sylv);
    print ('PP value is ', PP)
    print ('QQ value is ', QQ)

# perform simulation
eps = np.random.randn(nobs)*sigma
z = np.zeros((nobs+1,1))
y = np.zeros(nobs)
c = np.zeros(nobs)
i = np.zeros(nobs)
r = np.zeros(nobs)
w = np.zeros(nobs)

z[0] = eps[0]
for t in range(0, nobs):
    z[t+1,:] = phi*z[t,:] + eps[t]
    
k, blank = LinApp_SSL(kstart*kbar, z, kbar, logX, PP, QQ, RR, SS)
for t in range(0, nobs):
    y[t], c[t], i[t], r[t], w[t] = example_def(k[t+1], k[t], z[t], param)
    
# delete last observation for k and z
k = k[0:nobs]
z = z[0:nobs]

# plot data
t = range(0, nobs)
plt.plot(t, z, label='z')
plt.plot(t, k, label='k')
plt.plot(t, y, label='y')
plt.plot(t, c, label='c')
plt.plot(t, i, label='i')
plt.plot(t, r, label='r')
plt.plot(t, w, label='w')
plt.xlabel('time')
plt.legend(loc=9, ncol=7, bbox_to_anchor=(0., 1.02, 1., .102))
plt.show(3)