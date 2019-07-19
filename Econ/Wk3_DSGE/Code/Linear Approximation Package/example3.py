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

def example_def(kp, k, ell, z, x, param):
    # calculate definitions
    y = k**alpha*(np.exp(z)*ell)**(1-alpha)
    w = (1-alpha)*y/ell
    r = alpha*y/k
    c = w + (1+r-delta)*k - kp*(1+g)*(1+n)
    i = y - c
    return y, c, i, r, w

def example_dyn(invec, param):
    # unpack in
    kplus = invec[0]
    know = invec[1]
    kminus = invec[2]
    ellplus = invec[3]
    ellnow = invec[4]
    zplus = invec[5]
    xplus = invec[6]
    znow = invec[7]
    xnow = invec[8]

    # get definitions each period
    ynow, cnow, inow, rnow, wnow = \
        example_def(know, kminus, ellnow, znow, xnow, param)
    yplus, cplus, iplus, rplus, wplus = \
        example_def(kplus, know, ellplus, zplus, xplus, param)

    # calculate Gamma function
    Gamma1 = cnow**(-theta)*wnow - chi*ellnow**lambd
    Gamma2 = cnow**(-theta) - \
        beta*cplus**(-theta)*(1+rplus-delta)/((1+g)**theta*(1+rho))
        
    return np.array([Gamma1, Gamma2])

# set parameter values
#  model parameters
g = .025
n = .01
delta = .08
alpha = .33
theta = 2.5
lambd = .33
chi = 1.
rho = .05
phi_z = .9
sigma_z = .0075
phi_x = .9
sigma_x = .0075
beta = (1+g)**(1-theta)*(1+n)/(1+rho)
param = [g, n, delta, alpha, theta, lambd, chi, rho, phi_z, sigma_z, phi_x, \
         sigma_x, beta]
# soution algorithm paramteres
nx = 1
ny = 1
nz = 2
Zbar = np.array([0., 0.])
Z0 = Zbar
logX = 0   # do NOT log-linearize
Sylv = 0   # use the built-in sylvester solver

# program
nobs = 300      # number of periods in simulation
kstart = 1      # starting value for simulation (proportional to kbar)
solve = 1       # set to 1 to compute coeffs, 0 if coeffs already in memory

if solve == 1:
    # calculate steady state values
    guessXY = np.array([10., .5])
    bar = LinApp_FindSS(example_dyn, param, guessXY, Zbar, nx, ny)
    [kbar, ellbar] = bar
    print ('kbar value is ', kbar)
    print ('lbar value is ', ellbar)
    zbar = np.array([0.])
    xbar = np.array([0.])
    ybar, cbar, ibar, rbar, wbar = example_def(kbar, kbar, ellbar, zbar, \
                                               xbar, param)
    reportbar = np.array([[ybar], [cbar], [ibar], [kbar], [ellbar], [wbar], \
                          [rbar]])
    
    kbar = np.array([kbar])
    ellbar = np.array([ellbar])
    # check SS values
    invec = np.concatenate([kbar, kbar, kbar, ellbar, ellbar, zbar, xbar, \
                            zbar, xbar])
    check = example_dyn(invec, param)
    print ('SS check value is ', check)
    
    # find derivatives
    [AA, BB, CC, DD, FF, GG, HH, JJ, KK, LL, MM] = \
        LinApp_Deriv(example_dyn, param, invec, nx, ny, nz, logX);
    print ('AA value is ', AA)
    print ('BB value is ', BB)
    print ('CC value is ', CC)
    print ('DD value is ', DD)
    print ('FF value is ', FF)
    print ('GG value is ', GG)
    print ('HH value is ', HH)
    print ('JJ value is ', JJ)
    print ('KK value is ', KK)
    print ('LL value is ', LL)
    print ('MM value is ', MM)
    
    NN = np.array([[phi_z, 0.], [0., phi_x]])
    
    # find policy function coefficients
    PP, QQ, RR, SS = LinApp_Solve(AA, BB, CC, DD, FF, GG, HH, JJ, KK, LL, MM, \
                                  NN, Zbar, Sylv);
    print ('PP value is ', PP)
    print ('QQ value is ', QQ)

# perform simulation
eps_z = np.random.randn(nobs)*sigma_z
eps_x = np.random.randn(nobs)*sigma_x
z = np.zeros(nobs+1)
x = np.zeros(nobs+1)
y = np.zeros(nobs)
c = np.zeros(nobs)
i = np.zeros(nobs)
r = np.zeros(nobs)
w = np.zeros(nobs)

z[0] = eps_z[0]
x[0] = eps_x[0]
for t in range(0, nobs):
    z[t+1] = phi_z*z[t] + eps_z[t]
    x[t+1] = phi_x*x[t] + eps_x[t]
    
Z = np.vstack((z, x))
Z = Z.T
    
k, ell = LinApp_SSL(kstart*kbar, Z, bar, logX, PP, QQ, RR, SS)
for t in range(0, nobs):
    y[t], c[t], i[t], r[t], w[t] = example_def(k[t+1,0], k[t,0], ell[t,0], \
        z[t], x[t], param)
    
# delete last observation for X and Z
k = k[0:nobs]
ell = ell[0:nobs]
z = z[0:nobs]
x = x[0:nobs]

# plot data
t = range(0, nobs)
plt.plot(t, z, label='z')
plt.plot(t, x, label='x')
plt.plot(t, k, label='k')
plt.plot(t, ell, label='l')
plt.plot(t, y, label='y')
plt.plot(t, c, label='c')
plt.plot(t, i, label='i')
plt.plot(t, r, label='r')
plt.plot(t, w, label='w')
plt.xlabel('time')
#plt.legend(loc=9, ncol=8, bbox_to_anchor=(0., 1.02, 1., .102))
plt.show()