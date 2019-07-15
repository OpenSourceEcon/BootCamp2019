# -*- coding: utf-8 -*-
"""
Created on Mon Nov 21 21:12:03 2016

@author: Kerk Phillips
"""

import pandas as pd
import numpy as np

def calcmom(data, means = False, stds = True, relstds = True, 
            corrs = True, autos = True, cvars = False):
    '''
    This function calculates a user controlled set of moments from data
    It takes the following as inputs:
    1)  data: an nobs-by-nvar matrix ofdata with with nobs indexing the 
        time-series observations in rows and nvar indexing the variables in 
        columns.  If GDP is included it should be in the first column.
    
    The function caculates and returns the following outputs in a single
    matrix with  rows for each type of momemnt and nvar columns .  The moments,
    if calculated are reported in the following order.  Each column can be
    turned off or on by setting the appropriate flags.
    1) means -     means of the variables
    2) stdevs -    standard deviations of variables
    3) relstds -   standard deviations of variables relative to GDP
    4) correls -   correlations of variables with GDP
    5) autocoors - autocorrelations of variables
    6) coefvars -  coefficients of variation (stdev/mean)
        
    Notes:
   
    '''
    (nobs, nvar) = data.shape
    report = np.zeros(nvar)
    count = 0
    rindex = []
    
    if means or cvars:
        dmeans= np.mean(data,axis=0)
        report = np.vstack((report,dmeans))
        count = count + 1
        rindex.extend(['means'])
        
    if stds or cvars:
        dstds = np.std(data,axis=0)
        report = np.vstack((report,dstds))
        count = count + 1
        rindex.extend(['standard deviations'])
        
    if relstds:
        drelstds = dstds/dstds[0]
        report = np.vstack((report,drelstds))
        count = count + 1
        rindex.extend(['standard deviations relative to GDP'])
        
    if corrs:
        temp = np.corrcoef(np.transpose(data))
        dcorrs = temp[0:nvar,0]
        report = np.vstack((report,dcorrs))
        count = count + 1
        rindex.extend(['corrleations with GDP'])
        
    if autos:
        dlead = np.transpose(data[1:nobs,:])
        dlag = np.transpose(data[0:nobs-1,:])
        temp = np.corrcoef(dlead,dlag)
        dautos = np.diag(temp, -nvar)
        report = np.vstack((report,dautos))
        count = count + 1
        rindex.extend(['autocorrelations'])
        
    if cvars:
        dcvars = dstds/dmeans
        report = np.vstack((report,dcvars))
        count = count + 1
        rindex.extend(['coefficients of variation'])
        
    # appropriately truncate lists
    report = report[1:count+1,:]
    rindex = rindex[0:count+1]
    return report, rindex

#test teh function with unfiltered data from hw05
    
macrodata = pd.read_csv('DSGEtestdata.csv')
# cols = list(macrodata.columns.values)
y = np.log(macrodata.RGDP)
c = np.log(macrodata.RCON)
i = np.log(macrodata.RINV)
g = np.log(macrodata.RGOV)
w = np.log(macrodata.NWAGE/macrodata.DEFL)
e = np.log(macrodata.EMP)
data = np.stack((y,c,i,g,w,e),axis=1)
varnames = ['y','c','i','g','w','e']
# filter data using log first differences
(nobs, nvar) = data.shape
data = np.log(data)
data = data[1:nobs,:] - data[0:nobs-1,:]


moms, momnames = calcmom(data,means=True, cvars=True)
print(moms)
table = pd.DataFrame(moms, index = momnames, columns = varnames) 
table.to_csv('DSGEmoments.csv')