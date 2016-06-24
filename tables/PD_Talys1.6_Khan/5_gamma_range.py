from numpy import *
import os

# -----------------------------------------------
# Collect exclusive gamma cross sections from TALYS for remaining files after thinning
# -----------------------------------------------

data = genfromtxt('xs_gamma_thin.txt')
xssum = genfromtxt('xs_sum_daughter.txt')

for i,(Z,N,Zd,Nd) in enumerate(xssum[:,:4]):
    print Z, N, Zd, Nd
    idx = (data[:,0] == Z) * (data[:,1] == N) * (data[:,2] == Zd) * (data[:,3] == Nd)
    xs = data[:,5:][idx]
    xs_sum = xssum[i,4:]
    # calculate branching ratios, set to 0 when total cross-section = 0
    br = xs / xs_sum
    data[:,5:][idx][isinf(br)] = 0
    data[:,5:][idx][isnan(br)] = 0

savetxt('xs_gamma_thin.txt', data, fmt='%i\t%i\t%i\t%i\t%.4f' + '\t%g'*301)
