from numpy import *
import os

# -----------------------------------------------
# Set if cross section is larger than inelastic cross section cross section value to inelastic cross section value
# -----------------------------------------------

data = genfromtxt('xs_sum_daughter.txt')
data_thin = genfromtxt('xs_gamma_thin_renormalized.txt')

for i,channel in enumerate(data_thin):
    Z = channel[0]
    N = channel[1]
    Zd = channel[2]
    Nd = channel[3]
    print Z, N, Zd, Nd
    idx = (data[:,0] == Z) * (data[:,1] == N) * (data[:,2] == Zd) * (data[:,3] == Nd)
    xs = data[:,4:][idx]
    xs = xs[0,:]

    # renormalize cross sections
    for j in range(0,len(data_thin[0,5:])):
        if (data_thin[i,5+j] > xs[j]):
            data_thin[i,5+j] = xs[j]

savetxt('xs_gamma_thin_renormalized_optimized.txt', data_thin, fmt='%i\t%i\t%i\t%i\t%.4f' + '\t%g'*301)
