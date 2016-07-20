from numpy import *
import os

# -----------------------------------------------
# Renormalize the cross sections of the photon emission channels which remain after thinning
# that their sum before and after the thinning is equal 
# -----------------------------------------------

data = genfromtxt('xs_excl_gamma.txt')
daughter = vstack({tuple(r) for r in data[:,:4]})
data_thin = genfromtxt('xs_gamma_thin.txt')

for (Z,N,Zd,Nd) in daughter:
    print Z, N, Zd, Nd
    idx = (data[:,0] == Z) * (data[:,1] == N) * (data[:,2] == Zd) * (data[:,3] == Nd)
    idy = (data_thin[:,0] == Z) * (data_thin[:,1] == N) * (data_thin[:,2] == Zd) * (data_thin[:,3] == Nd)
    xs = data[:,5:][idx]
    xs_thin = data_thin[:,5:][idy]
    xs_sum = sum(xs, axis=0)
    xs_sum_thin = sum(xs_thin, axis=0)

    # renormalize cross sections
    renormalization = xs_sum/xs_sum_thin
    renormalization[isnan(renormalization)] = 0
    data_thin[:,5:][idy] *= renormalization 

savetxt('xs_gamma_thin_renormalized.txt', data_thin, fmt='%i\t%i\t%i\t%i\t%.4f' + '\t%g'*301)
