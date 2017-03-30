#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sun Mar 12 20:00:24 2017

@author: niels
"""

from simulations import *
factorial = np.vectorize(mpmath.factorial)

burst_color = '#E69F00'
obs_color = '#009E73'
poiss_color = '#56B4E9'
    
#obstimes = np.array([[0.,2.6/24.],[1.,1.+2.6/24.],[5.,5.+2.6/24.],[6.,6.+1.5/24.],[7.,7.+1.5/24.],[8.,8.+1.5/24.],[9.,9.+1.5/24.],[15.,15.+2.6/24.],[31.,31.+2.6/24.],[39.,39.+2.6/24.],[43.,43.+2.6/24.],[45.,45.+2.6/24.],[46.,46.+2.6/24.]])
obstimes = np.array([[0.,10./theta]])
start = obstimes[:,0]
end = obstimes[:,1]
samples = 100000
n = np.zeros(samples)
for j in range(samples):
    intlengths, intminlengths, n[j] = sim_finite(start,end)
bins = 30
a,b,c = pl.hist(n,range=(-0.5,bins - 0.5),bins=bins,normed=True)
pl.close()
fig = pl.figure(figsize=(5.,3.5))
ax = fig.add_axes([0.1,0.15,0.89,0.84])
x = np.arange(bins)
mu = (end - start).sum()/(lam*sp.gamma(1. + 1./k))
poiss = mu**x*np.exp(-mu)/factorial(x)
xpr = np.zeros(2*bins)
poisspr = np.zeros(2*bins)
weibpr = np.zeros(2*bins)
xpr[::2] = x - 0.5
poisspr[::2] = poiss
weibpr[::2] = a
xpr[1::2] = x + 0.5
poisspr[1::2] = poiss
weibpr[1::2] = a
pl.plot(xpr,weibpr,color=obs_color,lw=2)
pl.plot(xpr,poisspr,color=poiss_color,lw=2)
pl.xlim((-0.5,bins - 0.5))
pl.xlabel(r'$N$')
pl.ylabel(r'$P(N|k,r)$')
pl.yticks([0.,0.1,0.2],rotation='vertical')
#pl.savefig('n_dana.pdf')
pl.savefig('n_singleint.pdf')
