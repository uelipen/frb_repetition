#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sun Mar 12 16:15:41 2017

@author: niels
"""


import numpy as np
import pylab as pl
import scipy.special as sp
from nifty import *

np.random.seed(101)

k = 0.34
theta = 5.7
lam = 1./(theta*sp.gamma(1. + 1./k))


def get_Nbar(Delta_t,samples=10000):
    N = np.zeros(samples)
    for i in range(samples):
        #draw time of first burst:
        u = np.random.gamma(shape=1./k,scale=1.)
        t = lam*u**(1./k)
        while t < Delta_t:
            N[i] += 1
            t += lam*np.random.weibull(a=k)
    return N.mean()


def get_powspec(tmax=300./24.,nbins=300,samples=10000):
    grid = rg_space(num=nbins,dist=tmax/nbins)
    klen, rho, pindex, pundex = grid.get_codomain().get_power_indices(log=True)
    power = np.zeros(klen.shape)
    for i in range(samples):
        times = []
        u = np.random.gamma(shape=1./k,scale=1.)
        t = lam*u**(1./k)
        while t < tmax:
            times.append(t)
            t += lam*np.random.weibull(a=k)
        binned_counts = np.histogram(times,range=(0.,tmax),bins=nbins)[0]*1.
        binned_counts /= tmax/nbins
        n = field(grid,val=binned_counts)
        power += n.power(bare=True,log=True,pindex=pindex)/tmax
        if (i % (samples/10) == 0):
            print '%02i percent done'%((100*i)/samples)
    power /= samples
    return klen, power


klenls, powerls = get_powspec(tmax=100000.,nbins=2000000,samples=1000)
out = np.array([klenls.flatten(),powerls.flatten()])
np.save('power_k0.34_theta5.7_logls.npy',out)

klenss, powerss = get_powspec(tmax=10.,nbins=2000000,samples=1000)
out = np.array([klenss.flatten(),powerss.flatten()])
np.save('power_k0.34_theta5.7_logss.npy',out)

ls = np.load('power_k0.34_theta5.7_logls.npy')
ss = np.load('power_k0.34_theta5.7_logss.npy')

kls = ls[0]
kss = ss[0][ss[0] > kls.max()]
powls = ls[1]
powss = ss[1][ss[0] > kls.max()]
nu = np.append(kls,kss)
powspec = np.append(powls,powss)

poissval = theta
weibval = poissval*(sp.gamma(1. + 2./k)/sp.gamma(1. + 1./k)**2 - 1.)

burst_color = '#E69F00'
obs_color = '#009E73'
poiss_color = '#56B4E9'

fig = pl.figure(figsize=(5.,3.5))
ax = fig.add_axes([0.1,0.15,0.89,0.84])

#pl.loglog([nu[1],nu[-1]],[poissval,poissval],'--',lw=2,color=burst_color)
pl.loglog([nu[1],nu[-1]],[weibval,weibval],'--',lw=2,color=burst_color)
#pl.loglog([theta,theta],[poissval/2.,weibval*2.],'--',lw=2,color=burst_color)
pl.loglog([theta,theta],[2.e-1,2.e2],'--',lw=2,color=burst_color)
alpha = (theta*sp.gamma(1.+1./k))**k
#func = nu**(k-1.) #- alpha*nu**(2.*k - 1.) + 0.5*alpha**2*nu**(3.*k - 1.)
#pl.loglog(nu,func/func[-1]*(powspec[-1]-poissval),'--',lw=2,color=poiss_color)
pl.loglog(nu,powspec - poissval,lw=2,color=obs_color)

#pl.ylim(poissval/2.,weibval*2.)
pl.ylim(2.e-1,2.e2)
pl.xlim(nu[1],nu[-1])
pl.xlabel(r'$\nu/(\mathrm{day}^{-1})$')
pl.ylabel(r'$P(\nu)$')
pl.yticks(rotation='vertical')
pl.savefig('powspec_weibull.pdf')
pl.close()
