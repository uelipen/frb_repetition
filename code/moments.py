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


klenls, powerls = get_correlation_fct(tmax=300/24./10**i,nbins=2000000,samples=10)
out = np.array([klenls.flatten(),powerls.flatten()])
np.save('power_k0.34_theta5.7_logls.npy',out)
