#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Mar  8 16:33:21 2017

@author: niels
"""

from weibull_reparameterized import *


intlengths, intminlengths, start, end, times, nFRB = get_intervals()
end -= start.min()
start -= start.min()
k = 0.34
theta = 5.7
lam = 1./(theta*sp.gamma(1.+1./k))
#k = 1.
#lam = 3.3/17.


def sim_independent(start,end):
    Nobs = len(start)
    intlengths = []
    intminlengths = []
    obstimes = end - start
    nFRB = 0
    for i in range(Nobs):
        #draw time of first burst:
        u = np.random.gamma(shape=1./k,scale=1.)
        tnew = lam*u**(1./k)
        if tnew < obstimes[i]:
            intminlengths.append([tnew,0.,1])
            nFRB += 1
        else:
            intminlengths.append([obstimes[i],0.,0])
        #draw successive times, checking in each step whether still within interval
        while tnew < obstimes[i]:
            told = tnew
            tnew += lam*np.random.weibull(a=k)
            if tnew < obstimes[i]:
                intlengths.append(tnew - told)
                nFRB += 1
            else:
                intminlengths.append([obstimes[i],0.,2])
    return np.array(intlengths), np.array(intminlengths), nFRB


def sim_finite(start,end):
    Nobs = len(start)
    intlengths = []
    intminlengths = []
    nobs = np.zeros(Nobs)
    tFRBs = [[] for i in range(Nobs)]

    #draw time of first burst:
    u = np.random.gamma(shape=1./k,scale=1.)
    t = lam*u**(1./k)
    for c in range(Nobs):
        if ((start[c] < t) & (end[c] >= t)):
            nobs[c] += 1
            if nobs[c] == 1:
                intminlengths.append([t - start[c],0.,1])
            else:
                intlengths.append(t - tFRBs[c][-1])
            tFRBs[c].append(t)
    #draw successive bursts:
    while t < end.max():
        t += lam*np.random.weibull(a=k)
        for c in range(Nobs):
            if ((start[c] < t) & (end[c] >= t)):
                nobs[c] += 1
                if nobs[c] == 1:
                    intminlengths.append([t - start[c],0.,1])
                else:
                    intlengths.append(t - tFRBs[c][-1])
                tFRBs[c].append(t)
    for c in range(Nobs):
        if nobs[c] == 0:
            intminlengths.append([end[c] - start[c],0.,0])
        else:
            intminlengths.append([end[c] - tFRBs[c][-1],0.,2])
    return np.array(intlengths), np.array(intminlengths), nobs.sum()


if __name__ == '__main__':
    thetavals = np.linspace(-1.0,2.0,30)
    kvals = np.linspace(-1.0,0.5,30)
    
    np.random.seed(101)
    
    for i in range(300):
        intlengths, intminlengths, n = sim_independent(start,end)
        print i
        post, post_theta, post_k = get_posterior(thetavals,kvals,intlengths,intminlengths)
        np.save('../../sims/post_weibull_independent_%03i.npy'%i,post)
        np.save('../../sims/post_theta_weibull_independent_%03i.npy'%i,post_theta)
        np.save('../../sims/post_k_weibull_independent_%03i.npy'%i,post_k)
#        make_plot(thetavals,kvals,post,post_theta,post_k,ktrue=k,thetatrue=1./(lam*sp.gamma(1.+1./k)),save=True,name='../../sims/post_weibull_independent_%03i.png'%i)

    np.random.seed(101)
    
    for i in range(300):
        intlengths, intminlengths, n = sim_finite(start,end)
        print i
        post, post_theta, post_k = get_posterior(thetavals,kvals,intlengths,intminlengths)
        np.save('../../sims/post_weibull_finite_%03i.npy'%i,post)
        np.save('../../sims/post_theta_weibull_finite_%03i.npy'%i,post_theta)
        np.save('../../sims/post_k_weibull_finite_%03i.npy'%i,post_k)
#        make_plot(thetavals,kvals,post,post_theta,post_k,ktrue=k,thetatrue=1./(lam*sp.gamma(1.+1./k)),save=True,name='../../sims/post_weibull_finite_%03i.png'%i)
