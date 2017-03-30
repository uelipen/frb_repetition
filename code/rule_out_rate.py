#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 30 11:27:36 2017

@author: niels
"""

import numpy as np
import pylab as pl
import scipy.special as sp
import scipy.optimize as so
import mpmath

gammainc = np.vectorize(mpmath.gammainc)


k = 0.34
#theta = 5.7 #per day


def prob_of_nothing(k,theta,Delta):
    """
        Calculate the probability of seeing no burst in a continuous
        observation of duration Delta.
    """
    lam = 1./(theta*sp.gamma(1. + 1./k))
    out = gammainc(1./k,(Delta/lam)**k)
    out /= k
    out /= sp.gamma(1. + 1./k)
    return out


def rule_out_rate(k,Delta,alpha):
    """
        Calculate the rate that is ruled out at the alpha-percent level (for
        Weibull and Poisson).
    """
    
    def func(x):
        pon = prob_of_nothing(k,x,Delta)
        return pon - alpha/100.
    
    poissrate = -np.log(alpha/100.)/Delta
    rate = so.bisect(func,poissrate*1.e-5,poissrate*1.e5)
    return rate, poissrate


def make_plot(kvals,Delta,alpha):
    rateW = np.zeros(kvals.shape)
    rateP = np.zeros(kvals.shape)
    for i in range(len(kvals)):
        rateW[i], rateP[i] = rule_out_rate(kvals[i],Delta,alpha)
    pmrW, pmrP = rule_out_rate(k,Delta,alpha)
    pm_ratio = pmrW/pmrP
    print pm_ratio

    burst_color = '#E69F00'
    obs_color = '#009E73'
    poiss_color = '#56B4E9'

    fig = pl.figure(figsize=(5.,3.5))
    ax = fig.add_axes([0.13,0.15,0.86,0.84])
    pl.loglog(kvals,rateW/rateP,color=obs_color,lw=2)
    pl.plot(kvals,np.ones(len(kvals))*pm_ratio,'--',color=burst_color,lw=2)
    pl.plot([k,k],[1.e-1,1.e5],'--',color=burst_color,lw=2)
    pl.xlabel(r'$k$')
    pl.ylabel(r'$r_{5\%}^{\mathrm{(Weib)}}/r_{5\%}^{\mathrm{(Poiss)}}$')
    pl.xticks([0.2,0.34,1.0,3.],[r'0.2',r'      $\left<k\right>_{(k|N,t)}$',r'1.0',r'3.0'])
    pl.yticks([0.1,1,10,100,1000,10000],rotation='vertical')
    pl.savefig('rule_out_rate.pdf')
    return


kvals = np.logspace(-1.,1.,100)
Delta = 10.
alpha = 5.
make_plot(kvals,Delta,alpha)