#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sun Mar 12 17:48:03 2017

@author: niels
"""

from simulations import *

nobs = 10
tobs = 1./theta #days
#tgap = 24./24. #days

steps = 30
poneormore = np.zeros(steps)
nbar = np.zeros(steps)
tgaps = np.logspace(0.,np.log10(60),steps)/24.#np.arange(steps)/24.
samples = 1000000
for i in range(steps):
    print i
    tgap = tgaps[i]
    start = np.arange(nobs)*(tobs+tgap)
    end = start + tobs
    n = np.zeros(samples)
    for j in range(samples):
        intlengths, intminlengths, n[j] = sim_finite(start,end)
    poneormore[i] = float((n > 0).sum())/samples
    nbar[i] = n.mean()
#pl.hist(n,range=(-0.5,39.5),bins=40,normed=True)
#pl.show()
burst_color = '#E69F00'
obs_color = '#009E73'
poiss_color = '#56B4E9'
fig = pl.figure(figsize=(5.,3.5))
ax = fig.add_axes([0.1,0.15,0.89,0.84])
#pl.plot(tgaps,nbar)
#pl.show()
pl.semilogx(tgaps,poneormore,color=obs_color,lw=2)
#pl.show()
pl.xlabel(r'time between observations / days')
pl.ylabel(r'$P(N > 0)$')
pl.xlim((tgaps.min(),tgaps.max()))
pl.ylim(0.73,0.98)
pl.yticks([0.75,0.85,0.95],rotation='vertical')
pl.savefig('obs_strategy.pdf')

