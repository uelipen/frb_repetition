#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sun Mar 12 21:37:26 2017

@author: niels
"""

import numpy as np
import pylab as pl


def plot_interval(Delta,times):
    burst_color = '#E69F00'
    obs_color = '#009E73'
    
    fig = pl.figure(figsize=(5.0,1.5))
    ax = fig.add_axes([0.02,0.29,0.96,0.7])
    ax.plot([0.,times[-2] - 1.5],[1,1],color=obs_color,lw=2)
    ax.plot([times[-2] - 1.5,times[-2] - 0.5],[1,1],':',color=obs_color,lw=2)
    ax.plot([times[-2] - 0.5,Delta],[1,1],color=obs_color,lw=2)
    ax.plot(times,np.ones(len(times)),'.',color=burst_color,ms=8)
    pl.ylim(0,2)
    pl.xlabel(r'$t$')
    pl.tick_params(axis='y',labelleft=False,left=False,right=False)
    xt = [0]
    xtl = [0]
    for i in range(len(times) - 1):
        xt.append(times[i])
        if i == len(times) - 2:
            xtl.append(r'$t_N$')
        else:
            xtl.append(r'$t_%i$'%(i+1))
    xt.append(Delta)
    xtl.append(r'$\Delta$')
    xt.append(times[-1])
    xtl.append(r'$t_{N+1}$')
    pl.xticks(xt)
    ax.set_xticklabels(xtl)
    pl.xlim(-1,times[-1] + 1)
    pl.savefig('singleinterval.pdf',transparent=True)
#    pl.show()
    pl.close()
    return 0


plot_interval(6.,np.array([1.5,2.5,3.,5.2,7.5]))