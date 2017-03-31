#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sun Mar 12 14:51:08 2017

@author: niels
"""

import numpy as np
import scipy.special as sp
import scipy.interpolate as si
from scipy.optimize import bisect
import pylab as pl


def get_thresholds(post,thetavals,kvals):
    thetadiff = thetavals[1] - thetavals[0]
    kdiff = kvals[1] - kvals[0]

    def offset(threshold,prob,level):
        """
             Calcualte the fraction of points within an area defined by a threshold and
             return the difference between that value and the level defined above.
        """
        return prob[prob > threshold].sum()*thetadiff*kdiff - level
    
    thresh68 = bisect(offset,0.,post.max(),args=(post,0.68))
    thresh95 = bisect(offset,0.,post.max(),args=(post,0.95))
    thresh99 = bisect(offset,0.,post.max(),args=(post,0.99))
    
    return thresh68, thresh95, thresh99


def check_within_contours(post,thetavals,kvals,postval):
    thresh68, thresh95, thresh99 = get_thresholds(post,thetavals,kvals)
    return np.array([postval > thresh68, postval > thresh95, postval > thresh99]).reshape(3)


def make_plot(thetavals,kvals,postweib,postpoiss,ktrue=None,thetatrue=None,save=True,name='weibull_posterior_reparameterized.png'):
    
    thetadiff = thetavals[1] - thetavals[0]
    kdiff = kvals[1] - kvals[0]

    def offset(threshold,prob,level):
        """
             Calcualte the fraction of points within an area defined by a threshold and
             return the difference between that value and the level defined above.
        """
        return prob[prob > threshold].sum()*thetadiff*kdiff - level
    
    burst_color = '#E69F00'
    obs_color = '#009E73'
    poiss_color = '#56B4E9'

    fig = pl.figure(figsize=(5,5))
    ax = fig.add_axes((0.1,0.1,0.87,0.89))

    thresh68 = bisect(offset,0.,postweib.max(),args=(postweib,0.68))
    thresh95 = bisect(offset,0.,postweib.max(),args=(postweib,0.95))
    thresh99 = bisect(offset,0.,postweib.max(),args=(postweib,0.99))
    
    pl.contourf(10.**thetavals,10.**kvals,postweib.transpose(),[thresh99,post.max()],colors=obs_color,alpha=0.1)
    pl.contourf(10.**thetavals,10.**kvals,postweib.transpose(),[thresh95,post.max()],colors=obs_color,alpha=0.1)
    pl.contourf(10.**thetavals,10.**kvals,postweib.transpose(),[thresh68,post.max()],colors=obs_color,alpha=0.1)
    c_s = pl.contour(10.**thetavals,10.**kvals,postweib.transpose(),[thresh99,thresh95,thresh68],linewidths=1,colors=obs_color)

    thresh68 = bisect(offset,0.,postpoiss.max(),args=(postpoiss,0.68))
    thresh95 = bisect(offset,0.,postpoiss.max(),args=(postpoiss,0.95))
    thresh99 = bisect(offset,0.,postpoiss.max(),args=(postpoiss,0.99))
    
    pl.contourf(10.**thetavals,10.**kvals,postpoiss.transpose(),[thresh99,post.max()],colors=poiss_color,alpha=0.1)
    pl.contourf(10.**thetavals,10.**kvals,postpoiss.transpose(),[thresh95,post.max()],colors=poiss_color,alpha=0.1)
    pl.contourf(10.**thetavals,10.**kvals,postpoiss.transpose(),[thresh68,post.max()],colors=poiss_color,alpha=0.1)
    c_s = pl.contour(10.**thetavals,10.**kvals,postpoiss.transpose(),[thresh99,thresh95,thresh68],linewidths=1,colors=poiss_color)

    ax.set_xscale('log')
    ax.set_yscale('log')
    pl.xlabel(r'$r/(\mathrm{day}^{-1})$')
    pl.ylabel(r'$k$')
    pl.yticks([0.3,1.0,3.0],['0.3','1','3'],rotation='vertical')
    if (ktrue and thetatrue):
        pl.plot(10.**thetatrue,10.**ktrue,color=burst_color,marker='+',ms=10)

    if save:
        pl.savefig(name,transparent=True)
    else:
        pl.show()
    pl.close()
    return

    
thetavals = np.linspace(-1.0,2.0,30)
kvals = np.linspace(-1.0,0.5,30)
k = 0.345
lam = 0.038
theta = 1./(lam*sp.gamma(1. + 1./k))
k = np.log10(k)
theta = np.log10(theta)
within = np.zeros(3,dtype=int)
meanpostfinite = np.zeros((30,30))
for i in range(240):
    post = np.load('../../post_weibull_finite_%03i.npy'%i)
    postinterp = si.RectBivariateSpline(x=thetavals,y=kvals,z=post)
    postval = postinterp(theta,k)
    within += check_within_contours(post,thetavals,kvals,postval)
    meanpostfinite += post
#    thresh68, thresh95, thresh99 = get_thresholds(post,thetavals,kvals)
#    print postval, thresh68, thresh95, thresh99
#    c_s = pl.contour(thetavals,kvals,post.transpose(),[thresh99,thresh95,thresh68],linewidths=1)
#    pl.imshow(post.transpose()[::-1,:],interpolation='nearest',extent=[thetavals.min(),thetavals.max(),kvals.min(),kvals.max()])
#    pl.colorbar()
#    pl.xlabel(r'$\log \left({\theta}\,{\mathrm{day}}\right)$')
#    pl.yticks(rotation='vertical')
#    pl.ylabel(r'$\log k$')
#    pl.plot(theta,k,color='magenta',marker='+',ms=10)
#    pl.show()

meanpostfinite /= 240.

print 'finite sims:'
print 'within 68%:', within[0]/240.
print 'within 95%:', within[1]/240.
print 'within 99%:', within[2]/240.

within = np.zeros(3,dtype=int)
meanpostindependent = np.zeros((30,30))
for i in range(240):
    post = np.load('../../post_weibull_independent_%03i.npy'%i)
    postinterp = si.RectBivariateSpline(x=thetavals,y=kvals,z=post)
    postval = postinterp(theta,k)
    within += check_within_contours(post,thetavals,kvals,postval)
    meanpostindependent += post
#    thresh68, thresh95, thresh99 = get_thresholds(post,thetavals,kvals)
#    print postval, thresh68, thresh95, thresh99
#    c_s = pl.contour(thetavals,kvals,post.transpose(),[thresh99,thresh95,thresh68],linewidths=1)
#    pl.imshow(post.transpose()[::-1,:],interpolation='nearest',extent=[thetavals.min(),thetavals.max(),kvals.min(),kvals.max()])
#    pl.colorbar()
#    pl.xlabel(r'$\log \left({\theta}\,{\mathrm{day}}\right)$')
#    pl.yticks(rotation='vertical')
#    pl.ylabel(r'$\log k$')
#    pl.plot(theta,k,color='magenta',marker='+',ms=10)
#    pl.show()

meanpostindependent /= 240.

print 'independent sims:'
print 'within 68%:', within[0]/240.
print 'within 95%:', within[1]/240.
print 'within 99%:', within[2]/240.

make_plot(thetavals,kvals,meanpostindependent,meanpostfinite,ktrue=k,thetatrue=theta,save=True,name='finiteness.pdf')