# -*- coding: utf-8 -*-
"""
Created on Thu Dec 15 15:56:25 2016

@author: niels
"""


import numpy as np
import pylab as pl
#import jdcal
import astropy.time as at
import astropy.units as au
import astropy.coordinates as ac
from scipy.optimize import bisect
import scipy.special as sp
import mpmath


gammainc = np.vectorize(mpmath.gammainc)


FRB_loc = ac.SkyCoord('05:31:58','+33:08:04',unit=(au.hourangle,au.deg),equinox='J2000')
arecibo_loc = ac.EarthLocation.from_geodetic(lon=-66.7528,lat=18.3464) #index 0
effelsberg_loc = ac.EarthLocation.from_geodetic(lon='6:52:58',lat='50:31:29') #index 1
gbt_loc = ac.EarthLocation.from_geodetic(lon='-79.8398',lat='38.4322') #index 2
vla_loc = ac.EarthLocation.from_geodetic(lon='-107.6184',lat='34.0784') #index 3
lovell_loc = ac.EarthLocation.from_geodetic(lon='-2.3085',lat='53.2367') #index 4
locs = [arecibo_loc,effelsberg_loc,gbt_loc,vla_loc,lovell_loc] #locations of all involved observatories
FRB_DM = 559. #pc / cm^3
k_DM = 4.148808e3 #s MHz^2 cm^3 / pc


def log_weibull(x,theta,k):
    lam = 1./(theta*sp.gamma(1. + 1./k))
    return np.log(k/lam)[:,:,np.newaxis] + (k - 1.)[:,:,np.newaxis]*np.log(x[np.newaxis,np.newaxis,:]/lam[:,:,np.newaxis]) - (x[np.newaxis,np.newaxis,:]/lam[:,:,np.newaxis])**k[:,:,np.newaxis]


def log_cum_weibull(x,theta,k):
    lam = 1./(theta*sp.gamma(1. + 1./k))
    lcw = np.zeros((lam.shape[0],lam.shape[1],len(x)))
    for j in range(len(x)):
        print 'j =', j
        if x[j,2] == 0:
            lcw[:,:,j] = gammainc(1./k,(x[j,0]/lam)**k)/k/sp.gamma(1. + 1./k)
        elif x[j,2] == 1:
            lcw[:,:,j] = np.exp(-(x[j,0]/lam)**k)/(lam*sp.gamma(1. + 1./k))
        else:
            lcw[:,:,j] = np.exp(-(x[j,0]/lam)**k)
    return np.log(lcw)


def get_intervals():
    data = np.genfromtxt('observations.txt')
    times = np.genfromtxt('times.txt')[:,0]
    years = data[:,0]
    months = data[:,1]
    days = data[:,2]
    hours = data[:,3]
    minutes = data[:,4]
    seconds = data[:,5]
    duration = data[:,8]
    nFRB = data[:,9]
    tel = data[:,6]
    freq = data[:,7]
    starts = []
    ends = []
    intlengths = []
    intminlengths = []
    FRBcount = 0
    totFRB = int(nFRB.sum())
    for i in range(len(years)):
        startstr = '%04i-%02i-%02iT%02i:%02i:%02i'%(years[i],months[i],days[i],
                                                    hours[i],minutes[i],
                                                    seconds[i])
        start = at.Time(startstr,format='isot',scale='utc',
                        location=locs[int(tel[i])])
        dur = at.TimeDelta(duration[i],format='sec')
        dmcorr = at.TimeDelta(k_DM*FRB_DM/freq[i]**2,format='sec')
#        print dmcorr
        ltt_bary = start.light_travel_time(FRB_loc)
        start_bary = (start.tdb - dmcorr + ltt_bary)
        end_bary = (start.tdb - dmcorr + dur + ltt_bary)
        starts.append(start_bary.mjd)
        ends.append(end_bary.mjd)
        if int(nFRB[i]) == 0:
            intmin = end_bary.mjd - start_bary.mjd
            if ((FRBcount == 0) or (FRBcount == totFRB)):
                intmax = np.inf
            else:
                intmax = times[FRBcount] - times[FRBcount - 1]
            intminlengths.append([intmin,intmax,0])
        else:
            intmin = times[FRBcount] - start_bary.mjd
            if FRBcount > 0:
                intmax = times[FRBcount] - times[FRBcount - 1]
            else:
                intmax = np.inf
            intminlengths.append([intmin,intmax,1])
            FRBcount += 1
            for j in range(1,int(nFRB[i])):
                intlengths.append(times[FRBcount] - times[FRBcount - 1])
                FRBcount += 1
            intmin = end_bary.mjd - times[FRBcount - 1]
            if FRBcount < totFRB:
                intmax = times[FRBcount] - times[FRBcount - 1]
            else:
                intmax = np.inf
            intminlengths.append([intmin,intmax,2])
    return np.array(intlengths), np.array(intminlengths), np.array(starts), np.array(ends), times, nFRB


def plot_obs(start,end,times,nFRB):
    burst_color = '#E69F00'
    obs_color = '#009E73'
    
    nobs = len(start)
    burst = 0
    fig = pl.figure(figsize=(5,8.25))
    ax = fig.add_axes([0.02,0.07,0.96,0.92])
    for i in range(nobs):
        ax.plot(np.array([0.,end[i] - start[i]])*24.,[i,i],color=obs_color,lw=2)
        for j in range(int(nFRB[i])):
            ax.plot((times[burst] - start[i])*24.,i,'.',color=burst_color,ms=8)
            burst += 1
    pl.ylim(-0.5,nobs - 0.5)
    pl.xlabel(r'$t/\mathrm{hour}$')
    pl.tick_params(axis='y',labelleft=False,left=False,right=False)
    pl.xticks([0,1,2,3,4])
    pl.savefig('intervals.pdf',transparent=True)
#    pl.show()
    pl.close()
    return 0


def log_like(theta,k,intlengths,intminlengths):
    loglik = log_weibull(intlengths,10.**theta,10.**k).sum(axis=-1)
    loglik += log_cum_weibull(intminlengths,10.**theta,10.**k).sum(axis=-1)
    return loglik


def log_prior(theta,k):
#    val = -np.log(lam.max() - lam.min()) - np.log(k.max() - k.min())
#    return lam + k
    return 0.


def log_post(theta,k,intlengths,intminlengths):
    return log_like(theta,k,intlengths,intminlengths) + log_prior(theta,k)


def get_posterior(thetavals,kvals,intlengths,intminlengths):
#    diffs = get_differences()
#    print diffs
#    intlengths, intminlengths, start, end, times, nFRB = get_intervals()
    print intminlengths[:,0].min(), intminlengths[:,0].max()
    print intminlengths[:,1].min(), intminlengths[:,1].max()
    theta, k = np.meshgrid(thetavals,kvals,indexing='ij')
    post = np.exp(log_post(theta,k,intlengths,intminlengths))
    evidence = post.sum()*(thetavals[1] - thetavals[0])*(kvals[1] - kvals[0])
#    print 'evidence for Weibull:', evidence
    post /= evidence
    post_k = post.sum(axis=0)*(thetavals[1] - thetavals[0])
    post_theta = post.sum(axis=1)*(kvals[1] - kvals[0])
    print 'mean k:', (10.**kvals*post_k).sum()*(kvals[1] - kvals[0])
    print 'mean theta:', (10.**thetavals*post_theta).sum()*(thetavals[1] - thetavals[0])
    print 'mean interval:', (1./(10.**theta)*post).sum()*(thetavals[1] - thetavals[0])*(kvals[1] - kvals[0])
    return post, post_theta, post_k


def make_plot(thetavals,kvals,post,post_theta,post_k,ktrue=None,thetatrue=None,save=True,name='weibull_posterior_reparameterized.png'):
    
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
    
    burst_color = '#E69F00'
    obs_color = '#009E73'
    poiss_color = '#56B4E9'

    fig = pl.figure(figsize=(5,5))
    ax = fig.add_axes((0.1,0.1,0.6,0.6))
    pl.contourf(10.**thetavals,10.**kvals,post.transpose(),[thresh99,post.max()],colors=obs_color,alpha=0.1)
    pl.contourf(10.**thetavals,10.**kvals,post.transpose(),[thresh95,post.max()],colors=obs_color,alpha=0.1)
    pl.contourf(10.**thetavals,10.**kvals,post.transpose(),[thresh68,post.max()],colors=obs_color,alpha=0.1)
    c_s = pl.contour(10.**thetavals,10.**kvals,post.transpose(),[thresh99,thresh95,thresh68],linewidths=1,colors=obs_color)
    pl.xlabel(r'${r}/(\mathrm{day}^{-1})$')
    ax.set_xscale('log')
    ax.set_yscale('log')
    pl.yticks([0.2,0.3,0.4,0.5,0.6],['0.2','0.3','0.4','0.5','0.6'],rotation='vertical')
    pl.ylabel(r'$k$')
    if (ktrue and thetatrue):
        pl.plot(thetatrue,ktrue,color=burst_color,marker='+',ms=10)

    ax = fig.add_axes((0.1,0.7,0.6,0.29))
    ax.set_xscale('log')
    pl.plot(10.**thetavals,post_theta,color=obs_color)
    pl.tick_params(axis='x',labelbottom=False)
    pl.yticks([1,2,3],rotation='vertical')
    pl.ylabel(r'$\mathcal{P}(\log r|N,t)$')
#    pl.xlim(thetavals.min(),thetavals.max())
    poiss_post = (3.3*(10.**thetavals))**17.*np.exp(-3.3*(10.**thetavals))
    poiss_post /= poiss_post.sum()*(thetavals[1] - thetavals[0])
    pl.ylim(0.,1.1*poiss_post.max())
    pl.plot(10.**thetavals,poiss_post,color=poiss_color)
    if thetatrue:
        pl.plot([thetatrue,thetatrue],[0.,10.],color=burst_color)
    
    ax = fig.add_axes((0.7,0.1,0.29,0.6))
    ax.set_yscale('log')
    pl.plot(post_k,10.**kvals,color=obs_color)
    pl.tick_params(axis='y',labelleft=False)
    pl.xticks([1,3,5])
    pl.xlabel(r'$\mathcal{P}(\log k|N,t)$')
    pl.ylim(10.**kvals.min(),10.**kvals.max())
    pl.xlim(0.,1.1*post_k.max())
    if ktrue:
        pl.plot([0.,10.],[ktrue,ktrue],color=burst_color)
    
#    pl.show()
    if save:
        pl.savefig(name,transparent=True)
    else:
        pl.show()
    pl.close()
    return

    
def get_confident(cum,xvals):
    i = 0
    while cum[i] < 0.16:
        i += 1
    lower = xvals[i]
#    i = 0
#    while cum[i] < 0.5:
#        i += 1
#    print xvals[i]*conv
    i = 0
    while cum[i] < 0.84:
        i += 1
    upper = xvals[i]
    return np.array([lower,upper])


if __name__ == '__main__':
    intlengths, intminlengths, start, end, times, nFRB = get_intervals()
    
    print 'total observing time:', intlengths.sum() + intminlengths[:,0].sum()
    print 'total number of bursts:', nFRB.sum()
    print 'mean difference between bursts:', intlengths.mean()
    
    for i in range(len(start)):
        for j in range(len(start)):
            if ((start[i] > start[j]) & (start[i] < end[j])):
                print 'Observation %02i starts within observation %02i.'%(i,j)
            if ((end[i] > start[j]) & (end[i] < end[j])):
                print 'Observation %02i ends within observation %02i.'%(i,j)
        
    plot_obs(start,end,times,nFRB)
    thetavals = np.linspace(-1.0,2.0,100)
    kvals = np.linspace(-0.8,-0.2,100)
    post, post_theta, post_k = get_posterior(thetavals,kvals,intlengths,intminlengths)
    np.save('post.npy',post)
    np.save('post_theta.npy',post_theta)
    np.save('post_k.npy',post_k)
#    print post.shape
    post = np.load('post.npy')
    post_theta = np.load('post_theta.npy')
    post_k = np.load('post_k.npy')
    k = 10.**((kvals*post_k).sum()*(kvals[1] - kvals[0]))
    theta = 10.**((thetavals*post_theta).sum()*(thetavals[1] - thetavals[0]))
    print 'k =', k
    print 'theta =', theta
    theta_cum = post_theta.cumsum()*(thetavals[1] - thetavals[0])
    k_cum = post_k.cumsum()*(kvals[1] - kvals[0])
    print 'k interval:', 10.**get_confident(k_cum,kvals)
    print 'theta interval:', 10.**get_confident(theta_cum,thetavals)
    make_plot(thetavals,kvals,post,post_theta,post_k,ktrue=k,thetatrue=theta,save=True,name='2dpost.pdf')
    