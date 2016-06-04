# -*- coding: iso-8859-1 -*-
#
# Functions used to calculate variability parameters 
#
# needed by LC_var_stats.py 


from astroML.stats import median_sigmaG
import pandas as pd
import numpy as np

def gaussgauss_logL(xi, ei, mu, sigma):
    """Equation 5.63: gaussian likelihood with gaussian errors
    Source code https://github.com/astroML/astroML/blob/master/astroML/
    """
    ndim = len(np.broadcast(sigma, mu).shape)

    xi = xi.reshape(xi.shape + tuple(ndim * [1]))
    ei = ei.reshape(ei.shape + tuple(ndim * [1]))

    s2_e2 = sigma ** 2 + ei ** 2
    return -0.5 * np.sum(np.log(s2_e2) + (xi - mu) ** 2 / s2_e2,
                         -1 - ndim)


def approximate_mu_sigma(xi, ei, axis=None):
    """Estimates of mu0 and sigma0 via equations 5.67 - 5.68"""

    if axis is not None:
        xi = np.rollaxis(xi, axis)
        ei = np.rollaxis(ei, axis)
        axis = 0

    mu_approx, sigmaG = median_sigmaG(xi, axis=axis)
    e50 = np.median(ei, axis=axis)
    var_twiddle = (sigmaG ** 2 + ei ** 2 - e50 ** 2)
    sigma_twiddle = np.sqrt(np.maximum(0, var_twiddle))

    med = np.median(sigma_twiddle, axis=axis)
    mu = np.mean(sigma_twiddle, axis=axis)

    zeta = np.ones_like(mu)
    zeta[mu != 0] = med[mu != 0] / mu[mu != 0]

    var_approx = zeta ** 2 * sigmaG ** 2 - e50 ** 2
    sigma_approx = np.sqrt(np.maximum(0, var_approx))

    return mu_approx, sigma_approx



def get_mu_sigma(xi,ei, N_boot=1000):
    ''' A short function
    to calculate a full mu, sigma, based 
    on N_boot bootstraps of the given sample
    with approximate method, 
    and then calculating the 2D maximum of the
    log-likelihood calculated on the space spanning
    from the minimum to the maximum value of the bootstrapped 
    result. 
    
    Input:
    xi : array of measurement values (assumed flux, i.e. order 1e-27)
    ei ; array of measurement errors (same order of mag as xi)
    N_boot : integer, representing the number of bootstraps 
    Returns:
    mu_max : mu calculated as a maximum of the 2D log-likelihood 
    sig_mag : sigma calculated as a maximum of the 2D log-likelihood 
    '''

    # Calculate bootstrapped approximate.... 
    #print N_boot
    indices = np.random.randint(0, len(xi), (len(xi), N_boot))

    xi_boot = xi[indices]
    ei_boot = ei[indices]

    mu_boot, sigma_boot = approximate_mu_sigma(xi_boot, ei_boot, 0)


    # Calculate marginalized likelihood sigma and mu
    # using as boundaries the minimum and maximum result of the 
    # bootstrapped approximate calculation 
 
    max_factor = 1.0
    sigma = np.linspace(min(sigma_boot), max_factor*max(sigma_boot), 70)
    mu = np.linspace(min(mu_boot), max_factor*max(mu_boot), 70)
    
    if (len(mu) == 0) | (len(sigma) == 0) : 
        # for some reason we may be completely missing it... 
        return np.nan, np.nan 
    else : 
        logL = gaussgauss_logL(xi, ei, mu, sigma[:, np.newaxis])
        logL -= logL.max()
        ind = np.where(logL == np.max(logL))
    
        if len(ind) < 2  : 
            # for some reason, we may be unable to find the maximum...
            return np.nan, np.nan
              
        else : 
            # return the mu and sigma at the maximum of the likelihood
            # (note : I assume log-likelihood is smooth, and has only 
            # one maximum )
            return mu[ind[1]][0], sigma[ind[0]][0]

def calcChi2raw(y, yerr):
    """Compute simple reduced chi2  (mean-based) if more than 1 datapoints present
    chi2 = np.sum(((flux-meanFlux)**2.0) / (fluxErr ** 2.0)) / (N-1.0)
    """
    N = len(y)
    if N < 2:
        return np.nan
    else:
        chi2 = np.sum(((y-np.mean(y))/yerr)**2)
        return chi2/(N-1)
    
def calcChi2robust(y, yerr):
    """Compute simple robust reduced  chi2 (percentile-based) if more than 1 datapoints present
    """
    N = len(y)
    if N < 2:
        return np.nan
    else:
        Z = (y - np.median(y)) / yerr
        chi2R = 0.7414 * (np.percentile(Z,75)-np.percentile(Z,25))
        return chi2R

def calcMedian(y):
    ''' Calculate the median as 50-th percentile '''
    N = len(y)
    if N == 1 : 
        return float(y)    
    elif N == 0 : 
        return np.nan
    else : 
        return float(np.percentile(y,50))

def calcWeightedMean(y,yerr):
    ''' Calculate the weighted mean '''
    N = len(y)
    if N == 1 : 
        return float(y)    
    elif N == 0 : 
        return np.nan
    else: 
        # weights = 1 / (yerr ** 2.0)  
        # wMean = np.sum(weights * flux) / np.sum(weights)
        return float(np.add.reduce(y / (yerr * yerr)) / np.add.reduce((1/yerr)*(1/yerr)))

def calcWeightedMeanErr(yerr):
    ''' The error on the weighted mean '''
    N = len(yerr)
    if N == 1 : 
        return float(yerr)    
    elif N == 0 : 
        return np.nan
    else : 
        # weights = 1 / (yerr ** 2.0)
        # sigma_mean = 1 / sqrt(sum(weights))
        return np.sqrt(1. / np.add.reduce((1 / yerr)*(1/yerr)))

def calcWeightedStDev(y, yerr, yWmean):
    ''' Calculate the  weighted standard deviation
    Needs the weighted mean on y to be calculated 
    beforehand, using eg. calcWeightedMean()
    I'm using Bessel's correction to make it unbiased ...
    
    # calculate  weighted standard deviation corrected for intrinsic scatter 
    # using http://www.itl.nist.gov/div898/software/dataplot/refman2/ch2/weightsd.pdf
    # Yusra uses 1/N-1 instead of N/N-1.... calcWStdCorrAndMean
    # I'm pretty confused having read https://en.wikipedia.org/wiki/Bessel's_correction
    # and https://en.wikipedia.org/wiki/Weighted_arithmetic_mean#Weighted_sample_variance
    # after that http://stats.stackexchange.com/questions/6534/how-do-i-calculate-a-weighted-standard-deviation-in-excel I'm done. 

    '''
    N = len(y)
    if N == 1:
        return float(yerr)
    elif N == 0:
        return np.nan 
    else :     
        weights=1.0 / ( yerr *yerr)
        return np.sqrt((N / (N-1.0) ) * (np.sum(weights * ((y - yWmean) ** 2.0)) / np.sum(weights)))  

def calcSigmaG(y):
    ''' Calculate the  interquartile sigma.'''
    N = len(y)
    if N == 1:
        return float(y)
    elif N == 0:
        return np.nan  
    else: 
        q25,q75 = np.percentile(y, (25,75))
        return 0.7413 * (q75-q25)

def computeVarMetrics(group):
    ''' Variability metrics to compute for each object on full lightcurve  / 
    all points in a given season 
    
    '''
    # print diagnostic for figuring out error...
    #print 'objectId= ', group['objectId'].values[0]
    
    # even though I drop NaNs before, I do it here explicitly to save 
    # me from headaches 
    # eg. obj  216199180459189485    216199180459189485
    # have one row with NaN , not caught by other filters... 

    group.dropna(subset=['psfFlux', 'psfFluxErr'], inplace=True)

    # calculate range of dates in a given lightcurve    
    rangeMJD = group['mjd'].values.max() - group['mjd'].values.min() 
    
    
    Flux = group['psfFlux'].values
    FluxErr = group['psfFluxErr'].values
    
    # calculate Weighted  Mean
     
    FluxMean = calcWeightedMean(Flux,FluxErr)
    FluxMeanErr = calcWeightedMeanErr(FluxErr)
    psfFluxStDev = calcWeightedStDev(Flux,FluxErr, FluxMean)

    # calculate Median error : not necessary (since it's just const * meanErr)
    # medianErr = np.sqrt(np.pi / 2.0) * FluxMeanErr

    # note  : get_mu_sigma() gets digestion problems with NaN's - make sure 
    # that Flux and FluxErr are free of NaN's ! 
    # I multiply flux by a factor for stability issues    
    N = len(Flux)
    if N == 0 : 
        mu = np.nan
        sigma = np.nan
    elif N == 1  :
        mu, sigma = Flux[0], 0
    else : 
        mu, sigma = get_mu_sigma(Flux*1e27, FluxErr*1e27,1000)

    # set the flag about length...
    if N > 10 : 
        flagLtTenPts = np.nan
    else:
        flagLtTenPts = 1 
   
    
    return pd.Series({'N':group['psfFlux'].count(),
                      'psfFluxMean': FluxMean,
                      'psfFluxMeanErr' : FluxMeanErr,
                      'psfFluxMedian': calcMedian(Flux),
                      'psfFluxMedianErr': np.sqrt(np.pi / 2)*FluxMeanErr,
                      'psfFluxSkew' : group['psfFlux'].skew(),
                      'psfFluxSigG' : calcSigmaG(Flux),
                      'psfFluxStDev':psfFluxStDev,
                      'chi2DOF' : calcChi2raw(Flux,FluxErr),
                      'chi2R' : calcChi2robust(Flux,FluxErr),
                      'sigmaFull' :sigma,
                      'muFull' :mu,
                      'avgMJD' : group['mjd'].mean(),
                      'rangeMJD' : rangeMJD,
                      'flagLtTenPts' : flagLtTenPts
                     })


def ComputeVarFullBinned(group): 
    
    ''' A function to calculate averages for the full lightcurve binned into seasons 
    '''    
    Flux= group['psfFluxMean'].values
    FluxErr =group['psfFluxMeanErr'].values

    N = len(Flux)
    if N == 0 : 
        mu = np.nan
        sigma = np.nan
    elif N == 1  :
        mu, sigma = Flux, 0
    else : 
        mu, sigma = get_mu_sigma(Flux*1e27, FluxErr*1e27, 1000)

    FluxMean = calcWeightedMean(Flux,FluxErr)
    FluxMeanErr = calcWeightedMeanErr(FluxErr)
    psfFluxStDev = calcWeightedStDev(Flux,FluxErr, FluxMean)

    return pd.Series({'Nseasons':group['psfFluxMean'].count(),
                      'psfFluxMeanMean': FluxMean, 
                      'psfFluxMeanMeanErr' :FluxMeanErr,
                      'chi2DOFmean' : calcChi2raw(Flux,FluxErr),
                      'chi2DOFmedian' : calcChi2raw(group['psfFluxMedian'].values,group['psfFluxMedianErr'].values),
                      'chi2Rmean' : calcChi2robust(Flux,FluxErr),
                      'chi2Rmedian' : calcChi2robust(group['psfFluxMedian'].values,group['psfFluxMedianErr'].values),
                      'sigmaFull' :sigma,
                      'muFull' :mu
                     })


# http://stackoverflow.com/questions/18603270/progress-indicator-during-pandas-operations-python


def logged_apply(g, func, *args, **kwargs):
    step_percentage = 100. / len(g)
    import sys
    sys.stdout.write('apply progress:   0%')
    sys.stdout.flush()

    def logging_decorator(func):
        def wrapper(*args, **kwargs):
            progress = wrapper.count * step_percentage
            sys.stdout.write('\033[D \033[D' * 4 + format(progress, '3.0f') + '%')
            sys.stdout.flush()
            wrapper.count += 1
            return func(*args, **kwargs)
        wrapper.count = 0
        return wrapper

    logged_func = logging_decorator(func)
    res = g.apply(logged_func, *args, **kwargs)
    sys.stdout.write('\033[D \033[D' * 4 + format(100., '3.0f') + '%' + '\n')
    sys.stdout.flush()
    return res

from pandas.core.groupby import DataFrameGroupBy
DataFrameGroupBy.logged_apply = logged_apply

