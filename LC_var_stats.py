#
#  Chris Suberlak 04/30/2016 
#  
# OVERVIEW : 
# Code to calculate variability parameters on processed forced photometry
# lightcurves.  We calculate sigma, mu, chi2, based on full lightcurves,
# and only for those lightcurves where sigma>0 or chi2 > 1  we calculate 
# seasonal averages too . 
# 
# INPUT : processed forced photometry lightcurves, arranged by objectId , with 
# ['objectId', 'mjd', 'psfFlux', 'psfFluxErr', 'flag','faintMean','faintMedian','faintTwoSigma']
#
# OUTPUT: a file with full lightcurve variability characteristics : 
# [ u'N', u'avgMJD', u'chi2DOF', u'chi2R', u'muFull', u'psfFluxErrMean',
#  u'psfFluxMean', u'psfFluxMedian', u'psfFluxSigma', u'psfFluxSkew',
#        u'rangeMJD', u'sigmaFull' ]pwd
#
#  and a file with seasonal variability characteristics :
# 
import sys
sys.path.insert(0, '/astro/users/suberlak/S13Agg_analysis/packages/')

from astroML.stats import median_sigmaG
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os 
import variabilityFunctions as varF
from astropy.time import Time 
pd.options.mode.chained_assignment = None  # to avoid  the following warning:
# A value is trying to be set on a copy of a slice from a DataFrame.
# Try using .loc[row_indexer,col_indexer] = value instead
# see http://stackoverflow.com/questions/20625582/how-to-deal-with-this-pandas-warning




DirIn = '/astro/store/scratch/tmp/suberlak/s13_stripe82/forced_phot_lt_23/NCSA/'
DirOut  = '/astro/store/scratch/tmp/suberlak/s13_stripe82/forced_phot_lt_23/NCSA/Var/'

l = os.listdir(DirIn)


#for name in l : 
#        print('Processing filter_patch file %s' % name)
name = l[0]
print('Read  LC file %s'%name)
fp_data = pd.read_csv(DirIn+name ,delimiter=',', usecols=['objectId','mjd','psfFlux','psfFluxErr','flag','faintMean','faintMedian','faintTwoSigma'])

####  first drop all NaNs  in psfFlux...
fp_data = fp_data[np.isfinite(fp_data['psfFlux'])]

####  replace all psfFlux  where SNR < 2  with  faintMean
rows = fp_data['flag'] == 1
fp_data.ix[rows, 'psfFlux'] = fp_data.ix[rows, 'faintMean']

# group by objectId to calculate full LC variability characteristics 
grouped = fp_data.groupby('objectId')

# Calculate all variability stats  on full LC  (sigma_full,  mu_full,   chi2 )
varMetricsFull = grouped.apply(varF.computeVarMetrics)

# Calculate magnitudes based on psfFluxMean   psfFluxMedian , etc. 

def flux2absigma(flux, fluxsigma):
    """Compute AB mag sigma given flux and flux sigma"""
    return FIVE_OVER_2LOG10 * fluxsigma / flux;


def flux2ab(flux):
    """Compute AB mag given flux"""
    return -2.5 * np.log10(flux) - 48.6;


varMetricsFull['psfMean'] = flux2ab(varMetricsFull['psfFluxMean'])
varMetricsFull['psfMedian'] = flux2ab(varMetricsFull['psfFluxMedian'])
varMetricsFull['psfMeanErr'] = flux2absigma(varMetricsFull['psfFluxMean'],varMetricsFull['psfFluxErrMean'])
varMetricsFull['psfMedianErr'] = flux2absigma(varMetricsFull['psfFluxMedian'],varMetricsFull['psfFluxErrMedian'])



# Save the result 
path = DirOut + 'Var'+name
varMetricsFull.to_csv(path)


# Check how many sources  are variable... 
# i.e. have either  sigma > 0 or chi2 > 1

#http://docs.scipy.org/doc/numpy-1.10.0/reference/generated/numpy.ma.mask_or.html
m1 = (varMetricsFull['sigmaFull'] > 0).values
m2 =  (varMetricsFull['chi2DOF'] > 1).values
m3 = (varMetricsFull['chi2R'] > 1).values
m= np.ma.mask_or(m3, np.ma.mask_or(m1,m2))
print('Out of %d objects,  %d fulfill  sigma>0 or chi2R>1 or chi2DOF>1' % (len(m), np.sum(m)))

# save the diagnostics about variability to a file...
diagFile = path+'.diag'
file = open(diagFile, "w")
file.write('Input :  \n')
s = '    '+ DirIn + '\n' 
file.write(s)
s = '    '+ name + '\n\n'
file.write(s)
s = 'There are '+str(np.sum(m)) + ' points out of ' + str(len(m)) + ' that fulfill sigma>0 or chi2R>1 or chi2DOF>1 \n '
file.write(s)
file.write('Output :  \n')
s = '    '+ DirOut + '\n' 
file.write(s)
s = '    '+ 'Var'+name + '\n\n'  
file.write(s)
file.close()   


# Grab names of variable objects, that fulfill the criteria above... 
varObjectIds = varMetricsFull[m].index

# Grab only those rows that correspond to variable objects...
rows = np.in1d(fp_data['objectId'].values, varObjectIds)
fp_var = fp_data.ix[rows] 

# make a new column to designate seasons...
fp_var['season'] = np.nan

# I insert a very early date at the beginning of the list, so that all obs between
# 1990 and 2005 are clustered together, and 2005-2006, 2006-2007, etc are 
# averaged seasonally (season starts August 1st)

dates = [str(year)+'-08-01 00:00:00.000' for year in np.arange(2005,2009)]
dates.insert(0,'1990-08-01 00:00:00.000' )
cutDates = Time(dates, format='iso')
seasons = np.arange(len(cutDates))+0 

# Assign value of a season for each row...
for i in range(len(cutDates.mjd)-1):
    mask = (fp_var['mjd'].values > cutDates[i].mjd) * (fp_var['mjd'].values < cutDates[i+1].mjd)
    fp_var.ix[mask, 'season'] = seasons[i]  
 
# Calculate seasonal properties for objects that are variable (based on their full LC...) 
grouped = fp_var.groupby(['objectId','season'])

varMetricsSeasonal = grouped.apply(varF.ComputeSeasons)

path = DirOut + 'SeasVar'+name
varMetricsSeasonal.to_csv(path)
print('Saving Seasonal Variability statistics to %s'%path)


def ComputeVarMetricsSeasonal(group): 
    
    ''' A function to calculate averages for the full lightcurve based on lightcurve averages   
    '''
    #mu,sigma = get_mu_sigma(group['psfFlux'].values*1e27, group['psfFluxErr'].values*1e27)
    #rangeMJD = group['mjd'].values.max() - group['mjd'].values.min() 
    
    return pd.Series({'Nseasons':group['psfFluxMean'].count(),
                      'psfFluxMeanMean': group['psfFluxMean'].mean(),
                      'psfFluxErrMeanMean' : group['psfFluxErrMean'].mean(),
                      'chi2DOFmean' : calcChi2raw(group['psfFluxMean'].values,group['psfFluxErrMean'].values),
                      'chi2DOFmedian' : calcChi2raw(group['psfFluxMedian'].values,group['psfFluxErrMean'].values),
                      'chi2Rmean' : calcChi2robust(group['psfFluxMean'].values,group['psfFluxErrMean'].values),
                      'chi2Rmedian' : calcChi2robust(group['psfFluxMedian'].values,group['psfFluxErrMean'].values)
                     })


grouped  = varMetricsSeasonal.groupby(level=0)
#grouped.get_group(grouped.groups.keys()[0])['psfFluxMean'].values

varMetricsFullSeasonal = grouped.apply(ComputeVarMetricsSeasonal)



