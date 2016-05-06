#
#  Chris Suberlak 04/30/2016 
#  
# OVERVIEW : 
# Code to process the S82 forced  photometry lightcurves by calculating  
# 2sigma, mean, median, based on the truncated Gaussian distribution 
# 
# INPUT : raw forced photometry lightcurves, arranged by objectId , with 
# ['objectId', 'mjd', 'psfFlux', 'psfFluxErr' ] 
#
# OUTPUT: processed forced photometry lightcurves,  with unchanged flux, etc,
# but with additional columns  :
# ['objectId', 'mjd', 'psfFlux', 'psfFluxErr', 'flag','faintMean','faintMedian','faintTwoSigma']  



# Note : to run on magneto do
# %cpaste  in jupyter notebook 
#

# make all necessary imports....
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os 
pd.options.mode.chained_assignment = None
# to avoid  http://stackoverflow.com/questions/20625582/how-to-deal-with-this-pandas-warning

# faint source treatment 
from scipy.stats import norm
from scipy.special import erf

# variability 
import sys
sys.path.insert(0, '/astro/users/suberlak/S13Agg_analysis/packages/')
import variabilityFunctions as varF
reload(varF)
from astroML.stats import median_sigmaG
from astropy.time import Time 



def calculate_mean(psfFlux, psfFluxErr):
    xObs = psfFlux / psfFluxErr
    # Calculate my results
    xMean = (1/ (norm.sf(-xObs )*np.sqrt(2*np.pi))) * np.exp(-(xObs**2.0) / 2.0) + xObs
    F_mean = xMean * psfFluxErr
    return F_mean

def calculate_median(psfFlux, psfFluxErr):
    x0= -psfFlux / psfFluxErr
    return psfFluxErr*norm.ppf((1+norm.cdf(x0))/2.0) + psfFlux


def calculate_2sigma(psfFlux, psfFluxErr):
    x0= -psfFlux / psfFluxErr
    return psfFlux + psfFluxErr * norm.isf(0.05 * norm.sf(x0))

def calculate_rms(psfFlux, psfFluxErr):
    xObs = psfFlux / psfFluxErr
    xMean = (1/ (norm.sf(-xObs )*np.sqrt(2*np.pi))) * np.exp(-(xObs**2.0) / 2.0) + xObs
    delX = xObs - xMean
    I1 = norm.sf(-xObs)
    I0bysig2 = 0.5*erf(xObs/np.sqrt(2)) + (1.0/np.sqrt(2*np.pi))*np.exp(-(xObs**2.0) / 2.0)*(2*delX - xObs) + 0.5 + delX*delX*norm.sf(-xObs)
    xRMS = np.sqrt(I0bysig2 / I1)
    return  xRMS * psfFluxErr



#####
####  PROCESSING FILES IN   /astro/store/pogo4/s13_stripe82/forced_phot_lt_23/
####  


DirIn = '/astro/store/pogo4/s13_stripe82/forced_phot_lt_23/NCSA/'
DirOut = '/astro/store/scratch/tmp/suberlak/s13_stripe82/forced_phot_lt_23/NCSA/'


lProc = []
lProc += [each for each in os.listdir(DirOut) if each.endswith('.csv')]

lProc = [name[5:] for name in lProc]
lToDo = os.listdir(DirIn)
lToDoComp = [name[:-3] for name in lToDo]

maskDone = np.in1d(lToDoComp,lProc)

lToDoFilt =  np.array(lToDoComp)[~maskDone]

for name in lToDoFilt : 
	print('Processing filter_patch file %s' % name)

	# read in the raw lightcurve... 
	fp_data = pd.read_csv(DirIn+name+'.gz', compression='gzip',  usecols=['objectId', 'mjd', 'psfFlux', 'psfFluxErr'])
	#
	##########  STEP 1 : single-epoch data ###########  
	#

	####  first drop all NaNs  in psfFlux...      
	m1  = np.isnan(fp_data['psfFlux'])  # True if NaN  
	m2 =  ~np.isfinite(fp_data['psfFlux']) #  True if not finite  
	m  = m1 | m2  # a logical or 
	if np.sum(m) > 0 :  # only apply if there is anything to drop ... 
	    fp_data.drop(m.index[m], inplace=True)  # drop entire rows 
	 
	#### check to make sure that there are no NaN psfFluxErr... 
	m1  = np.isnan(fp_data['psfFluxErr'])  # True if NaN  
	m2 =  ~np.isfinite(fp_data['psfFluxErr']) #  True if not finite  
	m  = m1 | m2  # a logical or 
	if np.sum(m) > 0 :  # only apply if there is anything to drop ... 
	    fp_data.drop(m.index[m], inplace=True)

	# make a new column, fill with 0's
	fp_data['flagFaint'] = 0

	# mask those rows that correspond to SNR < 2
	mask = (fp_data['psfFlux'].values / fp_data['psfFluxErr'].values) < 2

	# print info how many points are affected
	print('There are %d points of %d that have SNR<2' %(np.sum(mask),len(mask)))

	# set flag at those rows to 1
	fp_data.ix[mask, 'flagFaint'] = 1

	# make new columns for  Mean  Median  2 sigma...
	fp_data['faintMean'] = np.nan
	fp_data['faintMedian'] = np.nan
	fp_data['faintTwoSigma'] = np.nan
	fp_data['faintRMS'] = np.nan
	# calculate the faint replacement only for faint points...
	fp_data.ix[mask, 'faintMean'] = calculate_mean(fp_data['psfFlux'][mask].values,fp_data['psfFluxErr'][mask].values)
	fp_data.ix[mask, 'faintMedian'] = calculate_median(fp_data['psfFlux'][mask].values,fp_data['psfFluxErr'][mask].values)
	fp_data.ix[mask, 'faintTwoSigma'] = calculate_2sigma(fp_data['psfFlux'][mask].values,fp_data['psfFluxErr'][mask].values)
	fp_data.ix[mask, 'faintRMS'] = calculate_rms(fp_data['psfFlux'][mask].values,fp_data['psfFluxErr'][mask].values)
	## If one wants, at this stage we may replace all fluxes by Mean or Median
	# I would keep both, and then use either, as appropriate.
	# This replaces all rows where flag == 1
	#rows = fp_data['flag'] == 1
	#fp_data.ix[rows, 'psfFlux'] = fp_data.ix[rows, 'faintMean']
	# 
	######################### SAVING OUTPUT 
	# as a new csv....
	#path = DirOut+'Proc_'+name
	#print('Saving processed file to %s'% path)
	#fp_data.to_csv(path)  

	# Save the diagnostics..
	path = DirOut+'Proc_'+name
	diagFile = path+'.diag'
	file = open(diagFile, "w")
	file.write('Input :  \n')
	s = '    '+ DirIn + '\n' 
	file.write(s)
	s = '    '+ name + '\n\n'
	file.write(s)
	s = 'There are '+str(np.sum(mask)) + ' points out of ' + str(len(mask)) + ' that have SNR < 2 \n '
	file.write(s)
	s = '( SNR = psfFlux / psfFluxErr ) \n \n'
	file.write(s)
	file.write('Output :  \n')
	s = '    '+ DirOut + '\n' 
	file.write(s)
	s = '    '+ 'Proc_'+name + '\n\n'  
	file.write(s)
	file.close()   

	#
	##########  STEP 2 : Derived Quantities ###########  
	#

	####  replace all psfFlux  where SNR < 2  with  faintMean  
	rows = fp_data['flagFaint'] == 1
	fp_data.ix[rows, 'psfFlux'] = fp_data.ix[rows, 'faintMean']

	# group by objectId to calculate full LC variability characteristics 
	grouped = fp_data.groupby('objectId')

	#
	#  Calculate low-order statistical properties of full lightcurves:
	# 
	#  psfFluxMean, psfFluxMeanErr, psfFluxMedian, psfFluxMedianErr
	#  psfFluxStDev, psfFluxSigG , psfFluxSkew, avgMJD, rangeMJD   
	#  chi2DOF , chi2R,  sigmaFull, muFull, flagLtTenPts
	# 
	varMetricsFull = grouped.apply(varF.computeVarMetrics)

	# Calculate magnitudes based on average fluxes :
	# psfMean  psfMedian  psfMeanErr  psfMedianErr 

	def flux2absigma(flux, fluxsigma):
	    """Compute AB mag sigma given flux and flux sigma"""
	    FIVE_OVER_2LOG10 = 1.085736204758129569
	    return FIVE_OVER_2LOG10 * fluxsigma / flux;


	def flux2ab(flux):
	    """Compute AB mag given flux"""
	    return -2.5 * np.log10(flux) - 48.6;



	varMetricsFull['psfMean'] = flux2ab(varMetricsFull['psfFluxMean'])
	varMetricsFull['psfMedian'] = flux2ab(varMetricsFull['psfFluxMedian'])
	varMetricsFull['psfMeanErr'] = flux2absigma(varMetricsFull['psfFluxMean'],varMetricsFull['psfFluxMeanErr'])
	varMetricsFull['psfMedianErr'] = flux2absigma(varMetricsFull['psfFluxMedian'],varMetricsFull['psfFluxMedianErr'])

	#
	######################### SAVING OUTPUT        ######################### 
	# 
	path = DirOut +'Var/'+ 'Var'+name
	varMetricsFull.to_csv(path)
	print('Saving Full, unbinned  LC statistics to %s'%path)

	#
	##########  STEP 3 : Variable Candidates ###########  
	# 

	#
	# for variable candidates we calculate  metrics for 
	# points per season, as well as bin the lightcurve 
	# into seasons, and calculate the metrics of binned 
	# lightcurve 
	# 

	# Check how many sources  are variable... 
	# i.e. have either  sigma > 0 or chi2 > 1
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
	 
	# Calculate seasonal metrics for objects that are variable (based on their full LC...) 
	grouped = fp_var.groupby(['objectId','season'])

	varMetricsSeasonal = grouped.apply(varF.computeVarMetrics)
	#
	######################### SAVING OUTPUT        ######################### 
	#
	path = DirOut + 'SeasVar'+name
	varMetricsSeasonal.to_csv(path)
	print('Saving Seasonal statistics to %s'%path)

	# Calculate binned lightcurve metrics (binned by seasons)
	grouped  = varMetricsSeasonal.groupby(level=0)
	#grouped.get_group(grouped.groups.keys()[0])['psfFluxMean'].values

	varMetricsFullSeasonal = grouped.apply(varF.ComputeVarFullBinned)
	#
	######################### SAVING OUTPUT        ######################### 
	#
	path = DirOut + 'FullSeasVar'+name
	varMetricsSeasonal.to_csv(path)
	print('Saving Full Seasonally binned LC statistics to %s'%path)
	   
	   






