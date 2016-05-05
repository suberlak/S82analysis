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

# faint source treatment 
from scipy.stats import norm
from scipy.special import erf

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

	fp_data = pd.read_csv(DirIn+name+'.gz', compression='gzip', usecols=['objectId', 'mjd', 'psfFlux', 'psfFluxErr'])

	# make a new column, fill with 0's
	fp_data['flag'] = 0

	# mask those rows that correspond to SNR < 2
	mask = (fp_data['psfFlux'].values / fp_data['psfFluxErr'].values) < 2

	# print info how many points are affected
	print('There are %d points of %d that have SNR<2' %(np.sum(mask),len(mask)))

	# set flag at those rows to 1
	fp_data.ix[mask, 'flag'] = 1

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

	# save as a new csv....
	path = DirOut+'Proc_'+name
	print('Saving processed file to %s'% path)
	fp_data.to_csv(path)  

	# Save the diagnostics..
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
	


