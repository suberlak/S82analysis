# -*- coding: iso-8859-1 -*-
import pandas as pd 
import numpy as np 

# faint source treatment 
import faintFunctions as faintF 
import variabilityFunctions as varF

name = lToDoFilt[0]

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
    print('Okay, we dropped %d rows where psfFlux is NaN or inf'%np.sum(m))

#### check to make sure that there are no NaN psfFluxErr... 
m1  = np.isnan(fp_data['psfFluxErr'])  # True if NaN  
m2 =  ~np.isfinite(fp_data['psfFluxErr']) #  True if not finite  
m  = m1 | m2  # a logical or 
if np.sum(m) > 0 :  # only apply if there is anything to drop ... 
    fp_data.drop(m.index[m], inplace=True)
    print('Okay, we dropped %d rows where psfFluxErr is NaN or inf'%np.sum(m))
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

fp_data.ix[mask, 'faintMean'] = faintF.calculate_mean(fp_data['psfFlux'][mask].values,fp_data['psfFluxErr'][mask].values)
fp_data.ix[mask, 'faintMedian'] = faintF.calculate_median(fp_data['psfFlux'][mask].values,fp_data['psfFluxErr'][mask].values)
fp_data.ix[mask, 'faintTwoSigma'] = faintF.calculate_2sigma(fp_data['psfFlux'][mask].values,fp_data['psfFluxErr'][mask].values)
fp_data.ix[mask, 'faintRMS'] = faintF.calculate_rms(fp_data['psfFlux'][mask].values,fp_data['psfFluxErr'][mask].values)
## If one wants, at this stage we may replace all fluxes by Mean or Median
# I would keep both, and then use either, as appropriate.
# This replaces all rows where flag == 1
#rows = fp_data['flag'] == 1
#fp_data.ix[rows, 'psfFlux'] = fp_data.ix[rows, 'faintMean']
# 
######################### SAVING OUTPUT (unneccesarily takes time...)
# as a new csv....
#path = DirOut+'Proc_'+name
#print('Saving processed file to %s'% path)
#fp_data.to_csv(path)  


#
##########  STEP 2 : Derived Quantities ###########  
#

####  replace all psfFlux  where SNR < 2  with  faintMean  
rows = fp_data['flagFaint'] == 1
fp_data.ix[rows, 'psfFlux'] = fp_data.ix[rows, 'faintMean']
fp_data = fp_data.dropna(subset=['psfFlux', 'psfFluxErr'])

# group by objectId to calculate full LC variability characteristics 
grouped = fp_data.groupby('objectId')

# 
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



Last few done objects  : 
objectId=  216357512281071831
objectId=  216357512281071832
objectId=  216357512281071875
objectId=  216357512281071876
objectId=  216357512281071878
objectId=  216357512281071929
objectId=  216357512281071969
objectId=  216357512281071970
objectId=  216357512281071972

last id : 39238240 

The error that threw all up : 
/astro/apps6/anaconda2.0/bin/python/astro/apps6/anaconda2.0/bin/python: : /astro/apps6/anaconda2.0/bin/pythonsymbol lookup errorsymbol lookup error: : /astro/apps6/anaconda2.0/bin/python/astro/apps6/anaconda2.0/bin/python: : symbol lookup error/astro/apps6/anaconda2.0/bin/python: symbol lookup error: /astro/apps6/anaconda2.0/lib/python2.7/site-packages/numexpr/../../../libmkl_vml_avx.so: undefined symbol: mkl_serv_getenv

--> so it came to the end of my file, and yet it could not proceed... 
--> What I can do is to insert some print statements, and hope that this error is not persistent... 

--> if this error is persistent, contact PACS!  


Again, I got to the end of analysis, and when trying to save I get an erro 

/astro/apps6/anaconda2.0/bin/python/astro/apps6/anaconda2.0/bin/python/astro/apps6/anaconda2.0/bin/python: : : symbol lookup errorsymbol lookup errorsymbol lookup error: : : /astro/apps6/anaconda2.0/lib/python2.7/site-packages/numexpr/../../../libmkl_vml_avx.so/astro/apps6/anaconda2.0/lib/python2.7/site-packages/numexpr/../../../libmkl_vml_avx.so/astro/apps6/anaconda2.0/lib/python2.7/site-packages/numexpr/../../../libmkl_vml_avx.so: : : undefined symbol: mkl_serv_getenvundefined symbol: mkl_serv_getenvundefined symbol: mkl_serv_getenv


/astro/apps6/anaconda2.0/bin/python: symbol lookup error: /astro/apps6/anaconda2.0/lib/python2.7/site-packages/numexpr/../../../libmkl_vml_avx.so
