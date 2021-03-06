# -*- coding: iso-8859-1 -*-
import pandas as pd
import numpy as np
import os 
import sys
sys.path.insert(0, '/astro/users/suberlak/S13Agg_analysis/packages/')
pd.options.mode.chained_assignment = None

# God, Relieve me of the bondage of self,
# that I may better do Thy will.
# Take away my difficulties,
# that victory over them may bear witness
# to those I would help of Thy Power,
# Thy Love, and Thy Way of life.
# May I do Thy will always!

# faint source treatment 
import faintFunctions as faintF 

# variability 
print('I see that line...') 
import variabilityFunctions as varF
reload(varF)

def process_patch(name, DirIn, DirOut, limitNrows=None):
    print('Processing filter_patch file %s' % name)
    #DirIn = '/astro/store/pogo4/s13_stripe82/forced_phot_lt_23/NCSA/'
    #DirOut = '/astro/store/scratch/tmp/suberlak/s13_stripe82/forced_phot_lt_23/NCSA/'

    # read in the raw lightcurve... 
    if limitNrows is not None:
        fp_data = pd.read_csv(DirIn+name+'.gz', compression='gzip',  
                     usecols=['objectId', 'mjd', 'psfFlux', 'psfFluxErr'], nrows=limitNrows)
    else : 
        fp_data = pd.read_csv(DirIn+name+'.gz', compression='gzip',  
                     usecols=['objectId', 'mjd', 'psfFlux', 'psfFluxErr'])
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

    #### check to make sure that there are no NaN or 0 psfFluxErr... 
    m1  = np.isnan(fp_data['psfFluxErr'])  # True if NaN  
    m2 =  ~np.isfinite(fp_data['psfFluxErr']) #  True if not finite
    m3 =   fp_data['psfFluxErr'].values == 0 # True if Err == 0  (IN2P3 problem...)
    m  = m1 | m2 | m3  # a logical or 
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
    print('Faint points treatment...')
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

    # Save the diagnostics..
    print('Saving the diagnostics...')
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
    s = '    '+ DirOut +'Var/'+ 'Var'+name+ '\n'  
    file.write(s)
    file.close()   

    #
    ##########  STEP 2 : Derived Quantities ###########  
    #

    ####  replace all psfFlux  where SNR < 2  with  faintMean  
    rows = fp_data['flagFaint'] == 1
    fp_data.ix[rows, 'psfFlux'] = fp_data.ix[rows, 'faintMean']

    # group by objectId to calculate full LC variability characteristics 
    print('Calculating the full LC statistics...')
    grouped = fp_data.groupby('objectId')

    #
    #  Calculate low-order statistical properties of full lightcurves:
    # 
    #  psfFluxMean, psfFluxMeanErr, psfFluxMedian, psfFluxMedianErr
    #  psfFluxStDev, psfFluxSigG , psfFluxSkew, avgMJD, rangeMJD   
    #  chi2DOF , chi2R,  sigmaFull, muFull, flagLtTenPts
    # 
    
    varMetricsFull = grouped.apply(varF.computeVarMetrics)
    print('Calculating metrics for full lightcurves is finished')

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
    print('Calculating magnitudes from fluxes is finished')
    #
    ######################### SAVING OUTPUT        ######################### 
    # 
    path = DirOut +'Var/'+ 'Var'+name
    print('Saving varMetricsFull...')
    varMetricsFull.to_csv(path)
    print('Saved Full, unbinned  LC statistics to %s'%path)

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
    #m1 = (varMetricsFull['sigmaFull'] > 0)  # 234883
    #m2 =  (varMetricsFull['chi2DOF'] > 1).values   # 346886
    #m3 = (varMetricsFull['chi2R'] > 1).values      # 172602
    #m= np.ma.mask_or(m3, np.ma.mask_or(m1,m2))
    #print('Out of %d objects,  %d fulfill  sigma>0 or chi2R>1 or chi2DOF>1' % (len(m), np.sum(m)))

    # The  '3 sigma' definition (more strict) ...
    m1 = (varMetricsFull['sigmaFull'] > 0)
    N = varMetricsFull['N'] #  number of pts per lightcurve 
    m4 = varMetricsFull['chi2DOF'] > (1 + 3.0 * np.sqrt(2 / N )) ##  170517
    m5 = varMetricsFull['chi2R'] > (1 + 3.0 * np.sqrt(2 / N))   ## 37464
    #np.sum(m4|m5)  # one or the other   170534#
    #np.sum(m1|m4|m5) # 265981   sigma OR  chi2R  OR chi2DOF 
    #np.sum(m1*(m4|m5)) # 139436     sigma AND  (chi2R OR chi2DOF)

    # save the diagnostics about variability to a file...
    diagFile = path+'.diag'
    file = open(diagFile, "w")
    file.write('Input :  \n')
    s = '    '+ DirIn + '\n' 
    file.write(s)
    s = '    '+ name + '\n\n'
    file.write(s)
    s = "chi2DOF is  varMetricsFull['chi2DOF'] > (1 + 3.0 * np.sqrt(2 / N )) \n" 
    file.write(s)
    s = "chi2R is  varMetricsFull['chi2R'] > (1 + 3.0 * np.sqrt(2 / N )) \n " 
    file.write(s)
    m = m1|m4|m5
    s = 'There are '+str(np.sum(m)) + ' points out of ' + str(len(m)) + \
        ' that fulfill sigma OR  chi2R  OR chi2DOF  \n '
    file.write(s)
    m = m1*(m4|m5)
    s = 'There are '+str(np.sum(m)) + ' points out of ' + str(len(m)) + \
        ' that fulfill sigma  AND  (chi2R OR chi2DOF)  \n '
    file.write(s)
    file.write('\n')
    file.write('Output :  \n')
    s = '    '+ DirOut + '\n' 
    file.write(s)
    s = '    '+ 'Var'+name + '\n\n'  
    file.write(s)
    file.close()   

    #
    ################ SEASONAL METRICS #################
    #
    #m = m1*(m4|m5) (already defined when saving the var stats above) 

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
    print('Calculating Seasonal metrics for %s ...'%name)
    varMetricsSeasonal = grouped.apply(varF.computeVarMetrics)
    #
    ######################### SAVING OUTPUT        ######################### 
    #
    path = DirOut + 'Var/'+ 'SeasVar'+name
    varMetricsSeasonal.to_csv(path)
    print('Saving Seasonal statistics to %s'%path)

    #
    ################ SEASONALLY-BINNED LIGHTCURVE METRICS #################
    #
    # Calculate binned lightcurve metrics (binned by seasons)
    grouped  = varMetricsSeasonal.groupby(level=0)
    #grouped.get_group(grouped.groups.keys()[0])['psfFluxMean'].values
    print('Calculating Seasonally Binned LC statistics')
    varMetricsFullSeasonal = grouped.apply(varF.ComputeVarFullBinned)
    #
    ######################### SAVING OUTPUT        ######################### 
    #
    path = DirOut + 'Var/'+ 'FullSeasVar'+name
    varMetricsFullSeasonal.to_csv(path)
    print('Saving Full Seasonally binned LC statistics to %s'%path)



def process_patch_minimal(name, DirIn, DirOut, limitNrows=None):
    ''' A variation of process_patch : ONLY calculate sigma, mu , 
    based on psfMag and psfMagErr,  i.e.  psfFlux, psfFluxErr 
    after faint point treatment, translated to magnitudes  
    '''
    print('Processing filter_patch file %s' % name)
    #DirIn = '/astro/store/pogo4/s13_stripe82/forced_phot_lt_23/NCSA/'
    #DirOut = '/astro/store/scratch/tmp/suberlak/s13_stripe82/forced_phot_lt_23/NCSA/'

    # read in the raw lightcurve... 
    if limitNrows is not None:
        fp_data = pd.read_csv(DirIn+name+'.gz', compression='gzip',  
                     usecols=['objectId', 'mjd', 'psfFlux', 'psfFluxErr'], nrows=limitNrows)
    else : 
        fp_data = pd.read_csv(DirIn+name+'.gz', compression='gzip',  
                     usecols=['objectId', 'mjd', 'psfFlux', 'psfFluxErr'])
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

    #### check to make sure that there are no NaN or 0 psfFluxErr... 
    m1  = np.isnan(fp_data['psfFluxErr'])  # True if NaN  
    m2 =  ~np.isfinite(fp_data['psfFluxErr']) #  True if not finite
    m3 =   fp_data['psfFluxErr'].values == 0 # True if Err == 0  (IN2P3 problem...)
    m  = m1 | m2 | m3  # a logical or 
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
    print('Faint points treatment...')
    fp_data.ix[mask, 'faintMean'] = faintF.calculate_mean(fp_data['psfFlux'][mask].values,fp_data['psfFluxErr'][mask].values)
    fp_data.ix[mask, 'faintMedian'] = faintF.calculate_median(fp_data['psfFlux'][mask].values,fp_data['psfFluxErr'][mask].values)
    fp_data.ix[mask, 'faintTwoSigma'] = faintF.calculate_2sigma(fp_data['psfFlux'][mask].values,fp_data['psfFluxErr'][mask].values)
    fp_data.ix[mask, 'faintRMS'] = faintF.calculate_rms(fp_data['psfFlux'][mask].values,fp_data['psfFluxErr'][mask].values)
    
    ##########  STEP 2 : Derived Quantities ###########  
    #

    ####  replace all psfFlux  where SNR < 2  with  faintMean  
    rows = fp_data['flagFaint'] == 1
    fp_data.ix[rows, 'psfFlux'] = fp_data.ix[rows, 'faintMean']
    
    def flux2absigma(flux, fluxsigma):
      """Compute AB mag sigma given flux and flux sigma"""
      FIVE_OVER_2LOG10 = 1.085736204758129569
      return FIVE_OVER_2LOG10 * fluxsigma / flux;


    def flux2ab(flux):
      """Compute AB mag given flux"""
      return -2.5 * np.log10(flux) - 48.6;
    
    # Calculate psfMag , psfMagErr from psfFlux,  psfFluxErr

    fp_data['psfMag'] = flux2ab(fp_data['psfFlux'])
    fp_data['psfMagErr'] = flux2absigma(fp_data['psfFlux'],fp_data['psfFluxErr'])


    # group by objectId to calculate full LC variability characteristics 
    print('Calculating the full LC statistics...')
    grouped = fp_data.groupby('objectId')

    # Calculate only sigmaFullMag  muFullMag   for full lightcurves from 
    # magnitudes     
    varMetricsFull = grouped.apply(varF.computeVarMetricsMagOnly)
    print('Calculating metrics for full lightcurves is finished')

   
    path = DirOut +'Var/'+ 'VarSigMag'+name
    print('Saving varMetricsFull...')
    varMetricsFull.to_csv(path)
    print('Saved Full, unbinned  LC statistics to %s'%path)





