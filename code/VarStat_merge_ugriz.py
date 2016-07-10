# -*- coding: iso-8859-1 -*-
# Necessary imports 
import pandas as pd
import numpy as np 
import matplotlib.pyplot as plt 
from pandas import compat

#
# A program to merge ugriz variability statistics for NCSA and IN2P3, 
# that are a result of running LC_processing.py,  stored in 
# /astro/store/scratch/tmp/suberlak/s13_stripe82/forced_phot_lt_23/
#
# It combines the photometry and variability parameters on full LCs 
# with the E(B-V), correcting for extinction, and 
# adding the RA, DEC, extendedness information from DeepSource files,
# located in  
# /astro/store/scratch/tmp/suberlak/S13Agg/repo_fls/'
 
data_source  = 'NCSA'
print('Combining ugriz variability results for forced photometry lightcurves from %s'%data_source)
dir_info_files = '/astro/store/scratch/tmp/suberlak/S13Agg/repo_fls/'
dir_var_stats =  '/astro/store/scratch/tmp/suberlak/s13_stripe82/forced_phot_lt_23/'+data_source+'/Var/'
dir_save_merged= '/astro/users/suberlak/S13Agg_analysis/data_products/'

# Get E(B-V) 
if data_source == 'NCSA':
    ebv_file = 'ebv_NCSA_lt235.dat.gz'
else:
    ebv_file = 'ebv_IN2P3_lt230.dat.gz'

ebv = pd.read_table(dir_info_files+ebv_file, delimiter=' ', usecols=[0,1])
ebv.columns = ['objectId','ebv']

# read in the  file listing objects with bright parents, to ignore at the end of merger...
bright_parents_file = 'DeepSource'+data_source+'_i_lt_235_parent_lt_17_too_bright_narrow.csv'
objects_with_bright_parents = pd.read_csv(dir_info_files+bright_parents_file, usecols=[1,2,3,4])

def process_patch(patch='00_21', ebv = ebv, varPatchesDF = None, dir_var_stats=dir_var_stats, ebv_file =  ebv_file):
    
    print('\nCombining ugriz variability results from patch %s'%patch)
    
    # Directory storing the results of full LC statistics...
    
    
    # only read columns that we use to speed up execution.... 
    columns = ['objectId','N','chi2DOF','muFull','sigmaFull','psfMean']

    # Read in all filters per patch ... 
    varPatch = {}
    for filter in 'ugriz':
        File = 'Var'+filter+patch+'.csv'
        varPatch[filter] = pd.read_csv(dir_var_stats+File, usecols = columns)
    
    # Read in sigma and mu in magnitude units...
    varPatchMag = {}
    for filter in 'ugriz':
        File = 'VarSigMag'+filter+patch+'.csv'
        varPatchMag[filter] = pd.read_csv(dir_var_stats+File)

    # Check if each patch-filter file has exactly the same number of objects... 
    for filter in 'ugriz':
        print('Number of unique objectId in %s is %d'%(filter, len(np.unique(varPatch[filter]['objectId'].values))))

    # add prefix for each filter, apart from the objectId column  
    for filter in 'ugriz':
        varPatch[filter].columns = [filter+col  if col != 'objectId' else col for col in varPatch[filter]]
        varPatchMag[filter].columns = [filter+col  if col != 'objectId' else col for col in varPatchMag[filter]]
    
    # merge sigma results in magnitudes
    howmerge='inner' 
    varPatchMagug = pd.merge(varPatchMag['u'], varPatchMag['g'], how=howmerge, on='objectId', copy=True, indicator=False)
    varPatchMagugr = pd.merge(varPatchMagug, varPatchMag['r'], how=howmerge, on='objectId', copy=True, indicator=False)
    varPatchMagiz = pd.merge(varPatchMag['i'], varPatchMag['z'], how=howmerge, on='objectId', copy=True, indicator=False)
    varPatchMagugriz = pd.merge(varPatchMagugr, varPatchMagiz , how=howmerge, on='objectId', copy=True, indicator=False)

    # merge ugriz  
    howmerge='inner' # to avoid those objects which were in one filter but not in the other...
    varPatchug = pd.merge(varPatch['u'], varPatch['g'], how=howmerge, on='objectId', copy=True, indicator=False)
    varPatchugr = pd.merge(varPatchug, varPatch['r'], how=howmerge, on='objectId', copy=True, indicator=False)
    varPatchiz = pd.merge(varPatch['i'], varPatch['z'], how=howmerge, on='objectId', copy=True, indicator=False)
    varPatchugriz = pd.merge(varPatchugr, varPatchiz , how=howmerge, on='objectId', copy=True, indicator=False)

    # merge ugriz across bands with sigma in mag which was calculated later on ....    
    varPatchugriz = pd.merge(varPatchugriz, varPatchMagugriz, how=howmerge, on='objectId', copy=True, indicator=False)

    # check how many objects have ugriz photometry and extinction information 
    withEBV = np.sum(np.in1d(varPatchugriz['objectId'].values, ebv['objectId'].values))
    allOBJ = len(varPatchugriz['objectId'].values)
    print('Of all %d objects with ugriz info, %d have E(B-V) values from %s'%(allOBJ, withEBV, ebv_file))

    # Now this can be a left merge - I only want objects that can be extinction-corrected 
    varPatchAll =pd.merge(varPatchugriz, ebv, how='inner', on='objectId', copy=True, indicator=False)

    # Correct for extinction 
    A = [5.155, 3.793, 2.751, 2.086, 1.479]
    filters = 'ugriz'

    for i in range(len(A)):
        label = filters[i] + 'psfMean'
        varPatchAll[label+'_corr'] = varPatchAll[label] +  varPatchAll['ebv'] * A[i]

    # Drop unnecessary columns with uncorrected magnitudes.... 
    compat.PY3 = True
    varPatchSave = varPatchAll.drop(['u'+'psfMean'], axis=1)

    for filter in 'griz':
        varPatchSave = varPatchSave.drop(filter+'psfMean', axis=1)

    # add a column saying which patch an object comes from...
    varPatchSave['patch'] = patch

    
    if varPatchesDF is not None : 
        varPatchesDF = varPatchesDF.append(varPatchSave)
        print('Now we have %d objects total'%len(varPatchesDF))
        
    else : 
        varPatchesDF = varPatchSave
        
    return varPatchesDF


if data_source == 'NCSA':  # NCSA patches (11)
    patches = ['00_21', '22_43', '44_65','66_87', '88_109','110_131', '132_153', 
               '154_175',  '176_181', '365_387', '388_409']

if data_source == 'IN2P3': # IN2P3 patches (11)
    patches = ['155_176', '176_197','197_218', '218_239', '239_260', '260_281', 
               '281_302',  '302_323','323_344', '344_365', '365_386']

#  
# Run the first patch to start the storage DF 
varPatchesDF=  process_patch(patch=patches[0], ebv = ebv, varPatchesDF = None)

# Loop over the rest of the patches to append to that DF 
for patch in patches[1:]:
    varPatchesDF=  process_patch(patch=patch, ebv = ebv, varPatchesDF = varPatchesDF)
    
compat.PY3 = False

# Select out those objects that had parents brighter than iPsfMag  (uncorrected for extinction)
# < 17 mag , because for those objects the deblender was not working properly,
# and variability may be spurious qq

mask_bright_objects = np.in1d(varPatchesDF['objectId'].values, objects_with_bright_parents['deepSourceId_primary'].values)
print('For a total of %d objects'%len(varPatchesDF))
print('There are %d objects with parents brighter than 17 mag which we do not keep'%np.sum(mask_bright_objects))
name_ignored = 'test_mag_'+data_source+'_all_patches_ugriz_var_bright_parents_ignored.csv' 
print('Their info is stored in %s'%(dir_save_merged+name_ignored))

print('We keep the remainder of %d objects '%np.sum(~mask_bright_objects))
name_keep = 'test_mag_'+data_source+'_all_patches_ugriz_var.csv'
print('Saving them as %s '%(dir_save_merged+name_keep))

varPatchesDF_with_bright_parents = varPatchesDF[mask_bright_objects]
varPatchesDF_keep =   varPatchesDF[~mask_bright_objects]

# Save the combined files ... 
varPatchesDF_keep.to_csv(dir_save_merged+name_keep)
varPatchesDF_with_bright_parents.to_csv(dir_save_merged+name_ignored)


