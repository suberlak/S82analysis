# -*- coding: iso-8859-1 -*-
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
# Note:  
# %cpaste is sensitive to tabs : 
# make sure that in whatever editor you do tabs
# it actually puts in 4 space

# make all necessary imports....
import os 
import numpy as np
#
# to avoid  http://stackoverflow.com/questions/20625582/how-to-deal-with-this-pandas-warning

# for all imports of my functions, 
# make python aware of my packages...
import sys
sys.path.insert(0, '/astro/users/suberlak/S13Agg_analysis/packages/')

# processing
import processPatch as procP

#####
####  PROCESSING FILES IN   /astro/store/pogo4/s13_stripe82/forced_phot_lt_23/
####  


DirIn = '/astro/store/pogo4/s13_stripe82/forced_phot_lt_23/NCSA/'
DirOut = '/astro/store/scratch/tmp/suberlak/s13_stripe82/forced_phot_lt_23/NCSA/'


#lProc = []
#lProc += [each for each in os.listdir(DirOut) if each.endswith('.csv')]

#lProc = [name[5:-5] for name in lProc]
#lToDo = os.listdir(DirIn)
#lToDoComp = [name[:-3] for name in lToDo]
#if len(lProc) > 0 : 
#    maskDone = np.in1d(lToDoComp,lProc)
#    lToDoFilt =  np.array(lToDoComp)[~maskDone]
#else:
#    lToDoFilt = lToDoComp

# Get only those done first, because then you can at least get working on colors, etc ! 

n = '66_87'  #'44_65' #'22_43'

for filter in 'ugriz':
    lToDoFilt.append(filter + n + '.csv')
    
#lToDoFilt = ['g00_21.csv', 'r00_21.csv','i00_21.csv','u00_21.csv','z00_21.csv']

for name in lToDoFilt:
    procP.process_patch(name, DirIn, DirOut)
 







