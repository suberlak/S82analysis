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

site = 'IN2P3'

DirIn = '/astro/store/pogo4/s13_stripe82/forced_phot_lt_23/'+site+'/'
DirOut = '/astro/store/scratch/tmp/suberlak/s13_stripe82/forced_phot_lt_23/'+site+'/'


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

lToDoFilt= []
#patches = ['00_21','22_43','44_65', '66_87' ,'88_109','110_131', '132_153', 
#        '154_175',  '176_181', '365_387', '388_409']  # NCSA 
patches = ['155_176', '176_197','197_218', '218_239', '239_260', '260_281', '281_302', 
         '302_323','323_344', '344_365', '365_386']  # IN2P3
# 
for patch in patches  :
    for filter in 'ugriz':
        lToDoFilt.append(filter + patch + '.csv')
    
#lToDoFilt = ['g00_21.csv', 'r00_21.csv','i00_21.csv','u00_21.csv','z00_21.csv']

# Processing full LC, but only calculating mu, sigma based on magnitudes 
# for now, not calculating seasonal quantities... 

for name in lToDoFilt[16:]:
    procP.process_patch_minimal(name, DirIn, DirOut)

# Full LC processing, where we calculate full-LC statiscics, seasonal, and seasonally-binned
# but both sigma and mu are based on fluxes.. 

#for name in lToDoFilt:
#    procP.process_patch(name, DirIn, DirOut)
 







