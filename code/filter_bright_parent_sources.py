# -*- coding: iso-8859-1 -*-
import pandas as pd
import numpy as np

##
# A code to find objects that have too bright parents 
# apparently the deblender does not work well if parent iPsfMag < 17 
##

Dir = '/astro/store/scratch/tmp/suberlak/S13Agg/repo_fls/'

# data processing center  (NCSA or IN2P3) 
data_source  = 'IN2P3'

# a file with primary sources, that has also parent sources listed 

if data_source == 'NCSA':
    primary_file = 'DeepSourceNCSA_i_lt235.csv.gz'
    not_primary_file = 'DeepSourceNCSA_i_lt300_narrow.csv.gz'
    too_bright_file = 'DeepSourceNCSA_i_lt_235_parent_lt_17_too_bright.csv'

if data_source == 'IN2P3':
    primary_file = 'DeepSourceIN2P3_i_lt235.csv.gz'
    not_primary_file = 'DeepSourceIN2P3_i_lt235_narrow_not_primary.csv.gz'
    too_bright_file = 'DeepSourceIN2P3_i_lt_235_parent_lt_17_too_bright.csv'

primary = pd.read_csv(Dir+primary_file, compression='gzip')
# IN2P3 has 4998901  rows  (5 mln) , each row is a unique deepSourceId 
# NCSA has 5474350 rows (5.4 mln) , each row is a unique deepSourceId 
not_primary = pd.read_csv(Dir+not_primary_file, compression='gzip')
# IN2P3 has 1882303 rows (1.8 mln), each row is a unique deepSourceId
# NCSA has 20258442 rows (20.2 mln) , because the file has all objects up to
# 30th magnitude...  3744255 of those (3.7 mln) have a flag detect_is_primary =0, i.e. 
# are not primaries... 

merged = pd.merge(primary, not_primary, how='left', left_on = 'parentDeepSourceId', 
                  right_on = 'deepSourceId', suffixes=('_primary','_parent')) 
# parent Ids are 
# np.unique(merged['parentDeepSourceId_primary'].values)

# apart from -1 , which means that there is no parent, there may be a single parent to various 
# primary detections.  (eg. parentDeepSourceId_primary = 1398579193184483  has three 
# deepSourceId_primary


# add info about parent magnitude to remove those primary sources that have a parent with
# iPsfMag < 17 

# len(merged) == len(french_primary)
# remove NaN rows (this would correspond to eg. primary sources which were not blended, i.e.
# have no parent , and so during merge, their parentDeepSourceId = -1  would match no rows 
# with the french_not_primary  
# 
mask_finite = np.isfinite(merged['psfMag_parent'].values)
# How many have a parent brighter than 17 mag ? 
mask_bright_parents = merged['psfMag_parent'].values[mask_finite] < 17
np.sum(mask_bright_parents)
# IN2P3 : 395970 ,  NCSA: 585920 -  that many rows, which corresponds to 
len(np.unique(merged['deepSourceId_primary'][mask_finite].values[mask_bright_parents]))
# the same number of deepSourceIds  that have too bright parent source to be reliably 
# handled by the deblender...  

# It's as much as   0.4 mln out of 4.9 mln  -  8 %  (IN2P3),  
# and 0.58 mln out of 5.4 mln - 10 % (NCSA)...

#  Save those objects (to be ignored in all future analysis as unreliable)
#  to a file "IN2P3_lt235_primaries_with_parents_lt_17.csv" 

merged_too_bright = merged[mask_finite][mask_bright_parents]
merged_too_bright.to_csv(Dir+too_bright_file)

merged_too_bright_narrow = merged_too_bright[['deepSourceId_primary','parentDeepSourceId_primary','psfMag_primary', 'psfMag_parent']]
merged_too_bright_narrow.to_csv(Dir+too_bright_file[:-4] + '_narrow.csv')