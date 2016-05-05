import pandas as pd
import numpy as np 
import matplotlib.pyplot as plt 
from matplotlib.pyplot import cm 
from matplotlib.colors import ListedColormap

# DIR with all DATA 
DirIn =  '/astro/users/suberlak/S13Agg/'

######################################################
##############  ONLY NCSA DATA !!!!! #################
######################################################



# DETECTION DATA 
coadd_data = pd.read_csv(DirIn+'iCoaddPhotometryAll.csv')

# MAKE AN EXTENDEDNESS COLUMN... (same threshold as I found 
# for DeepSourceNCSA_i_lt300.csv , 
# by merging coadd_data with deep_src,  and finding what is 
# the maximum value of iPsf - iModel   for extentedness=0 )

# Initialize all with zeros, and those that are bigger than the threshold, become 1 
threshold = 0.084645219604293942
coadd_data['extendedness'] = 0
coadd_data.ix[coadd_data['iPsfMag'] - coadd_data['iModelMag'] > threshold, 'extendedness'] = 1 


###############################################
####   HOW EXTENDEDNESS PARAMETER WORKS...
###############################################

 
# EXTENDEDNESS : a derivative of DETECTION DATA for sources less than 30th mag 
# Instead of reading it in , may as well just calculate the extendedness
# as  a binary depending on whether  iPsfMag - iModelMag > or < threshold...
# In DeepSourceNCSA   the threshold used was 0.084645 

deep_src = pd.read_csv(DirIn+'DeepSourceNCSA_i_lt300.csv')


print(' How many sources do we have in coadd phot  ? ')
print(len(np.unique(coadd_data['deepSourceId'].values)))
# 16520093    (16 mln  srcs) 

nmatch = np.sum(np.in1d(coadd_data['deepSourceId'].values, deep_src['deepSourceId'].values))
print('Matched sources with extendedness info: %d' % nmatch)
# 16514187     (almost 16 mln : 5906 orphan sources, that are not here... )

# Merge the iPsf - iModel    and  extendedness information for the plot .... 
coadd_merged = pd.merge(coadd_data, deep_src, how='left', on='deepSourceId', suffixes=('_x', '_y'), copy=True, indicator=False)

# Maximum value of the  iPsf - iModel   for compact sources   is the cutoff value 
msk = coadd_merged['extendedness'] == 0
diff = iPsf[msk]-iModel[msk]
print('The threshold  iPs - iModel for extendedness is %f' % max(diff.values))


# I use  coadd_merged  to show that the extendedness parameter is not too bad ... 
# I make plot of iPsfMag - iModelMag  vs iModelMag  , and then plot histograms of 
# iPsfMag - iModelMag in several iModMag bins  

cmap = ListedColormap(['red', 'green'])


iPsf = coadd_merged['iPsfMag']
iModel = coadd_merged['iModelMag']
color = coadd_merged['extendedness']
# limit before plotting 
xlim = [15,25]
ylim = [-0.5,2.5]
mask_x = (iModel > xlim[0]) & (iModel < xlim[1])
mask_y = ((iPsf-iModel) < ylim[1]) & ((iPsf-iModel) > ylim[0])
mask = mask_x * mask_y
plt.clf()
fig, ax = plt.subplots(1,1, figsize=(8,6))
ax.set_ylabel('iPsfMag-iModelMag', fontsize=15)
ax.set_xlabel('iModelMag', fontsize=15)
ax.tick_params(axis='both', which='major', labelsize=15)
ax.set_xlim(xlim[0], xlim[1])
ax.set_ylim(ylim[0], ylim[1])
ax.scatter(iModel[mask],  iPsf[mask]-iModel[mask], s=0.01, c=color[mask], cmap=cmap,lw = 0)
fig.tight_layout()
plt.savefig('Extendedness_coadd_data_'+str(len(iPsf))+'_srcs_lim.png')
#plt.show()


# Divide into 5 bins by iModelMag,  and plot histograms... 


n = 5
color=iter(cm.rainbow(np.linspace(0,1,n)))

ylim = [-0.5,2.5]
mask_y = ((iPsf-iModel) < ylim[1]) & ((iPsf-iModel) > ylim[0])

bin_start = 22
bin_end = 22.5
nbins=50
fig,ax = plt.subplots(2,3,  figsize=(12,8),sharex = False)
axs = np.ravel(ax)
for i in range(5):
    
    print('Range iModelMag %.1f to %.1f' %(bin_start, bin_end))
    
    mask_bin = (iModel > (bin_start)) & (iModel < (bin_end))
    mask = mask_bin * mask_y 
    print np.sum(mask)
    
    hist, bin_edges = np.histogram(iPsf[mask]-iModel[mask], bins=nbins, density=False)
    bin_cen = (bin_edges[:-1] + bin_edges[1:])/2
    axs[i].plot(bin_cen, hist, color = 'red', ls='steps', c =next(color), lw=2, label='bins:'+str(bin_start)+'-'+str(bin_end))
    axs[i].legend(fontsize=12)
    axs[i].tick_params(axis='both', which='major', labelsize=12)
    bin_start += 0.5 
    bin_end += 0.5 
axs[i+1].axis('off')
fig.text(0.5, 0.04, 'iPsfMag-iModelMag', ha='center', va='center',fontsize=20)
fig.text(0.03, 0.5, 'Counts', ha='center', va='center', rotation='vertical',fontsize=20)

fig.tight_layout()
fig.subplots_adjust(wspace=0.36, hspace=0.26, left=0.12, right=0.94, bottom=0.10, top=0.95)
plt.savefig('Extendedness_coadd_histograms.png')



# Convinced that extendedness indeed is a good parameter, 
# I grab only extendedness and location (ra,dec)  from deep_src  to merge with the photometric information below...
deep_src_ext = deep_src[['deepSourceId', 'extendedness', 'ra', 'decl']]


###
###   COLORS  AND  EXTINCTION CORRECTION  
###

# COADD PHOTOMETRY 
# Use ONLY FOR E(B-V)  correction   --> only grab the objectId and ebv 
med_coadd_phot = pd.read_csv(DirIn+'medianPhotometry.csv', usecols=['objectId','ebv'])
# There are  12373162  : 12.3 mln values of objectId 

# From this file,  grab the objectId and ebv,  and join  on objectId (which is unique to both 
# medianPhotometry.csv and ugrizMetrics.csv)
#med_coadd_ebv = med_coadd_phot[['objectId', 'ebv']]


# CATALOG PHOTOMETRY 	
# Use for COLORS   !! 
med_cat_phot = pd.read_csv(DirIn+'lightcurveMetrics/metricsFiles/ugrizMetrics.csv')
# There are 5892054  : 5.8 mln values of objectId here 

# MERGE CATALOG and COADD for ebv 
# Join with med_coadd_phot  from above for ebv info ... Check that all objectIds are unique, and exist only once in each 
# dataset !  
med_cat_merged = pd.merge(med_cat_phot, med_coadd_ebv, how='left', on='objectId', suffixes=('_x', '_y'), copy=True, indicator=False)
# Since it was a left-join ,  it has the same objects as in med_cat_phot : 5892054 : 5.8 mln values of objectId

# CORRECT FOR EXTINCTION using the ebv 

A = [5.155, 3.793, 2.751, 2.086, 1.479]
filters = ['u','g','r','i','z']

for i in range(len(filters)):
    label = 'median_mag_'+filters[i]
    med_cat_merged[label+'_corr'] = med_cat_merged[label] +  med_cat_merged['ebv'] * A[i]

# Merge with EXTENDEDNESS information from 

# How many elements are shared between DeepSrcIDs and ObjectIDs?  
nshared = np.sum(np.in1d(med_coadd_ebv['objectId'].values, deep_src_ext['deepSourceId'].values)) 
print('Shared IDs between  DeepSourceNCSA_i_lt300.csv (extendedness)  and ugrizMetrics.csv')
print(' %f of %f ' % (nshared, len(med_coadd_ebv['objectId'].values)))

# 11998309 : 11.9 mln  objectIDs out of 12373162: 12.3 mln  are shared by the deepSourceIDs ( of which there is 20258442 : 20.2 mln values) 

# Only add extended

med_cat_merged_all = pd.merge(med_cat_merged, deep_src_ext, how='left', left_on='objectId', right_on='deepSourceId' ,suffixes=('_x', '_y'), copy=True, indicator=False)
# len(med_cat_merged_all['objectId'].values) = 5892054
# Since it was a left-join, it has the same objects as med_cat_merged: 5892054 : 5.8 mln values of objectId 

