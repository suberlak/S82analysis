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


###############################################
##### REDO SESAR PLOT  23  
############################################### 


#raMin = [0, 10,20,30,40,50,60, 280,290,300,310,320,330,340, 350]
#raMax = [10,20,30,40,50,60,280,290,300,310,320,330,340,350, 360]



#ct = 0
#for j in range(len(raMin)):
#    maskRa =  (raMin[j] < ra)  & (ra < raMax[j])
#    num = np.sum(maskRa.values)
#    print('For %d<ra<%d,  we have %d points'%(raMin[j],raMax[j], num))
#    ct += num  

# The outcome of above function : how many objects do we have in each bin...

# expected, because  ra ranges are
#  +320 < NCSA < +10   and    +5 < IN2P3 < +55

#For 50<ra<60,  we have 0 points
#For 60<ra<280,  we have 0 points
#For 280<ra<290,  we have 0 points
#For 290<ra<300,  we have 0 points
#For 300<ra<310,  we have 0 points

#For 310<ra<320,  we have 24970 points
#For 320<ra<330,  we have 1419807 points
#For 330<ra<340,  we have 1225997 points
#For 340<ra<350,  we have 1125263 points
#For 350<ra<360,  we have 1046624 points
#For 0<ra<10,  we have 1049393 points

#For 10<ra<20,  we have 0 points
#For 20<ra<30,  we have 0 points
#For 30<ra<40,  we have 0 points
#For 40<ra<50,  we have 0 points
# ttl: 0<ra<360 : 5892054  points ( == len(ra.values)) 


#   DO THE 3 PANEL PLOTS 

raMin = [350, 330, 310]
raMax = [360, 340, 320]
fig, ax = plt.subplots(1,3, figsize=(18,6))
#axs = ax.ravel()
# Select a range of gi,  i  for the histogram.... for ra 350-360 
ylim = [15,23.5]
xlim = [0,1.3]
nbinsArr = [100,100,50]

print('Number of point sources is %d'%np.sum(med_cat_merged_all['extendedness'] == 0))
# 1908043
print('Number of extended sources is %d'%np.sum(med_cat_merged_all['extendedness'] == 1))
# 


for j in range(len(raMin)):
    print('j=%d'%j)
    nbins = nbinsArr[j]
    maskRa =  (raMin[j] < med_cat_merged_all['ra'])  & (med_cat_merged_all['ra'] < raMax[j])
    maskExtend = med_cat_merged_all['extendedness'] == 0 
    maskTot = maskRa * maskExtend
    gi = med_cat_merged_all['median_mag_g_corr'][maskTot] - med_cat_merged_all['median_mag_i_corr'][maskTot]
    i = med_cat_merged_all['median_mag_i_corr'][maskTot]
    maskCol = (gi < xlim[1])&(gi>xlim[0])&(i<ylim[1])&(i>ylim[0])
    H, xedges, yedges = np.histogram2d(gi[maskCol],i[maskCol],bins=nbins)
    # H needs to be rotated and flipped
    H = np.rot90(H)
    H = np.flipud(H)
    # Mask zeros
    Hmasked = np.ma.masked_where(H==0,H) # Mask pixels with a value of zeropltssssss
    # Plot 2D histogram using pcolor
     
    ax[j].pcolormesh(xedges,yedges,np.log10(Hmasked), cmap='jet')
    ax[j].set_xlabel('g-i', fontsize=15)
    ax[j].set_ylabel('i', fontsize=15)
    ax[j].set_xlim(xlim[0],xlim[1])
    ax[j].set_title(str(raMin[j])+'< RA <' + str(raMax[j]))
    ax[j].invert_yaxis()
    #plt.gca().invert_yaxis()
    ax[j].tick_params(axis = 'both', labelsize=15) 

#plt.tick_params(axis='both', which='major', labelsize=15)
plt.tight_layout()
#cbar = plt.colorbar()
#cbar.ax.set_ylabel(r'$\log_{10}{Counts}$')
plt.savefig('Fig_g-i_vs_i_ra_'+str(min(raMin))+'-'+str(max(raMax))+'hist_n_'+str(nbins)+'row_ext_0.png')


# DO 3 SEPARATE PLOTS....

raMin = [350, 330, 310]
raMax = [360, 340, 320]

# Select a range of gi,  i  for the histogram.... for ra 350-360 
ylim = [15,23.5]
xlim = [0,1.3]
nbins = 100

for j in range(len(raMin)):
    print('j=%d'%j)
    plt.clf()
    fig2 = plt.figure()
    maskRa =  (raMin[j] < ra)  & (ra < raMax[j])
    gi = med_cat_merged_all['median_mag_g_corr'][maskRa] - med_cat_merged_all['median_mag_i_corr'][maskRa]
    i = med_cat_merged_all['median_mag_i_corr'][maskRa]
    maskCol = (gi < xlim[1])&(gi>xlim[0])&(i<ylim[1])&(i>ylim[0])
    H, xedges, yedges = np.histogram2d(gi[maskCol],i[maskCol],bins=nbins)
    # H needs to be rotated and flipped
    H = np.rot90(H)
    H = np.flipud(H)
    # Mask zeros
    Hmasked = np.ma.masked_where(H==0,H) # Mask pixels with a value of zero
    # Plot 2D histogram using pcolormesh
    plt.pcolormesh(xedges,yedges,np.log10(Hmasked), cmap='jet')
    plt.xlabel('g-i', fontsize=15)
    plt.ylabel('i', fontsize=15)
    plt.xlim(xlim[0],xlim[1])
    plt.title(str(raMin[j])+'< RA <' + str(raMax[j]))
    plt.gca().invert_yaxis()
    plt.tick_params(axis='both', which='major', labelsize=15)
    plt.tight_layout()
    cbar = plt.colorbar()
    cbar.ax.set_ylabel(r'$\log_{10}{Counts}$')
    plt.savefig('Fig_g-i_vs_i_ra_'+str(raMin[j])+'-'+str(raMax[j])+'hist_n_'+str(nbins)+'_ext_0.png')




plt.scatter(gi,i, s=0.001)
plt.xlim(0,2)
plt.ylim(14,25)
plt.xlabel('g-i', fontsize=15)
plt.ylabel('i', fontsize=15)
plt.gca().invert_yaxis()
plt.tick_params(axis='both', which='major', labelsize=15)
plt.tight_layout()
plt.savefig('Fig_g-i_vs_i_ra_'+str(raMin[j])+'-'+str(raMax[j])+'hist.png')
plt.show()

#######################################################
###   PHOTOMETRIC DATA : RMS G-BAND VS MEAN U-G COLOR #
#######################################################

# Photometry in each small square patch was saved to  
#  /rawDataFP/

# When it was  sorted according to objID, the result was saved to 
#  /rawDataFPSplit/ 
#  g176_181  means   g filter,  from patch 176 to 181.  

# Use this for LIGHTCURVES   because they are here 
# 


fp_data = pd.read_csv(DirIn+'rawDataFPSplit/g176_181.csv')

len(np.unique(fp_data['objectId'].values))
len(med_cat_merged_all['objectId'].values)

