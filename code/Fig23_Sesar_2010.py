# -*- coding: iso-8859-1 -*-
# make all necessary imports....
import pandas as pd
import numpy as np 
import matplotlib.pyplot as plt 

####################
# META META META   #
####################
# Code here uses lightcurve - derived as well as detection data. 
# Detection data from coadded images iCoaddPhotometryAll.csv 
# provides ra dec location of each object, as well as iPsfFlux, iModelFlux
# 
# Lightcurve-derived data that Yusra provided (ugrizMetrics.csv) is based 
# on the median of all epochs of raw lightcurves (not corrected for faint 
# flux, etc.) .  Some of that is almost completely obsolete, made by the sql database
# for which I don't even have code (it never ended up in the github repo 
# because median over single-epoch detections is biased). I use 
# medianPhotometry.csv because it has E(B-V) data,  but that gma
####################
# A MINIMUM NEEDED TO MAKE A  PLOT LIKE THAT OF SESAR+2010, FIG.23 

# DIR with all DATA 
DirIn =  '/astro/users/suberlak/S13Agg/'

#
#
# COADD PHOTOMETRY - ONLY FOR E(B-V)  correction  
med_coadd_ebv = pd.read_csv(DirIn+'medianPhotometry.csv', usecols=['objectId','ebv'])
# There are  12373162  : 12.3 mln values of objectId 

# CATALOG PHOTOMETRY - USE FOR COLORS 
med_cat_phot = pd.read_csv(DirIn+'lightcurveMetrics/metricsFiles/ugrizMetrics.csv')
# There are 5892054  : 5.8 mln values of objectId here 

# DETECTION DATA  - FOR PSF and MODEL magnitudes 
coadd_det_mag = pd.read_csv(DirIn+'iCoaddPhotometryAll.csv', usecols=['deepSourceId','parentDeepSourceId','iPsfMag','iPsfMagSigma','iModelMag','iModelMagSigma'])




# MERGE  COADD and CATALOG for ebv 
med_cat_coadd = pd.merge(med_cat_phot, med_coadd_ebv, how='left', on='objectId', suffixes=('_x', '_y'), copy=True, indicator=False)
# Since it was a left-join ,  it has the same objects as in med_cat_phot : 5892054 : 5.8 mln values of objectId

# MERGE in DETECTION for iPdf and iModel mags 
med_cat_merged = pd.merge(med_cat_coadd, coadd_det_mag, how='left', left_on ='objectId', right_on='deepSourceId', suffixes=('_x', '_y'), copy=True, indicator=False)


# CORRECT MEDIAN UGRIZ BANDS FOR EXTINCTION using the ebv 
A = [5.155, 3.793, 2.751, 2.086, 1.479]
filters = ['u','g','r','i','z']

for i in range(len(filters)):
    label = 'median_mag_'+filters[i]
    med_cat_merged[label+'_corr'] = med_cat_merged[label] +  med_cat_merged['ebv'] * A[i]


# CORRECT DETECTION i-MAG  FOR EXTINCTION
labels = ['Psf', 'Model']
for label in labels:
    med_cat_merged['i'+label+'Mag'+'_corr'] = med_cat_merged['i'+label+'Mag']+  med_cat_merged['ebv'] * 2.086



# MERGE with EXTENDEDNESS information

# EXTENDEDNESS INFO , take only necessary columns 
deep_src_ext =  pd.read_csv(DirIn+'DeepSourceNCSA_i_lt300.csv', usecols=['deepSourceId', 'extendedness', 'ra', 'decl'])

# ADD EXTENDEDNESS 
med_cat_merged_all = pd.merge(med_cat_merged, deep_src_ext, how='left', left_on='objectId', right_on=s ,suffixes=('_x', '_y'), copy=True, indicator=False)
# len(med_cat_merged_all['objectId'].values) = 5892054
# Since it was a left-join, it has the same objects as med_cat_merged: 5892054 : 5.8 mln values of objectId 


# How many elements are shared between DeepSrcIDs and ObjectIDs?  
nshared = np.sum(np.in1d(med_coadd_ebv['objectId'].values, deep_src_ext['deepSourceId'].values)) 
print('Shared IDs between  DeepSourceNCSA_i_lt300.csv (extendedness)  and ugrizMetrics.csv')
print(' %f of %f ' % (nshared, len(med_coadd_ebv['objectId'].values)))
# 11998309 : 11.9 mln  objectIDs out of 12373162: 12.3 mln  are shared by the deepSourceIDs ( of which there is 20258442 : 20.2 mln values) 

print('Number of point sources is %d'%np.sum(med_cat_merged_all['extendedness'] == 0))
# 1908043
print('Number of extended sources is %d'%np.sum(med_cat_merged_all['extendedness'] == 1))
# 3984011  


# How many points we have from NCSA ...

#For 310<ra<320,  we have 24970 points
#For 320<ra<330,  we have 1419807 points
#For 330<ra<340,  we have 1225997 points
#For 340<ra<350,  we have 1125263 points
#For 350<ra<360,  we have 1046624 points
#For 0<ra<10,  we have 1049393 points


# https://github.com/olgabot/prettyplotlib/wiki/scatter-and-pcolormesh:-Motivating-examples 
#  http://matplotlib.org/basemap/users/examples.html 
# http://pyhogs.github.io/colormap-examples.html 

# COMPARE  iModelMag and iPsfMag to median_mag_i ... 
iMed = med_cat_merged_all['median_mag_i_corr']
iPsf = med_cat_merged_all['iPsfMag_corr']
iModel =med_cat_merged_all['iModelMag_corr'] 

nbins=100

ylim = [12,26]
xlim = [12,26]

fig, axs = plt.subplots(2,2, figsize=(12,12))
#ax = np.ravel(axs)

for i in range(2):
    # left : ext=0,  right : ext=1 
    ext=i # choose extendedness 

    # top plot : iMed vs iModel
    x = iMed
    y = iModel
    maskLim = (x < xlim[1])&(x>xlim[0])&(y<ylim[1])&(y>ylim[0])
    maskExt = med_cat_merged_all['extendedness'] == ext
    maskCol = maskLim*maskExt
    H, xedges, yedges = np.histogram2d(x[maskCol],y[maskCol], bins=nbins)
    # H needs to be rotated and flipped
    H = np.rot90(H)
    H = np.flipud(H)
    # Mask zeros
    Hmasked = np.ma.masked_where(H==0,H) # Mask pixels with a value of zero
    # Plot 2D histogram using pcolor
    axs[0,i].plot(xlim, ylim, 'k-', alpha=0.75, zorder=0)
    axs[0,i].pcolormesh(xedges,yedges,np.log10(Hmasked), cmap='jet')
    axs[0,i].tick_params(axis = 'both', labelsize=15) 
    #axs[0,i].set_xlabel('iMedian', fontsize=15)
    #axs[0,i].set_ylabel('iModel', fontsize=15)

    # bottom plot : iMed vs iPsf
    x = iMed
    y = iPsf
    maskLim = (x < xlim[1])&(x>xlim[0])&(y<ylim[1])&(y>ylim[0])
    maskExt = med_cat_merged_all['extendedness'] == ext
    maskCol = maskLim*maskExt
    H, xedges, yedges = np.histogram2d(x[maskCol], y[maskCol], bins=nbins)
    # H needs to be rotated and flipped
    H = np.rot90(H)
    H = np.flipud(H)
    # Mask zeros
    Hmasked = np.ma.masked_where(H==0,H) # Mask pixels with a value of zero
    # Plot 2D histogram using pcolor
    axs[1,i].plot(xlim, ylim, 'k-', alpha=0.75, zorder=0)
    axs[1,i].pcolormesh(xedges,yedges,np.log10(Hmasked), cmap='jet')
    axs[1,i].tick_params(axis = 'both', labelsize=15) 
    #axs[1,i].set_xlabel('iMedian', fontsize=15)
    #axs[1,i].set_ylabel('iPsf', fontsize=15)

# fig.text(x-coord, y-coord, **kwargs)
fig.text(0.25, 0.96, 'Extendedness=0', ha='center', va='center',fontsize=20)
fig.text(0.75, 0.96, 'Extendedness=1', ha='center', va='center',fontsize=20)

fig.text(0.5, 0.04, 'iMedian', ha='center', va='center',fontsize=20)
fig.text(0.03, 0.25, 'iPsf', ha='center', va='center', rotation='vertical',fontsize=20)
fig.text(0.03, 0.75, 'iModel', ha='center', va='center', rotation='vertical',fontsize=20)
fig.tight_layout()
fig.subplots_adjust(wspace=0.2, hspace=0.2, left=0.12, right=0.94, bottom=0.10, top=0.9)
#plt.tight_layout()
plt.savefig('Fig_iMedian_iPsf_iModel_all_RA.png')
plt.close(fig)


#   DO THE 3 PANEL PLOTS Fig23 Sesar_2010 : 
# g-i vs i  

# in the ipython, use 
# %cpaste   magick !!!! 
# http://stackoverflow.com/questions/10886946/how-does-ipythons-magic-paste-work

raMin = [320, 330, 340,350,0]
raMax = [330, 340, 350,360,10]


#####   g-i vs i   

# Select a range of gi,  i  for the histogram....  
ylim = [15,25.5]
xlim = [0,1.3]
#nbinsArr = [100,100,50]

nbins = 100

for i in range(2):
    ext=i # choose extendedness 
    fig, axs = plt.subplots(2,3, figsize=(18,12))
    ax = np.ravel(axs) 
    
    for j in range(len(raMin)):
        print('j=%d'%j)
        maskRa =  (raMin[j] < med_cat_merged_all['ra'])  & (med_cat_merged_all['ra'] < raMax[j])
        maskExtend = med_cat_merged_all['extendedness'] == ext
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
        
        ax[j].tick_params(axis = 'both', labelsize=15) 

    ax[j+1].axis('off')
    plt.tight_layout()
    plt.savefig('Fig_g-i_vs_i_extendedness_'+str(ext)+'.png')

#plt.tick_params(axis='both', which='major', labelsize=15)
#plt.gca().invert_yaxis()
#cbar = plt.colorbar()
#cbar.ax.set_ylabel(r'$\log_{10}{Counts}$')




#######   g-r   vs  r-i : choosing magnitude range (after Easter, NEW)

raMin = [320, 330, 340,350,0]
raMax = [330, 340, 350,360,10]

# Select a range of r-i,  g-r  for the histogram.... 
xlim = [-1,3]
ylim = [-1,3]
#nbinsArr = [100,100,50]

mag_lims = [5,22.4,22.6,22.8,23,23.2,23.4,23.6,24,35]

nbins = 100
count=0
for i in range(2):
    ext=i # choose extendedness 
    print('extendedness=%d'%ext)
    for k in range(len(mag_lims)-1):
        maskIband = (mag_lims[k] < med_cat_merged_all['median_mag_i_corr'])*(med_cat_merged_all['median_mag_i_corr'] < mag_lims[k+1])
        print('mag_min=%.2f, mag_max=%.2f'%(mag_lims[k], mag_lims[k+1]))
        # make a plot for each extendedness and mag range 
        count+=1
        fig, axs = plt.subplots(2,3, figsize=(18,12))
        ax = np.ravel(axs) 

        for j in range(len(raMin)):
            #print('j=%d'%j)
            
            maskRa =  (raMin[j] < med_cat_merged_all['ra'])  & (med_cat_merged_all['ra'] < raMax[j])
            maskExtend = med_cat_merged_all['extendedness'] == ext

            maskTot = maskRa * maskExtend * maskIband
            ri = med_cat_merged_all['median_mag_r_corr'][maskTot] - med_cat_merged_all['median_mag_i_corr'][maskTot]
            gr = med_cat_merged_all['median_mag_g_corr'][maskTot] - med_cat_merged_all['median_mag_r_corr'][maskTot]
            grri = pd.concat([gr,ri], axis=1)
            grri_dropna = grri.dropna()
            gr = grri_dropna[0]
            ri = grri_dropna[1]
            maskCol = (ri < xlim[1])&(ri>xlim[0])&(gr<ylim[1])&(gr>ylim[0])
            H, xedges, yedges = np.histogram2d(gr[maskCol],ri[maskCol],bins=nbins)
            #H, xedges, yedges = np.histogram2d(grri_dropna[1],grri_dropna[0],bins=nbins)
            # H needs to be rotated and flipped
            H = np.rot90(H)
            H = np.flipud(H)
            # Mask zeros
            Hmasked = np.ma.masked_where(H==0,H) # Mask pixels with a value of zero
            # Plot 2D histogram using pcolor
            ax[j].pcolormesh(xedges,yedges,np.log10(Hmasked), cmap='jet')
            #ax[j].set_ylabel('r-i', fontsize=15)
            #ax[j].set_xlabel('g-r', fontsize=15)
            ax[j].set_xlim(xlim[0],xlim[1])
            ax[j].set_ylim(ylim[0],ylim[1])
            ax[j].set_title(str(raMin[j])+'< RA <' + str(raMax[j]))
            #ax[j].invert_yaxis()

            ax[j].tick_params(axis = 'both', labelsize=15) 

        #fig.subplots_adjust(wspace=0, hspace=0, left=0.1, bottom=0.20, top=0.95)

        ax[j+1].axis('off')
        fig.text(0.5, 0.96, 'Extendedness='+str(ext)+' , '+str(mag_lims[k])+' < i < '+str(mag_lims[k+1]), ha='center', va='center',fontsize=20)
        fig.text(0.5, 0.04, 'g-r', ha='center', va='center',fontsize=20)
        fig.text(0.03, 0.5, 'r-i', ha='center', va='center', rotation='vertical',fontsize=20)

        fig.tight_layout()
        fig.subplots_adjust(wspace=0.26, hspace=0.26, left=0.12, right=0.94, bottom=0.10, top=0.9)

        plt.savefig('test'+str(count)+'_Fig_g-r_vs_r-i_ext_'+str(ext)+'_'+str(mag_lims[k])+'-i-'+str(mag_lims[k+1])+'.png')
        plt.close(fig)
#plt.show()



#######   g-r   vs  r-i :full range (old) 
#-->  CHECK FOR NANS.... 

# Select a range of r-i,  g-r  for the histogram.... 
xlim = [-1,3]
ylim = [-1,3]
#nbinsArr = [100,100,50]

nbins = 100

for i in range(2):
    ext=i # choose extendedness 
    fig, axs = plt.subplots(2,3, figsize=(18,12))
    ax = np.ravel(axs) 
    
    for j in range(len(raMin)):
        print('j=%d'%j)
        maskRa =  (raMin[j] < med_cat_merged_all['ra'])  & (med_cat_merged_all['ra'] < raMax[j])
        maskExtend = med_cat_merged_all['extendedness'] == ext
        maskTot = maskRa * maskExtend
        ri = med_cat_merged_all['median_mag_r_corr'][maskTot] - med_cat_merged_all['median_mag_i_corr'][maskTot]
        gr = med_cat_merged_all['median_mag_g_corr'][maskTot] - med_cat_merged_all['median_mag_r_corr'][maskTot]
        grri = pd.concat([gr,ri], axis=1)
        grri_dropna = grri.dropna()
        gr = grri_dropna[0]
        ri = grri_dropna[1]
        maskCol = (ri < xlim[1])&(ri>xlim[0])&(gr<ylim[1])&(gr>ylim[0])
        H, xedges, yedges = np.histogram2d(gr[maskCol],ri[maskCol],bins=nbins)
        #H, xedges, yedges = np.histogram2d(grri_dropna[1],grri_dropna[0],bins=nbins)
        # H needs to be rotated and flipped
        H = np.rot90(H)
        H = np.flipud(H)
        # Mask zeros
        Hmasked = np.ma.masked_where(H==0,H) # Mask pixels with a value of zeropltssssss
        # Plot 2D histogram using pcolor
        ax[j].pcolormesh(xedges,yedges,np.log10(Hmasked), cmap='jet')
        ax[j].set_ylabel('r-i', fontsize=15)
        ax[j].set_xlabel('g-r', fontsize=15)
        #ax[j].set_xlim(xlim[0],xlim[1])
        ax[j].set_title(str(raMin[j])+'< RA <' + str(raMax[j]))
        #ax[j].invert_yaxis()
        
        ax[j].tick_params(axis = 'both', labelsize=15) 

    #fig.subplots_adjust(wspace=0, hspace=0, left=0.1, bottom=0.20, top=0.95)
    #fig.text(0.5, 0.04, 'Right Ascension [degrees]', ha='center', va='center',fontsize=20)
    ax[j+1].axis('off')
    plt.tight_layout()
    plt.savefig('Fig_g-r_vs_r-i_extendedness_'+str(ext)+'_lims.png')
plt.close(fig)

# i vs RA : new way (after meeting Z before Easter,  
# plot all on one long stripe,  plot RA vs i ) 


# http://www.astropy.org/astropy-tutorials/plot-catalog.html
#import astropy.coordinates as coord
#import astropy.units as u
#ra = coord.Angle(table['ra']*u.degree)  # it should work but it takes so much time! 
# perhaps there  is a more workable solution... 
#ra = ra.wrap_at(180*u.degree)
#dec = coord.Angle(table['declination']*u.degree)
#ax.scatter(ra.radian, med_cat_merged_all['median_mag_i_corr']) 

raMin = [320, 330, 340,350,0]
raMax = [330, 340, 350,360,10]

nbins = 100
ylim = [12,26]

for i in range(2):
    print('ext = %d'% i)
    ext=i # choose extendedness 
    fig, axs = plt.subplots(1,5, figsize=(30,6))
    ax = np.ravel(axs) 
    
    for j in range(len(raMin)):
        print('j=%d'%j)
        maskRa =  (raMin[j] < med_cat_merged_all['ra'] )  & (med_cat_merged_all['ra']  < raMax[j])
        maskExtend = med_cat_merged_all['extendedness'] == ext
        maskTot = maskRa * maskExtend
        ra = med_cat_merged_all['ra'][maskTot]
        iBand = med_cat_merged_all['median_mag_i_corr'][maskTot]
        rai = pd.concat([ra,iBand], axis=1)
        rai_dropna = rai.dropna()
        ra = rai_dropna['ra'] 
        iBand = rai_dropna['median_mag_i_corr']
        #H, xedges, yedges = np.histogram2d(i[maskCol],ra[maskCol],bins=nbins)
        H, xedges, yedges = np.histogram2d(ra, iBand,bins=nbins)
        # H needs to be rotated and flipped
        H = np.rot90(H)
        H = np.flipud(H)
        # Mask zeros
        Hmasked = np.ma.masked_where(H==0,H) # Mask pixels with a value of zeropltssssss
        # Plot 2D histogram using pcolor
        ax[j].pcolormesh(xedges,yedges,np.log10(Hmasked), cmap='jet')
        #ax[j].set_xlabel('RA', fontsize=15)
        if j == 0 : 
            ax[j].set_ylabel('SDSS i', fontsize=20)
        ax[j].set_ylim(ylim[0],ylim[1])
        #ax[j].set_title(str(raMin[j])+'< RA <' + str(raMax[j]), fontsize=20)
        #ax[j].invert_yaxis()
        ax[j].tick_params(axis = 'both', labelsize=20) 
        if j != 0 : 
            ax[j].set_yticks([])
            xticks = ax[j].xaxis.get_major_ticks()
            xticks[0].label1.set_visible(False)
            ax[j].spines['left'].set_visible(False)
        #xticks[-1].label1.set_visible(False)
    #ax[j+1].axis('off')
    fig.subplots_adjust(wspace=0, hspace=0, left=0.1, bottom=0.20, top=0.95)
    fig.text(0.5, 0.04, 'Right Ascension [degrees]', ha='center', va='center',fontsize=20)
    #plt.tight_layout()
    plt.savefig('test_Fig_RA_vs_i_extendedness_'+str(ext)+'.png')
plt.close(fig)


#  i  vs RA   : old way (6 panels split by RA, and plot i vs RA) 
xlim = [-1,3]
ylim = [-1,3]

for i in range(2):
    print('ext = %d'% i)
    ext=i # choose extendedness 
    fig, axs = plt.subplots(2,3, figsize=(18,12))
    ax = np.ravel(axs) 
    
    for j in range(len(raMin)):
        print('j=%d'%j)
        maskRa =  (raMin[j] < med_cat_merged_all['ra'])  & (med_cat_merged_all['ra'] < raMax[j])
        maskExtend = med_cat_merged_all['extendedness'] == ext
        maskTot = maskRa * maskExtend
        ra = med_cat_merged_all['ra'][maskTot]
        iBand = med_cat_merged_all['median_mag_i_corr'][maskTot]
        rai = pd.concat([ra,iBand], axis=1)
        rai_dropna = rai.dropna()
        ra = rai_dropna['ra'] 
        iBand = rai_dropna['median_mag_i_corr']
	#H, xedges, yedges = np.histogram2d(i[maskCol],ra[maskCol],bins=nbins)
	H, xedges, yedges = np.histogram2d(iBand,ra,bins=nbins)
	# H needs to be rotated and flipped
	H = np.rot90(H)
	H = np.flipud(H)
	# Mask zeros
	Hmasked = np.ma.masked_where(H==0,H) # Mask pixels with a value of zeropltssssss
	# Plot 2D histogram using pcolor
	ax[j].pcolormesh(xedges,yedges,np.log10(Hmasked), cmap='jet')
	ax[j].set_xlabel('i', fontsize=15)
        ax[j].set_ylabel('RA', fontsize=15)
	#ax[j].set_xlim(xlim[0],xlim[1])
	ax[j].set_title(str(raMin[j])+'< RA <' + str(raMax[j]))
	#ax[j].invert_yaxis()
        ax[j].tick_params(axis = 'both', labelsize=15) 

    ax[j+1].axis('off')
    plt.tight_layout()
    plt.savefig('Fig_i_vs_RA_extendedness_'+str(ext)+'.png')
plt.close(fig)

