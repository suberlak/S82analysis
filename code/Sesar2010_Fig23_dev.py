# -*- coding: iso-8859-1 -*-
# make all necessary imports....
import pandas as pd
import numpy as np 
import matplotlib.pyplot as plt 


dir_info_files = '/astro/store/scratch/tmp/suberlak/S13Agg/repo_fls/'
dir_save_merged= '/astro/users/suberlak/S13Agg_analysis/data_products/'
dir_deep_source = '/astro/store/pogo4/s13_stripe82/deep_source/'

# the ugriz photometry is already extinction-corrected  
# we read in only those objects whose parents are brighter than 17 mag 
data_source = 'NCSA'
name_keep = 'test_'+data_source+'_all_patches_ugriz_var.csv'
name_detection = 'DeepSource'+data_source+'_i_lt235.csv.gz'

ncsa_photometry_data = pd.read_csv(dir_save_merged+name_keep)
ncsa_detection_data = pd.read_csv(dir_deep_source+name_detection, usecols=['deepSourceId', 'ra', 'decl', 'extendedness'])

# merge in the extendedness, ra, dec, extendedness   (there is psfFlux and ModelFlux,  but not modelMag.  To translate 
# that to modelMag would need to convert the flux to magnitudes - as long as I'm, using the same threshold for extendedness 
# as Yusra, it is not needed)

merged_ncsa = pd.merge(ncsa_photometry_data, ncsa_detection_data, how='left', left_on='objectId', right_on='deepSourceId')


data_source = 'IN2P3'
name_keep = 'test_'+data_source+'_all_patches_ugriz_var.csv'
name_detection = 'DeepSource'+data_source+'_i_lt235.csv.gz'
intpt_photometry_data = pd.read_csv(dir_save_merged+name_keep)
intpt_detection_data = pd.read_csv(dir_deep_source+name_detection)

# merge in the extendedness, ra, dec, extendedness  
merged_intpt = pd.merge(intpt_photometry_data, intpt_detection_data,how='left', left_on='objectId', right_on='deepSourceId')

# make sure that both have the same columns 
assert merged_ncsa.columns == merged_intpt.columns 

# combine data from NCSA and IN2P3  
merged_all = merged_ncsa.append(merged_intpt)


# make the plot ... 
raMin = [30, 10, 350,330]
raMax = [40, 20, 360,340]


#####   g-i vs r (like Fig.23)   

# Select a range of gi,  r  for the histogram....  
gi_lim = [0,1.3]
r_lim = [15,23.5] 
#nbinsArr = [100,100,50]

nbins = 100

for i in range(2):
    i = 0 
    ext=i # choose extendedness 
    fig, axs = plt.subplots(2,3, figsize=(18,12))
    ax = np.ravel(axs) 

    for j in range(len(raMin)):
        print('j=%d'%j)
        maskRa =  (raMin[j] < merged_all['ra'])  & (merged_all['ra'] < raMax[j])
        maskExtend = merged_all['extendedness'] == ext
        maskTot = maskRa * maskExtend
        gi = merged_all['gpsfMean_corr'][maskTot] - merged_all['ipsfMean_corr'][maskTot]
        r =  merged_all['rpsfMean_corr'][maskTot]
        maskCol = (gi < gi_lim[1])&(gi>gi_lim[0])&(r<r_lim[1])&(r > r_lim[0])
        H, xedges, yedges = np.histogram2d(gi[maskCol],r[maskCol],bins=nbins)
        # H needs to be rotated and flipped
        H = np.rot90(H)
        H = np.flipud(H)
        # Mask zeros
        Hmasked = np.ma.masked_where(H==0,H) # Mask pixels with a value of zeropltssssss
        # Plot 2D histogram using pcolor
        ax[j+1].pcolormesh(xedges,yedges,np.log10(Hmasked), cmap='jet')
        ax[j+1].set_xlabel('g-i', fontsize=15)
        ax[j+1].set_ylabel('r', fontsize=15)
        ax[j+1].set_xlim(gi_lim[0],gi_lim[1])
        ax[j+1].set_ylim(r_lim[0], r_lim[1])
        ax[j+1].set_title(str(raMin[j])+'< RA <' + str(raMax[j]))
        ax[j+1].invert_yaxis()
        
        ax[j+1].tick_params(axis = 'both', labelsize=15) 
    ax[0].axis('off')
    ax[5].axis('off')
    plt.tight_layout()
    plt.savefig('Sesar2010_Fig_23_g-i_vs_r_ext_'+str(ext)+'_NCSA_IN2P3_comb.png')


#####  g-r vs r-i   (the new plot) 

raMin = [320, 330, 340,350,0, 10,20,30,40]
raMax = [330, 340, 350,360,10,20,30,40,55]

# Select a range of r-i,  g-r  for the histogram.... 
xlim = [-1,3]
ylim = [-1,3]

mag_lims = [5,19.0,19.2,19.4,19.6,19.8, 20.0, 20.2, 20.4, 20.6, 20.8, 21.0, 21.2, 21.4, 21.6, 21.8, 22.0, 22.4,22.6,22.8,23,23.2,23.4,23.6,24]
#mag_lims = [22.0, 22.4, 22.6]
dir_save = '/astro/users/suberlak/S13Agg_analysis/data_products/Sesar_2010/'

nbins = 100
count=0
for i in range(2):
    ext=i # choose extendedness 
    print('extendedness=%d'%ext)

    for k in range(len(mag_lims)-1):
        mask_band = (mag_lims[k] < merged_all['ipsfMean_corr']) & (merged_all['ipsfMean_corr'] < mag_lims[k+1])
        print('i_mag_min=%.2f, i_mag_max=%.2f'%(mag_lims[k], mag_lims[k+1]))
        # make a plot for each extendedness and mag range 
        count+=1
        fig, axs = plt.subplots(3,3, figsize=(12,12))
        ax = np.ravel(axs) 

        for j in range(len(raMin)):
            #print('j=%d'%j)
            mask_ra =  (raMin[j] < merged_all['ra'])  & (merged_all['ra'] < raMax[j])
            mask_extendedness = merged_all['extendedness'] == ext

            mask_total = mask_ra & mask_extendedness & mask_band
            ri = merged_all['rpsfMean_corr'][mask_total] - merged_all['ipsfMean_corr'][mask_total]
            gr = merged_all['gpsfMean_corr'][mask_total] - merged_all['rpsfMean_corr'][mask_total]
            grri = pd.concat([gr,ri], axis=1)
            grri_dropna = grri.dropna()
            gr = grri_dropna[0]
            ri = grri_dropna[1]
            maskCol = (gr < xlim[1])&(gr>xlim[0])&(ri<ylim[1])&(ri>ylim[0])
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

        #ax[j+1].axis('off')
        fig.text(0.5, 0.96, 'Extendedness='+str(ext)+' , '+str(mag_lims[k])+' < iPsfMean < '+str(mag_lims[k+1]), ha='center', va='center',fontsize=20)
        fig.text(0.5, 0.04, 'g-r', ha='center', va='center',fontsize=20)
        fig.text(0.03, 0.5, 'r-i', ha='center', va='center', rotation='vertical',fontsize=20)

        fig.tight_layout()
        fig.subplots_adjust(wspace=0.26, hspace=0.26, left=0.12, right=0.94, bottom=0.10, top=0.9)

        plt.savefig(dir_save+ 'test'+str(count).zfill(2)+'_Fig_g-r_vs_r-i_ext_'+str(ext)+'_'+str(mag_lims[k])+'-i-'+str(mag_lims[k+1])+'.png')
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

raMin = [320, 330, 340,350,0, 10,20,30,40]
raMax = [330, 340, 350,360,10,20,30,40,55]

nbins = 100
ylim = [12,26]

for i in range(2):
    print('ext = %d'% i)
    ext=i # choose extendedness 
    fig, axs = plt.subplots(1,9, figsize=(30,6))
    ax = np.ravel(axs) 
    
    for j in range(len(raMin)):
        print('j=%d'%j)
        maskRa =  (raMin[j] < merged_all['ra'] )  & (merged_all['ra']  < raMax[j])
        maskExtend = merged_all['extendedness'] == ext
        maskTot = maskRa * maskExtend
        ra = merged_all['ra'][maskTot]
        iBand = merged_all['ipsfMean_corr'][maskTot]
        rai = pd.concat([ra,iBand], axis=1)
        rai_dropna = rai.dropna()
        ra = rai_dropna['ra'] 
        iBand = rai_dropna['ipsfMean_corr']
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


