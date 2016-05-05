bin_start = 22
bin_end = 22.5
from matplotlib.pyplot import cm 
n = 5
color=iter(cm.rainbow(np.linspace(0,1,n)))

ylim = [-0.5,2.5]
mask_y = ((iPsf-iModel) < ylim[1]) & ((iPsf-iModel) > ylim[0])


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
