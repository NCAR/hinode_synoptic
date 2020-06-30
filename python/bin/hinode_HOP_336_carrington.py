from scipy.io import readsav
from astropy.time import Time
from astropy.io import fits
from astropy.coordinates import SkyCoord
from astropy import units as u
from sunpy.coordinates import sun
import matplotlib.pyplot as plt
import numpy as np
import glob, os

# Read in the Carrington data from an IDL sav file
filename = 'HMImagbfly.sav'
sav_data = readsav(filename)

# Full Carrington map. For whatever reason, the time axis is for future Carrington
# rotations at the end, so cutting those off (last index 139).
carr_full = sav_data['bfly'][:,:140]

# Stretch by a factor of 10 on the time axis for a finer Carrington rotation axis
carr_full_stretched = np.repeat(carr_full, 10, axis=1)

# We are only interested in the CRs 2187 to 2231, which corresponds to the indices 95 to 139
carr_of_interest = sav_data['bfly'][:,95:140]

# Stretch by a factor of 10 on the time axis for a finer Carrington rotation axis
carr_of_interest_stretched = np.repeat(carr_of_interest, 10, axis=1)

# The data directories for the HOP336 data
data_directory = ['/Users/Christian/Desktop/HAO_Hinode_work/HOP/north/hinode/',   \
                  '/Users/Christian/Desktop/HAO_Hinode_work/HOP/equator/hinode/', \
                  '/Users/Christian/Desktop/HAO_Hinode_work/HOP/south/hinode/']

# Path and names of the plots
plotfile1 = '/Users/Christian/Desktop/Hinode_carrington_plots/Carrington_plot_full'
plotfile2 = '/Users/Christian/Desktop/Hinode_carrington_plots/HOP336_carrington_plot_full'
plotfile3 = '/Users/Christian/Desktop/Hinode_carrington_plots/Carrington_plot_zoom'
plotfile4 = '/Users/Christian/Desktop/Hinode_carrington_plots/HOP336_carrington_plot_zoom'

# Plot with secondary axis? 0=no, 1=yes
plot_sec_axis = 1

current_dir = os.getcwd() 

# Go through all of the Hinode data and save Carrington rotation and latitudes
# (bottom and top of the map) in a list 
HOP336_carr = []
for jj in range(0,3):
    os.chdir(data_directory[jj])
    for fitsfile in glob.glob('*.fits'):
        # Open the fits file
        hdulist_hinode = fits.open(fitsfile)

        # Create astropy time objects of the start and end times of the Hinode scan in UTC.
        start_time = Time(hdulist_hinode[0].header['TSTART'], format='isot', scale='utc')
        end_time   = Time(hdulist_hinode[0].header['TEND'], format='isot', scale='utc')

        # Determine time in the middle of the scan (still in UTC). 
        middle_of_scan = start_time+(end_time-start_time)/2

        # Determine the Carrington rotation of this scan, round to the first significant
        # figure (which corresponds to the stretching of 10 we made above)
        this_carr_number = np.round(sun.carrington_rotation_number(t=middle_of_scan), decimals=1)

        # Here I am sorting out bad coordinates with quantiles, the example below assumes that a
        # maximum of 2% of the coordinates are bad (1% at the lower end, 1% at the upper end). 
        # That fixes most of the maps, but not all of them. For now, I am just adding/subtracting 
        # an arcsec on either end to counteract the (probable) overcompensation, i.e. 
        # throwing away good coordinate values as well. If the HMI cutout needs to be
        # *really* precise, this needs to be done more carefully. Right now, the HMI map
        # might be a bit too small or too large depending on how many coordinate outliers
        # there are.
        hinode_xcoords = [np.quantile(hdulist_hinode[38].data,0.01)-1., np.quantile(hdulist_hinode[38].data,0.99)+1.]
        hinode_ycoords = [np.quantile(hdulist_hinode[39].data,0.01)-1., np.quantile(hdulist_hinode[39].data,0.99)+1.]

        # Take care of some problematic files in the list manually
        if fitsfile in ['20170402_210600.fits', '20170403_003405.fits', '20170828_102201.fits', \
                        '20171009_181600.fits', '20171127_175734.fits', '20180401_215905.fits', \
                        '20180611_161204.fits', '20190708_221322.fits']:
            tmp_xcoords = hdulist_hinode[38].data
            tmp_ycoords = hdulist_hinode[39].data
            if fitsfile == '20170402_210600.fits':
                hinode_xcoords = [tmp_xcoords[tmp_xcoords > (-250)].min()-5., tmp_xcoords[tmp_xcoords < 200].max()]
                hinode_ycoords = [tmp_ycoords[tmp_ycoords > 500].min()-5., tmp_ycoords[tmp_ycoords < 800].max()]
            if fitsfile == '20170403_003405.fits':
                hinode_xcoords = [tmp_xcoords[tmp_xcoords > (-250)].min()-5., tmp_xcoords[tmp_xcoords < 400].max()]
                hinode_ycoords = [tmp_ycoords[tmp_ycoords > 500].min()-5., tmp_ycoords[tmp_ycoords < 800].max()]
            if fitsfile == '20170828_102201.fits':
                hinode_xcoords = [tmp_xcoords[tmp_xcoords > (-200)].min()-5., tmp_xcoords[tmp_xcoords < 200].max()]
                hinode_ycoords = [tmp_ycoords[tmp_ycoords > 500].min()-5., tmp_ycoords[tmp_ycoords < 800].max()]
            if fitsfile == '20171009_181600.fits':
                hinode_xcoords = [tmp_xcoords[tmp_xcoords > (-400)].min()-10., tmp_xcoords[tmp_xcoords < 100].max()]
                hinode_ycoords = [tmp_ycoords[tmp_ycoords > (-800)].min()-5., tmp_ycoords[tmp_ycoords < (-600)].max()]
            # The following files I am simply skipping. The coordinates seem completely off, and searching
            # the entire HMI field-of-view for a specific patch of Quiet Sun is next to impossible.
            if fitsfile in ['20171127_175734.fits', '20180401_215905.fits', \
                            '20180611_161204.fits', '20190708_221322.fits']:
                print('Skipping '+fitsfile)
                continue
            
        # Calculate Heliographic-Stonyhurst coordinates for bottom and top 
        # of the Hinode scan. x is at the center of the scan.
        hin_bottom_hpc = SkyCoord(((hinode_xcoords[0]+hinode_xcoords[1])/2)*u.arcsec, \
                                    hinode_ycoords[0]*u.arcsec, \
                                    frame='helioprojective', observer='earth', obstime=middle_of_scan) 
        hin_bottom_hgs = hin_bottom_hpc.transform_to('heliographic_stonyhurst')

        hin_top_hpc = SkyCoord(((hinode_xcoords[0]+hinode_xcoords[1])/2)*u.arcsec, \
                                hinode_ycoords[1]*u.arcsec, \
                                frame='helioprojective', observer='earth', obstime=middle_of_scan) 
        hin_top_hgs = hin_top_hpc.transform_to('heliographic_stonyhurst')
       
        # Add to list: Carrington rotation number, sine of latitude of bottom and top
        HOP336_carr.append([this_carr_number, np.sin(np.radians(hin_bottom_hgs.lat.value)), \
                            np.sin(np.radians(hin_top_hgs.lat.value))])

os.chdir(current_dir)

# And now plot everything
minmax = [-10.,10.]

# Full range
fig, ax = plt.subplots(constrained_layout=True, figsize=(8,5))
ax.imshow(carr_full_stretched, cmap="gray", vmin=minmax[0], vmax=minmax[1], origin='lower', \
              extent=[2092.,2231.,-1,1], aspect='auto')
ax.set_xlabel('Carrington rotation')
ax.set_ylabel('sin(Latitude)')

if plot_sec_axis != 0:
    from matplotlib.ticker import FuncFormatter
    ax2 = ax.twiny()
    formatter_x = FuncFormatter(lambda x, pos: Time(sun.carrington_rotation_time(x), scale='utc', format='isot').value[0:10])
    ax2.xaxis.set_major_formatter(formatter_x)
    ax2.set_xlim(ax.get_xlim())
    ax2.set_xlabel('Date')
    ax2.xaxis.set_tick_params(labelsize=8)

    ax3 = ax.twinx()
    formatter_y = FuncFormatter(lambda x, pos: '{:0.2f}'.format((np.arcsin(x)/np.pi)*180.))
    ax3.yaxis.set_major_formatter(formatter_y)
    ax3.set_ylim(ax.get_ylim())
    ax3.set_ylabel('Latitude [degrees]')

    plotfile1 = plotfile1+'_sec_axis'
    plotfile2 = plotfile2+'_sec_axis'
    
plt.savefig(plotfile1+'.jpg', dpi=200, format='jpg', quality=95)

# Overplot the lines on the Carrington map. obs[0] is the Carrington rotation number,
# obs[1] the sine of the latitude at the bottom, obs[2] the sine of the latitude at the top.
for obs in HOP336_carr:
    ax.plot([obs[0], obs[0]], [obs[1], obs[2]], '-', linewidth=1, color='firebrick')

#ax.set_title('Hinode HOP336 observations')    
plt.savefig(plotfile2+'.jpg', dpi=200, format='jpg', quality=95)

# Now the zoom in
fig, ax = plt.subplots(constrained_layout=True, figsize=(8,5))
ax.imshow(carr_of_interest_stretched, cmap="gray", vmin=minmax[0], vmax=minmax[1], origin='lower', \
              extent=[2187.,2231.,-1,1], aspect='auto')
ax.set_xlabel('Carrington rotation')
ax.set_ylabel('sin(Latitude)')

if plot_sec_axis != 0:
    ax2 = ax.twiny()
    ax2.xaxis.set_major_formatter(formatter_x)
    ax2.set_xlim(ax.get_xlim())
    ax2.set_xlabel('Date')
    ax2.xaxis.set_tick_params(labelsize=8)

    ax3 = ax.twinx()
    ax3.yaxis.set_major_formatter(formatter_y)
    ax3.set_ylim(ax.get_ylim())
    ax3.set_ylabel('Latitude [degrees]')

    plotfile3 = plotfile3+'_sec_axis'
    plotfile4 = plotfile4+'_sec_axis'
    
plt.savefig(plotfile3+'.jpg', dpi=200, format='jpg', quality=95)

# Overplot the lines on the Carrington map. obs[0] is the Carrington rotation number,
# obs[1] the sine of the latitude at the bottom, obs[2] the sine of the latitude at the top.
for obs in HOP336_carr:
    ax.plot([obs[0], obs[0]], [obs[1], obs[2]], '-', linewidth=1, color='firebrick')

#ax.set_title('Hinode HOP336 observations')    
plt.savefig(plotfile4+'.jpg', dpi=200, format='jpg', quality=95)
