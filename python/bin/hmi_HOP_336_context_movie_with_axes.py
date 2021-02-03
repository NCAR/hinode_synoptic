from astropy.io import fits
from astropy.time import Time, TimeDelta
from astropy.coordinates import SkyCoord
from astropy import units as u
from scipy import ndimage
from sunpy.image import coalignment
import datetime as dt
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
import numpy as np
import drms, glob, os

def get_image_shift(image1, image2):
    # slightly modified code snippet from bitbucket to determine the shift
    # between two images using FFT:
    
    # discrete fast fourier transformation and complex conjugation of image2 
    image1FFT = np.fft.fft2(image1)
    image2FFT = np.conjugate( np.fft.fft2(image2) )

    # inverse fourier transformation of product -> equal to cross correlation
    imageCCor = np.real( np.fft.ifft2( (image1FFT*image2FFT) ) )

    # shift the zero-frequency component to the center of the spectrum
    imageCCorShift = np.fft.fftshift(imageCCor)

    # determine the distance of the maximum from the center
    row, col = image1.shape

    yshift, xshift = np.unravel_index( np.argmax(imageCCorShift), (row,col) )

    xshift -= int(col/2)
    yshift -= int(row/2)

    return xshift, yshift

# some machine- and user-dependent settings
drms_email     = 'cbethge@usra.edu'

data_directory = ['/Users/Christian/Desktop/HAO_Hinode_work/HOP/north/hinode/',   \
                  '/Users/Christian/Desktop/HAO_Hinode_work/HOP/equator/hinode/', \
                  '/Users/Christian/Desktop/HAO_Hinode_work/HOP/south/hinode/']

plot_directory = '/Users/Christian/Desktop/HAO_Hinode_work/HOP/context_movie_plots/'

show_plots       = 0
cadence_in_hours = 3

HOP_336_start_date = dt.date(2017, 3, 7)
HOP_336_end_date   = dt.date(2020, 1, 21)

# Calculate the duration of the HOP 336 campaign
HOP_336_duration = HOP_336_end_date - HOP_336_start_date

# Query the DRMS for the entire HOP 336 campaign at the given cadence
c = drms.Client(email=drms_email, verbose=True)
res = c.query('HMI.ME_720s_fd10['+str(HOP_336_start_date.year)+'.'+str(HOP_336_start_date.month)+'.'+     \
              str(HOP_336_start_date.day)+'/'+str(HOP_336_duration.days)+'d@'+str(cadence_in_hours)+'h]', \
                        key='T_REC, CRPIX1, CRPIX2, CRVAL1, CRVAL2, CDELT1, CDELT2, CROTA2', \
                        seg='inclination, field')

# Save the HMI times from the query in a list
hmi_times = res[0]['T_REC'].to_list()

# Make a list of time strings in ISO format
hmi_times_iso = []
for timestamp in hmi_times:
    timestamp_iso = timestamp[0:4]+'-'+timestamp[5:7]+'-'+timestamp[8:10]+'T'+timestamp[11:19]
    hmi_times_iso.append(timestamp_iso)

# Convert that list to an astropy time object 
hmi_times_iso = Time(hmi_times_iso, format='isot', scale='tai')
    
# Now we get all the information we need from the Hinode HOP 336 files
current_dir = os.getcwd() 

HOP_336_info = []
print('Getting times and coordinates from Hinode HOP 336 files...')
for jj in range(0,3):
    os.chdir(data_directory[jj])
    for fitsfile in glob.glob('*.fits'):
        # Read the necessary extensions from Hinode fits file
        # and compute the field.
        hdulist_hinode = fits.open(fitsfile)

        # Create astropy time objects of the start and end times of the Hinode scan in UTC.
        start_time = Time(hdulist_hinode[0].header['TSTART'], format='isot', scale='utc')
        end_time   = Time(hdulist_hinode[0].header['TEND'], format='isot', scale='utc')

        # Determine time in the middle of the scan. Convert start_time and end_time to TAI,
        # because that is what HMI uses.
        middle_of_scan = start_time.tai+(end_time.tai-start_time.tai)/2

        # For this Hinode file, find the closest HMI time and index at the given cadence
        tdeltas = TimeDelta(hmi_times_iso - middle_of_scan, format='sec')
        closest_hmi_index = np.abs(tdeltas.value).argmin()
        closest_hmi_time  = hmi_times_iso.value[closest_hmi_index]

        #Here I am sorting out bad coordinates with quantiles, the example below assumes that a
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
            else:    
                # add to list: HMI index, HMI timestamp, Path to Hinode file, Hinode timestamp, Hinode xcoords, Hinode ycoords 
                HOP_336_info.append([closest_hmi_index, closest_hmi_time, data_directory[jj]+fitsfile, \
                                         middle_of_scan.value, hinode_xcoords, hinode_ycoords])
        else:    
            # add to list: HMI index, HMI timestamp, Path to Hinode file, Hinode timestamp, Hinode xcoords, Hinode ycoords 
            HOP_336_info.append([closest_hmi_index, closest_hmi_time, data_directory[jj]+fitsfile, \
                                     middle_of_scan.value, hinode_xcoords, hinode_ycoords])

os.chdir(current_dir)

# Now go through all of the HMI images from the query and plot them
for ii in res[0].index:
    print('Processing index '+str(ii)+' of '+str(res[0].index.stop))
    
    inclination_url = 'http://jsoc.stanford.edu'+res[1]['inclination'][ii]
    field_url       = 'http://jsoc.stanford.edu'+res[1]['field'][ii]
        
    hmi_inclination_hdu = fits.open(inclination_url)
    hmi_field_hdu = fits.open(field_url)

    hmi_inclination = hmi_inclination_hdu[1].data
    hmi_field = hmi_field_hdu[1].data

    # Now compute the flux density for HMI and rotate it with CROTA2. 
    magnetic_flux_density_hmi = hmi_field*np.cos(np.radians(hmi_inclination))
    magnetic_flux_density_hmi_rot = ndimage.rotate(np.nan_to_num(magnetic_flux_density_hmi), \
                                                    res[0]['CROTA2'][ii], reshape=False)

    # Great, HMI is rotated by ~180 degrees, but only sometimes ¯\_(ツ)_/¯
    # Introduce a very simple switch for the calculation of the coordinates
    # here because I am too lazy to do this properly right now.
    if res[0]['CROTA2'][ii] > 170.:
        tmp_hmi_xcoords = res[0]['CRVAL1'][ii] + (res[0]['CRPIX1'][ii] - \
                                                      np.arange(hmi_inclination.shape[0], dtype=float))*res[0]['CDELT1'][ii]
        tmp_hmi_ycoords = res[0]['CRVAL2'][ii] + (res[0]['CRPIX2'][ii] - \
                                                      np.arange(hmi_inclination.shape[1], dtype=float))*res[0]['CDELT2'][ii]
                                                        
        hmi_xcoords = tmp_hmi_xcoords*np.cos(np.radians(res[0]['CROTA2'][ii])) - \
          tmp_hmi_ycoords*np.sin(np.radians(res[0]['CROTA2'][ii])) 
        hmi_ycoords = tmp_hmi_xcoords*np.sin(np.radians(res[0]['CROTA2'][ii])) + \
          tmp_hmi_ycoords*np.cos(np.radians(res[0]['CROTA2'][ii])) 
    else:
        hmi_xcoords = res[0]['CRVAL1'][ii] + (np.arange(hmi_inclination.shape[0], dtype=float) - \
                                                  res[0]['CRPIX1'][ii])*res[0]['CDELT1'][ii]
        hmi_ycoords = res[0]['CRVAL2'][ii] + (np.arange(hmi_inclination.shape[1], dtype=float) - \
                                                  res[0]['CRPIX2'][ii])*res[0]['CDELT2'][ii]
    
    # Now plot the HMI image  
    fig, ax = plt.subplots(1, 1, constrained_layout=True, figsize=(10,10))
    ax.imshow(magnetic_flux_density_hmi_rot, cmap="gray", vmin=-10., vmax=10., origin='lower', \
                  extent=[hmi_xcoords.min(),hmi_xcoords.max(),hmi_ycoords.min(),hmi_ycoords.max()], aspect='auto')
    # ax.get_xaxis().set_visible(False)
    # ax.get_yaxis().set_visible(False)
    ax.title.set_text('HMI '+hmi_times_iso.value[ii])
    
    # If one of the HOP 336 maps has this HMI image as the nearest one,
    # plot a rectangle with the Hinode coordinates. First we need to
    # establish the coordinates though and get the data for the cross-
    # correlation...
    for jj in range(0,len(HOP_336_info)):
        if ii == HOP_336_info[jj][0]:
            # Read the necessary extensions from Hinode fits file
            # and compute the field.
            hdulist_hinode = fits.open(HOP_336_info[jj][2])

            B     = hdulist_hinode[1].data
            theta = hdulist_hinode[2].data
            alpha = hdulist_hinode[12].data
        
            magnetic_flux_density = alpha*B*np.cos(np.radians(theta))

            hinode_xcoords = HOP_336_info[jj][4]
            hinode_ycoords = HOP_336_info[jj][5]
            
            # Find out which indices in HMI the Hinode coordinates correspond to.
            hmi_index_x = [np.argmin(np.abs(np.array(hmi_xcoords)-np.floor(hinode_xcoords[0]))), \
                            np.argmin(np.abs(np.array(hmi_xcoords)-np.ceil(hinode_xcoords[1])))]
                        
            hmi_index_y = [np.argmin(np.abs(np.array(hmi_ycoords)-np.floor(hinode_ycoords[0]))), \
                            np.argmin(np.abs(np.array(hmi_ycoords)-np.ceil(hinode_ycoords[1])))]        
        
            # Degrade the Hinode data with a Gaussian and rescale to the same pixel size as HMI.
            # The 1.3 sigma for the Gaussian is not actually chosen based on any physical parameters,
            # just such that the Hinode map looks ~ish like HMI.
            tmp_hinode_data = ndimage.gaussian_filter(magnetic_flux_density,1.3)
            tmp_hinode_data = ndimage.zoom(tmp_hinode_data, \
                                            [(hmi_index_y[1]-hmi_index_y[0])/tmp_hinode_data.shape[0], \
                                                (hmi_index_x[1]-hmi_index_x[0])/tmp_hinode_data.shape[1]])
            
            # Now add 150 pixels padding on either side to allow for the Hinode offset wrt HMI. 
            # Technically this could hit the HMI array boundaries and there is no failsafe 
            # implemented, but HMI maps are so large that for the HOP maps we will probably
            # not hit the sides. At least it does not for the 153 HOP north maps...
            pix_pad = 150 
            hmi_data_template = magnetic_flux_density_hmi_rot[hmi_index_y[0]-pix_pad:hmi_index_y[1]+pix_pad, \
                                                                hmi_index_x[0]-pix_pad:hmi_index_x[1]+pix_pad]

            # Use the sunpy coalignment routine to get an initial idea where the degraded Hinode data
            # should be on the padded HMI map.  
            thisxyshift = coalignment.calculate_shift(hmi_data_template, tmp_hinode_data, repair_nonfinite=False)

            # Use that information to create a hopefully more precise HMI cutout.
            hmi_index_x_new = np.around(hmi_index_x+(thisxyshift[0].value-pix_pad)).astype(int)
            hmi_index_y_new = np.around(hmi_index_y+(thisxyshift[1].value-pix_pad)).astype(int)
            hmi_data = magnetic_flux_density_hmi_rot[hmi_index_y_new[0]:hmi_index_y_new[1],hmi_index_x_new[0]:hmi_index_x_new[1]]

            # Now use the FFT image shift routine for better precision and
            # create the final HMI cutout. 
            xshift, yshift = get_image_shift(hmi_data, tmp_hinode_data)

            hmi_index_x_new += xshift
            hmi_index_y_new += yshift

            rect = Rectangle((hmi_xcoords[hmi_index_x_new[0]],hmi_ycoords[hmi_index_y_new[0]]), \
                              hmi_xcoords[hmi_index_x_new[1]]-hmi_xcoords[hmi_index_x_new[0]],  \
                              hmi_ycoords[hmi_index_y_new[1]]-hmi_ycoords[hmi_index_y_new[0]],  \
                              linewidth=1, edgecolor='r', facecolor='none')
            
            ax.add_patch(rect)
            # Plot the Hinode timestamp next to the rectangle 
            #ax.text(hmi_index_x_new[0], hmi_index_y_new[1]+15, HOP_336_info[jj][3], fontsize=10, color='red')
            # Wasn't really legible, so put it at the top next to the title instead
            plt.gcf().text(0.68, 0.983, 'Hinode: '+HOP_336_info[jj][3], fontsize=12)


    fig.savefig(plot_directory+'{:06d}'.format(ii)+'.jpg', dpi=107, format='jpg', quality=95)        
    if show_plots != 0:
        plt.show(block=False)
        plt.pause(1)

    plt.close(fig)
