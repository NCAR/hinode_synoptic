from astropy.io import fits
from astropy.time import Time, TimeDelta
from astropy.coordinates import SkyCoord
from astropy import units as u
from scipy import ndimage
from sunpy.image import coalignment
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
import numpy as np
import drms, glob, os

# some machine- and user-dependent settings
drms_email     = 'cbethge@usra.edu'

data_directory = ['/Users/Christian/Desktop/HAO_Hinode_work/HOP/north/hinode/',   \
                  '/Users/Christian/Desktop/HAO_Hinode_work/HOP/equator/hinode/', \
                  '/Users/Christian/Desktop/HAO_Hinode_work/HOP/south/hinode/']

fits_directory = ['/Users/Christian/Desktop/HAO_Hinode_work/HOP/north/hmi_cycle25_cutouts/',   \
                  '/Users/Christian/Desktop/HAO_Hinode_work/HOP/equator/hmi_cycle25_cutouts/', \
                  '/Users/Christian/Desktop/HAO_Hinode_work/HOP/south/hmi_cycle25_cutouts/']

show_plots = 0

for jj in range(0,3):
    os.chdir(data_directory[jj])
    for fitsfile in glob.glob('*.fits'):
        # Read the necessary extensions from Hinode fits file
        # and compute the field.
        hdulist_hinode = fits.open(fitsfile)

        B     = hdulist_hinode[1].data
        theta = hdulist_hinode[2].data
        alpha = hdulist_hinode[12].data

        # For one file in the HOP south list, the filling factor was written as
        # a binary table instead of an image extension in the fits file, and
        # caused the code to crash. If that happens, just skip the file. 
        if isinstance(alpha, fits.fitsrec.FITS_rec):
            print('Skipping '+fitsfile)
            continue
        
        magnetic_flux_density = alpha*B*np.cos(np.radians(theta))

        # Create astropy time objects of the start and end times of the Hinode scan in UTC.
        start_time = Time(hdulist_hinode[0].header['TSTART'], format='isot', scale='utc')
        end_time   = Time(hdulist_hinode[0].header['TEND'], format='isot', scale='utc')

        # Determine time in the middle of the scan. Convert start_time and end_time to TAI,
        # because that is what HMI uses.
        middle_of_scan = start_time.tai+(end_time.tai-start_time.tai)/2

        # Convert to a time string that the DRMS understands.
        query_middle_of_scan = middle_of_scan.value.replace('-','.').replace('T','_')[:19]+'_TAI/12m'

        # Now query the DRMS for the HMI.ME_720s_fd10 datasets. Download the inclination and 
        # the field and get the fits header keywords that we need for the coordinates.
        c = drms.Client(email=drms_email, verbose=True)

        res = c.query('HMI.ME_720s_fd10['+query_middle_of_scan+']', \
                        key='T_REC, CRPIX1, CRPIX2, CRVAL1, CRVAL2, CDELT1, CDELT2, CROTA2', \
                        seg='inclination, field')

        # HMI timestamp in fits/isot format
        hmi_timestamp = Time(res[0]['T_REC'][0].replace('.','-').replace('_TAI','').replace('_','T'), scale='tai')

        # HMI timestamp in jd format
        hmi_timestamp_jd = Time(hmi_timestamp, format='jd', scale='tai')
        
        # Calculate time difference and give a warning if it is more than 12 minutes.
        time_diff = TimeDelta(hmi_timestamp-middle_of_scan, format='sec') 

        print('-----------------------------------------------------------')
        print('Hinode time (TAI, middle of scan): ', middle_of_scan)
        print('HMI time (TAI):                    ', hmi_timestamp)
        if time_diff.value > 720:
            print('*WARNING*: HMI/Hinode time more than 12 minutes apart.')
    
        export_request = 'HMI.ME_720s_fd10['+query_middle_of_scan+']{inclination, field}'

        r = c.export(export_request)
        
        hmi_inclination_hdu = fits.open(r.urls.url[0])
        hmi_field_hdu = fits.open(r.urls.url[1])

        hmi_inclination = hmi_inclination_hdu[1].data
        hmi_field = hmi_field_hdu[1].data

        # Now compute the flux density for HMI and rotate it with CROTA2. 
        magnetic_flux_density_hmi = hmi_field*np.cos(np.radians(hmi_inclination))
        magnetic_flux_density_hmi_rot = ndimage.rotate(np.nan_to_num(magnetic_flux_density_hmi), \
                                                        res[0]['CROTA2'][0], reshape=False)
            
        # Great, HMI is rotated by ~180 degrees, but only sometimes ¯\_(ツ)_/¯
        # Introduce a very simple switch for the calculation of the coordinates
        # here because I am too lazy to do this properly right now.
        if res[0]['CROTA2'][0] > 170.:
            tmp_hmi_xcoords = res[0]['CRVAL1'][0] + (res[0]['CRPIX1'][0] - \
                                                        np.arange(hmi_inclination.shape[0], dtype=float))*res[0]['CDELT1'][0]
            tmp_hmi_ycoords = res[0]['CRVAL2'][0] + (res[0]['CRPIX2'][0] - \
                                                        np.arange(hmi_inclination.shape[1], dtype=float))*res[0]['CDELT2'][0]
                                                        
            hmi_xcoords = tmp_hmi_xcoords*np.cos(np.radians(res[0]['CROTA2'][0])) - \
              tmp_hmi_ycoords*np.sin(np.radians(res[0]['CROTA2'][0])) 
            hmi_ycoords = tmp_hmi_xcoords*np.sin(np.radians(res[0]['CROTA2'][0])) + \
              tmp_hmi_ycoords*np.cos(np.radians(res[0]['CROTA2'][0])) 
        else:
            hmi_xcoords = res[0]['CRVAL1'][0] + (np.arange(hmi_inclination.shape[0], dtype=float) - \
                                                    res[0]['CRPIX1'][0])*res[0]['CDELT1'][0]
            hmi_ycoords = res[0]['CRVAL2'][0] + (np.arange(hmi_inclination.shape[1], dtype=float) - \
                                                    res[0]['CRPIX2'][0])*res[0]['CDELT2'][0]


        # Find out where the middle of the FOV should be by finding the
        # HMI timestamp (near the middle of the Hinode scan) in the
        # cycle 25 band, which is roughly defined by the decimal years
        # [2011.80, 2031.83] and Latitudes [55.0, 0.0].

        # First define astropy times (decimalyear, tai) for the cycle 25 band
        time1 = Time(2011.80, format='decimalyear', scale='tai')
        time2 = Time(2031.83, format='decimalyear', scale='tai')

        # Convert to JD
        time1.format = 'jd'
        time2.format = 'jd'
        
        # The latitudes of cycle 25
        lat1  = 55.
        lat2  = 0.
        
        # How many seconds in cycle 25?
        cyc25_in_sec = TimeDelta(time2-time1, format='sec') 
        cyc25_in_sec = np.around(cyc25_in_sec.value).astype(int)

        # Now build a latitude ramp with this many number of seconds on the x-axis
        latitudes = lat1 + np.linspace(0, cyc25_in_sec, cyc25_in_sec+1) * ((lat2-lat1)/cyc25_in_sec)

        # Find the second that corresponds to the middle of the Hinode scan
        which_sec_in_cyc25 = TimeDelta(hmi_timestamp_jd-time1, format='sec')
        which_sec_in_cyc25 = np.around(which_sec_in_cyc25.value).astype(int)
        
        # Take 300 arcseconds in x and 82 arcseconds in y to standardize the map size.
        # Calculate what that means in pixels.
        #hmi_cutsize_in_pix = [300./res[0]['CDELT1'][0], 82./res[0]['CDELT2'][0]]
        # Actually, just use 596 and 162 pixels, which corresponds to ~300.39 arcsec in x
        # and ~81.65 arcsec in y. CDELT does not change for HMI, and we probably want the
        # same map size for all cutouts. If we are rounding to the nearest pixel value
        # later, the map size could change.
        hmi_cutsize_in_pix = [596,162]
        
        # Heliographic-Stonyhurst coordinates for HMI FOV center on cycle 25 band
        hmi_cent_hgs = SkyCoord(0.*u.deg, latitudes[which_sec_in_cyc25]*u.deg, \
                                frame='heliographic_stonyhurst', observer='earth', obstime=hmi_timestamp) 
        
        # Convert to HPC
        this_hpc_coord = hmi_cent_hgs.transform_to("helioprojective")

        # The HMI array indices for the middle of the desired FOV
        this_hmi_index_center = [np.argmin(np.abs(np.array(hmi_xcoords)-np.around(this_hpc_coord.Tx.value))), \
                                 np.argmin(np.abs(np.array(hmi_ycoords)-np.around(this_hpc_coord.Ty.value)))]

        # The HMI array indices for the cutout
        this_hmi_cutout_index_x = [(this_hmi_index_center[0]-hmi_cutsize_in_pix[0]/2).astype(int), \
                                   (this_hmi_index_center[0]+hmi_cutsize_in_pix[0]/2).astype(int)]
        this_hmi_cutout_index_y = [(this_hmi_index_center[1]-hmi_cutsize_in_pix[1]/2).astype(int), \
                                   (this_hmi_index_center[1]+hmi_cutsize_in_pix[1]/2).astype(int)]                    

        # The actual cutout data
        this_hmi_data = magnetic_flux_density_hmi_rot[this_hmi_cutout_index_y[0]:this_hmi_cutout_index_y[1], \
                                                      this_hmi_cutout_index_x[0]:this_hmi_cutout_index_x[1]]

        if show_plots != 0:
            fig, ax = plt.subplots(1, 2, constrained_layout=True, figsize=(14,7))
            ax[0].imshow(magnetic_flux_density_hmi_rot, cmap="gray", vmin=-10., vmax=10., origin='lower')
            rect = Rectangle((this_hmi_cutout_index_x[0],this_hmi_cutout_index_y[0]), \
                              hmi_cutsize_in_pix[0], hmi_cutsize_in_pix[1],linewidth=1, edgecolor='r', facecolor='none')
            ax[0].add_patch(rect)                 
            ax[1].imshow(this_hmi_data, cmap="gray", vmin=-10., vmax=10., origin='lower')
            for ii in range(0,2):
                ax[ii].get_xaxis().set_visible(False)
                ax[ii].get_yaxis().set_visible(False)
            plt.show(block=False)
            plt.pause(1)
            plt.close(fig)
            
        # write HMI cutouts as fits files, add some info to the header
        hdr = fits.Header()
        hdr.set('hmi_xlow', hmi_xcoords[this_hmi_cutout_index_x[0]], ' lower HMI x-coordinate')
        hdr.set('hmi_xup',  hmi_xcoords[this_hmi_cutout_index_x[1]], ' upper HMI x-coordinate')
        hdr.set('hmi_ylow', hmi_ycoords[this_hmi_cutout_index_y[0]], ' lower HMI y-coordinate')
        hdr.set('hmi_yup',  hmi_ycoords[this_hmi_cutout_index_y[1]], ' upper HMI y-coordinate')
        hdr.set('hin_time', middle_of_scan.value[:-4], ' TAI timestamp middle of Hinode scan')
        hdr.set('hmi_time', hmi_timestamp.value[:-4], ' TAI timestamp of the HMI file')
    
        hdu = fits.PrimaryHDU(this_hmi_data, header=hdr)
        hdu.writeto(fits_directory[jj]+'HMI_'+hmi_timestamp.value[:-4].replace(':','').replace('-','').replace('T','_')+'.fits', \
                        overwrite=True)

            
