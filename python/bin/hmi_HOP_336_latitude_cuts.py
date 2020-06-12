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
                  '/Users/Christian/Desktop/HAO_Hinode_work/HOP/south/hinode/']

fits_directory = ['/Users/Christian/Desktop/HAO_Hinode_work/HOP/north/hmi_latitude_cutouts/',   \
                  '/Users/Christian/Desktop/HAO_Hinode_work/HOP/south/hmi_latitude_cutouts/']

show_plots      = 0
along_longitude = 1

for jj in range(0,2):
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

        # Show the full disk HMI to see if it looks alright. 
        #plt.figure()
        #plt.imshow(magnetic_flux_density_hmi_rot, cmap="gray", vmin=-10., vmax=10., origin='lower') 
        #plt.show()

        # Initially I just used the minimum and maximum 
        # coordinate values to get the cutout: 
        #hinode_xcoords = [hdulist_hinode[38].data.min(), hdulist_hinode[38].data.max()]
        #hinode_ycoords = [hdulist_hinode[39].data.min(), hdulist_hinode[39].data.max()]

        # However, turns out that often there seem to be outliers in the coordinates. I tried
        # sorting those out with quantiles, the example below assumes that a maximum of 2%
        # of the coordinates are bad (1% at the lower end, 1% at the upper end). That fixes
        # most of the maps, but not all of them. For now, I am just adding/subtracting 
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

        # Mark the datasets as suspicious where the full-res Hinode and the degraded Hinode
        # have different aspect ratios. These are the ones where the quantiles did not sort
        # out all of the bad coordinates. Should not happen anymore though after the 
        # problematic files are taken care of manually above. 
        suspect = 0
        if magnetic_flux_density.shape[0]/magnetic_flux_density.shape[1] >= 1.:
            if tmp_hinode_data.shape[0]/tmp_hinode_data.shape[1] < 1.:
                print(fitsfile+' suspicious.')
                suspect = 1
        else:
            if tmp_hinode_data.shape[0]/tmp_hinode_data.shape[1] >= 1.:
                print(fitsfile+' suspicious.')
                suspect = 1
            
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

        # Now use the FFT image shift routine for better precision
        xshift, yshift = get_image_shift(hmi_data, tmp_hinode_data)

        hmi_index_x_new += xshift
        hmi_index_y_new += yshift

        # The arcsec center of the HMI FOV
        hmi_center = [(hmi_xcoords[hmi_index_x_new[0]]+hmi_xcoords[hmi_index_x_new[1]])/2, \
                      (hmi_ycoords[hmi_index_y_new[0]]+hmi_ycoords[hmi_index_y_new[1]])/2]

        # The index center in x of the HMI FOV, in case we don't go along the Longitude
        if along_longitude == 0:
            hmi_center_index_x = np.argmin(np.abs(np.array(hmi_xcoords)-np.around(hmi_center[0]))).astype(int)
            
        # Take 300 arcseconds in x and 82 arcseconds in y to standardize the map size.
        # Calculate what that means in pixels.
        #hmi_cutsize_in_pix = [300./res[0]['CDELT1'][0], 82./res[0]['CDELT2'][0]]
        # Actually, just use 596 and 162 pixels, which corresponds to ~300.39 arcsec in x
        # and ~81.65 arcsec in y. CDELT does not change for HMI, and we probably want the
        # same map size for all cutouts. If we are rounding to the nearest pixel value
        # later, the map size could change.
        hmi_cutsize_in_pix = [596,162]
        
        # Calculate Heliographic-Stonyhurst coordinates for HMI FOV center
        hmi_cent_hpc = SkyCoord(hmi_center[0]*u.arcsec, \
                                hmi_center[1]*u.arcsec, \
                                frame='helioprojective', observer='earth', obstime=hmi_timestamp) 
        hmi_cent_hgs = hmi_cent_hpc.transform_to('heliographic_stonyhurst')

        # create primary HDU for the fits file
        hdr = fits.Header()
        hdr.set('hin_time', middle_of_scan.value[:-4], ' TAI timestamp middle of Hinode scan')
        hdr.set('hmi_time', hmi_timestamp.value[:-4], ' TAI timestamp of the HMI file')
        hdr.set('suspect',  suspect, ' probably wrong HMI cutout (0=no, 1=yes)')
        hdr.set('alonglon', along_longitude, ' following the longitude (0=no, 1=yes)')
        
        hdu = fits.PrimaryHDU(header=hdr)
        this_filename = fits_directory[jj]+'HMI_latitude_cuts_'+ \
                        hmi_timestamp.value[:-4].replace(':','').replace('-','').replace('T','_')+'.fits'
        hdu.writeto(this_filename, overwrite=True)
        
        # Now we want 21 cutouts at that Longitude, starting at Lat 50 deg North going down
        # every 2.5 degrees (or -50 deg South going up every 2.5 degrees)
        for cc in range(0,21):
            if jj == 0:
                this_y_hgs_coord = (50.-2.5*cc)*u.deg
            else:
                this_y_hgs_coord = (-50.+2.5*cc)*u.deg

            # The middle of the FOV for the desired cutout in HGS coordinates
            this_hgs_coord = SkyCoord(hmi_cent_hgs.lon, \
                                      this_y_hgs_coord, \
                                      frame='heliographic_stonyhurst', observer='earth', obstime=hmi_timestamp)

            # Convert to HPC
            this_hpc_coord = this_hgs_coord.transform_to("helioprojective")

            # Distinguish if we are going along the Longitude or along a straight line in y
            if along_longitude != 0:
                # The HMI array indices for the middle of the desired FOV
                this_hmi_index_center = [np.argmin(np.abs(np.array(hmi_xcoords)-np.around(this_hpc_coord.Tx.value))), \
                                         np.argmin(np.abs(np.array(hmi_ycoords)-np.around(this_hpc_coord.Ty.value)))]
            else:
                # The HMI array indices for the middle of the desired FOV
                this_hmi_index_center = [hmi_center_index_x, \
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
            
            # Append the data to the fits file as an extension.
            # The extension name is the current latitude.
            this_hdr = fits.Header()
            this_hdr.set('EXTNAME', str(this_y_hgs_coord.value))
            this_hdr.set('hmi_xlow', hmi_xcoords[this_hmi_cutout_index_x[0]], ' lower HMI x-coordinate')
            this_hdr.set('hmi_xup',  hmi_xcoords[this_hmi_cutout_index_x[1]], ' upper HMI x-coordinate')
            this_hdr.set('hmi_ylow', hmi_ycoords[this_hmi_cutout_index_y[0]], ' lower HMI y-coordinate')
            this_hdr.set('hmi_yup',  hmi_ycoords[this_hmi_cutout_index_y[1]], ' upper HMI y-coordinate')
            fits.append(this_filename, this_hmi_data, this_hdr)

            
