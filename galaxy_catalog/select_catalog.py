import os
import numpy as np
import healpy as hp
import fitsio
from astropy.table import Table, vstack, hstack
import numpy.lib.recfunctions as rfn
import pickle

# Mask of 100% Survey Completeness
mask_file = "./mask/loa_mask.fits"
smoothed_mask_file = "./mask/loa_mask_smoothed.fits"


# Redshift limits
zmin = 0.
zmax = 0.24
zmax_buffer = 0.25 #an extra buffer of galaxies at higher redshift for testing purposes

# Magnitude limit
mag_lim= -20.
mag_lim_buffer = -19.5 #an extra buffer of dim galaxies for testing purposes


# Read in mask
# There are two mask versions. The first ("mask") covers all regions of BGS that have been completed (for up to 4 passes).
# The second ("smoothed_mask") removes small non-contiguous regions from the mask for the purpose of void-finding on a contigious mask
# From now on, we will use "mask" to refer to "smoothed mask"
mask=fitsio.read(smoothed_mask_file)
#smoothed_mask=fitsio.read(smoothed_mask_file)
nside=hp.npix2nside(len(mask)) #healpix nside parameter

def get_spec_file_name (healpix_index):

    formated_index = str(healpix_index).zfill(2)
    
    file_name = f'/global/cfs/cdirs/desi/vac/dr2/fastspecfit/loa/v1.0/catalogs/fastspec-loa-main-bright-nside1-hp{formated_index}.fits'

    return file_name

# Save information about the mask and masked galaxies
catalog_info = {}
catalog_info['file'] = get_spec_file_name('XX')
catalog_info['mask_file'] = mask_file
catalog_info['mask'] = mask
catalog_info['zmin'] = zmin 
catalog_info['zmax'] = zmax 
catalog_info['zmax_buffer'] = zmax_buffer 
catalog_info['mag_lim'] = mag_lim
catalog_info['mag_lim_buffer'] = mag_lim_buffer

catalog_info['file_count'] = 0
catalog_info['bgs_bright_count'] = 0
catalog_info['zlim_count'] = 0
catalog_info['masked_count'] = 0
catalog_info['vlim_count'] = 0
catalog_info['final_count'] = 0
catalog_info['final_count_ngc'] = 0
catalog_info['final_count_sgc'] = 0
catalog_info['final_count_buffer'] = 0

ngc_catalog = Table(names = ['targetID','ra','dec','redshift',
                          'uabsmag','gabsmag', 'rabsmag', 'zabsmag',
                          'logmstar', 'sfr', 'halphaew'])

sgc_catalog = Table(names = ['targetID','ra','dec','redshift',
                          'uabsmag','gabsmag', 'rabsmag', 'zabsmag',
                          'logmstar', 'sfr', 'halphaew'])

#Open the FastSpecFit VAC
for idx in range (12):
    
    with fitsio.FITS(get_spec_file_name(idx)) as full_catalog:

        #TODO: ADD MORE COLORS AND GALAXY PROPERTIES FOR DR2 ENV ANALYSIS
        
        print("reading data", idx)
        metadata = full_catalog[1]['TARGETID','Z','ZWARN','DELTACHI2','SPECTYPE','RA','DEC','BGS_TARGET', 'SURVEY', 'PROGRAM'
                           ][:]
        specphot = full_catalog[2][
                           'ABSMAG01_SDSS_U', 'ABSMAG01_SDSS_G',
                           'ABSMAG01_SDSS_R', 'ABSMAG01_SDSS_Z',
                           'LOGMSTAR', 'SFR',
                          ][:]

        fastspec = full_catalog[3]['HALPHA_EW',][:]

        calculated_data = np.recarray( (len(metadata), 1), dtype=[('IN_SAMPLE', '<i8'), ('ABSMAG_R_ECORR', '<f8')])

        rabsmag_ecorr = specphot['ABSMAG01_SDSS_R'] + 0.97*(metadata['Z']-0.1)
        calculated_data['ABSMAG_R_ECORR'] = rabsmag_ecorr.reshape((len(metadata),1))

        calculated_data['IN_SAMPLE'] = 1

        catalog = rfn.merge_arrays([metadata, specphot, fastspec, calculated_data], flatten=True, usemask=False)

        del metadata, specphot, fastspec, calculated_data, rabsmag_ecorr

        catalog_info['file_count'] += len(catalog) #target observations read in
        print(len(catalog), "target observations read in")
        
        # Select BGS Bright galaxies
        select = np.where(
                   (catalog['SPECTYPE']=='GALAXY') 
                   & (catalog['BGS_TARGET'] & 2**1 != 0) #BGS bright target
                   & (catalog['SURVEY']=='main')
                   & (catalog['PROGRAM']=='bright')
                 )
        
        catalog=catalog[select]
            
        #check for duplicate targets
        _, select = np.unique(catalog['TARGETID'], return_index=True)
        if len(catalog) != len(select):
            raise ValueError(f'Duplicate galaxies detected. {len(select)} out of {len(catalog)} are unique')
                
        catalog_info['bgs_bright_count'] += len(catalog) #BGS Bright galaxies
        print(len(catalog), "BGS Bright galaxies")
    
        # Impose redshift limits
        print("Imposing redshift limits")
        select = np.where((catalog['Z']>zmin)  # > zmin and not >=zmin to avoid galaxies at origin
                        & (catalog['Z']<=zmax_buffer) 
                     )

        catalog=catalog[select]

        select = catalog['Z'] > zmax   # > zmin and not >=zmin to avoid galaxies at origin
        
        catalog['IN_SAMPLE'][select] = 0 # mark buffer galaxies

        num_selected = np.sum(catalog['IN_SAMPLE'])
        
        catalog_info['zlim_count'] += num_selected #BGS Bright galaxies
        print(num_selected, "galaxies in redshift limits")
        
        # Cut on mask
        #The pixel IDs for every (ra, dec) position
        pxid=hp.ang2pix(nside, catalog['RA'], catalog['DEC'], nest=True,lonlat=True)
        #Select the galaxies that fall within the mask
        select = np.isin(pxid, mask[mask["DONE"]==1]["HPXPIXEL"])
        catalog=catalog[select]

        num_selected = np.sum(catalog['IN_SAMPLE'])
        
        catalog_info['masked_count'] += num_selected #targets within angular mask
        print(num_selected, "galaxies within smoothed angular mask")
        
        #Absolute magnitude cut
        select = catalog['ABSMAG_R_ECORR'] <= mag_lim_buffer

        catalog=catalog[select]
        
        select = catalog['ABSMAG_R_ECORR'] > mag_lim

        catalog['IN_SAMPLE'][select] = 0 # mark buffer galaxies

        num_selected = np.sum(catalog['IN_SAMPLE'])
        
        catalog_info['vlim_count'] += num_selected #Volume limited Catalog
        print(num_selected, "galaxies within magnitude cut")
        
        #Quality cuts 
        #made to match Ross 2024, The Dark Energy Spectroscopic Instrument: Construction of Large-scale Structure Catalogs
        select = np.where(
                   (catalog['ZWARN']==0) 
                   & (catalog['DELTACHI2']>40) 
        )   
        
        catalog=catalog[select]

        num_selected = np.sum(catalog['IN_SAMPLE'])
        
        catalog_info['final_count'] += num_selected #Quality cuts  
        print(num_selected, "galaxies in final catalog")
        
        
        #split into NGC and SGC
        select = (catalog['RA'] < 304) * (catalog['RA'] > 83)

        num_selected = np.sum(select*catalog['IN_SAMPLE'])
        catalog_info['final_count_ngc'] += num_selected
        print(num_selected,'galaxies in NGC')
        num_selected = np.sum(~select*catalog['IN_SAMPLE'])
        catalog_info['final_count_sgc'] += num_selected
        print(num_selected,'galaxies in SGC')
        
        print(len(catalog),'galaxies in final catalog with magnitude and redshift buffer')
        catalog_info['final_count_buffer'] += len(catalog)

        """
        #count how many galaxies are in smoothed mask
        #The pixel IDs for every (ra, dec) position
        pxid=hp.ang2pix(nside, catalog2['RA'], catalog2['DEC'], nest=True,lonlat=True)
        #Select the galaxies that fall within the mask
        select3 = np.isin(pxid, smoothed_mask[smoothed_mask["DONE"]==1]["HPXPIXEL"])
        catalog_info['smoothed_final_count'] += len(redshifts2[select3])
        catalog_info['smoothed_final_count_ngc'] += len(redshifts2[select3*select2])
        catalog_info['smoothed_final_count_sgc'] +=  len(redshifts2[select3*~select2])
        pxid=hp.ang2pix(nside, redshifts['RA'], redshifts['DEC'], nest=True,lonlat=True)
        select3 = np.isin(pxid, smoothed_mask[smoothed_mask["DONE"]==1]["HPXPIXEL"])
        """

        #save catalog
        out=Table([catalog['TARGETID'],
                   catalog['RA'],
                   catalog['DEC'],
                   catalog['Z'],
                   catalog['ABSMAG01_SDSS_U'],
                   catalog['ABSMAG01_SDSS_G'],
                   catalog['ABSMAG01_SDSS_R'],
                   catalog['ABSMAG01_SDSS_Z'],
                   catalog['LOGMSTAR'], 
                   catalog['SFR'], 
                   catalog["HALPHA_EW"]
                  ],
                   
                   names=['targetID','ra','dec','redshift',
                          'uabsmag','gabsmag', 'rabsmag', 'zabsmag',
                          'logmstar', 'sfr', 'halphaew'])

        ngc_catalog = vstack([ngc_catalog, out[select]])
        sgc_catalog = vstack([sgc_catalog, out[~select]])

with open(f'loa_info.pickle', 'wb') as catalog_info_file:
    #save catalog metadata
    pickle.dump(catalog_info, catalog_info_file)


"""#unsmoothed mask catalog
out[select].write('loa_ngc.fits', format='fits', overwrite=True) 
out[~select].write('loa_sgc.fits', format='fits', overwrite=True) 

#smoothed mask catalog
out[select3*select].write('loa_smoothed_ngc.fits', format='fits', overwrite=True) 
out[select3*~select].write('loa_smoothed_sgc.fits', format='fits', overwrite=True)"""

ngc_catalog.write('loa_ngc.fits', format='fits', overwrite=True) 
sgc_catalog.write('loa_sgc.fits', format='fits', overwrite=True) 

