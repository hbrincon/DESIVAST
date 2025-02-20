import os
import numpy as np
import healpy as hp
import fitsio
from astropy.table import Table
import pickle

# Path to Iron VAC
# This path has been changed before, so make sure it is up to date when running the code for the first time in a while
iron_spec_path = "/global/cfs/cdirs/desi/public/dr1/vac/dr1/fastspecfit/iron/v2.1/catalogs/fastspec-iron.fits"

# Mask of 100% Survey Completeness
mask_file = "./mask/iron_mask.fits"
smoothed_mask_file = "./mask/iron_mask_smoothed.fits"


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
mask=fitsio.read(mask_file)
smoothed_mask=fitsio.read(smoothed_mask_file)
nside=hp.npix2nside(len(mask)) #healpix nside parameter

# Save infomration about the mask and masked galaxies
iron_info = {}
iron_info['file'] = iron_spec_path
iron_info['mask_file'] = mask_file
iron_info['mask'] = mask
iron_info['zmin'] = zmin 
iron_info['zmax'] = zmax 
iron_info['zmax_buffer'] = zmax_buffer 
iron_info['mag_lim'] = mag_lim
iron_info['mag_lim_buffer'] = mag_lim_buffer

#Open the Iron VAC
with fitsio.FITS(iron_spec_path) as iron, open(f'iron_info.pickle', 'wb') as iron_info_file:
    
    print("reading data")
    redshifts = iron[2]['TARGETID','Z','ZWARN','DELTACHI2','SPECTYPE','RA','DEC','BGS_TARGET', 'SURVEY', 'PROGRAM'
                       ][:]
    fastspec = iron[1][
                       'ABSMAG01_SDSS_R','ABSMAG01_IVAR_SDSS_R',
                      ][:]
    #E-correction
    fastspec['ABSMAG01_SDSS_R'] += 0.97*(redshifts['Z']-0.1)
    
    
    iron_info['file_count'] = len(fastspec) #target observations read in
    print(len(fastspec), "target observations read in")
    
    # Select BGS Bright galaxies
    select = np.where(
               (redshifts['SPECTYPE']=='GALAXY') 
               & (redshifts['BGS_TARGET'] & 2**1 != 0) #BGS bright target
               & (redshifts['SURVEY']=='main')
               & (redshifts['PROGRAM']=='bright')
             )
    
    redshifts=redshifts[select]
    fastspec=fastspec[select]
        
    #check for duplicate targets
    _, select = np.unique(redshifts['TARGETID'], return_index=True)
    if len(redshifts) != len(select):
        raise ValueError(f'Duplicate galaxies detected. {len(select)} out of {len(redshifts)} are unique')
            
    iron_info['bgs_bright_count'] = len(fastspec) #BGS Bright galaxies
    print(len(fastspec), "BGS Bright galaxies")

    # Impose redshift limits
    print("Imposing redshift limits")
    select = np.where((redshifts['Z']>zmin)  # > zmin and not >=zmin to avoid galaxies at origin
                    & (redshifts['Z']<=zmax_buffer) 
                 )
    select2 = np.where((redshifts['Z']>zmin)  # > zmin and not >=zmin to avoid galaxies at origin
                & (redshifts['Z']<=zmax) 
             )
    
    # Make two versions of the catalog, with and without the buffer values
    redshifts2=redshifts[select2]
    fastspec2=fastspec[select2]
    
    redshifts=redshifts[select]
    fastspec=fastspec[select]
    
    iron_info['zlim_count'] = len(fastspec2) #BGS Bright galaxies
    print(len(fastspec2), "galaxies in redshift limits")
    
    # Cut on (un-smoothed) mask
    #The pixel IDs for every (ra, dec) position
    pxid=hp.ang2pix(nside, redshifts['RA'], redshifts['DEC'], nest=True,lonlat=True)
    #Select the galaxies that fall within the mask
    select = np.isin(pxid, mask[mask["DONE"]==1]["HPXPIXEL"])
    redshifts=redshifts[select]
    fastspec=fastspec[select]
    
    #The pixel IDs for every (ra, dec) position
    pxid=hp.ang2pix(nside, redshifts2['RA'], redshifts2['DEC'], nest=True,lonlat=True)
    #Select the galaxies that fall within the mask
    select2 = np.isin(pxid, mask[mask["DONE"]==1]["HPXPIXEL"])
    redshifts2=redshifts2[select2]
    fastspec2=fastspec2[select2]
    
    
    iron_info['masked_count'] = len(fastspec2) #targets within angular mask
    print(len(fastspec2), "galaxies within angular mask")
    
    #Absolute magnitude cut
    select = np.where(fastspec['ABSMAG01_SDSS_R'] <= mag_lim_buffer) 
    select2 = np.where(fastspec2['ABSMAG01_SDSS_R'] <= mag_lim) 
    
    # Preserve two versions of the catalog, with and without the buffer values
    redshifts=redshifts[select]
    fastspec=fastspec[select]
    redshifts2=redshifts2[select2]
    fastspec2=fastspec2[select2]
    
    iron_info['vlim_count'] = len(fastspec2) #Volume limited Catalog
    print(len(fastspec2), "galaxies within magnitude cut")
    
    #Quality cuts 
    #made to match Ross 2024, The Dark Energy Spectroscopic Instrument: Construction of Large-scale Structure Catalogs
    select = np.where(
               (redshifts['ZWARN']==0) 
               & (redshifts['DELTACHI2']>40) 
    )   
    select2 = np.where(
               (redshifts2['ZWARN']==0) 
               & (redshifts2['DELTACHI2']>40) 
    ) 
    
    # Preserve two versions of the catalog, with and without the buffer values
    redshifts=redshifts[select]
    fastspec=fastspec[select]
    redshifts2=redshifts2[select2]
    fastspec2=fastspec2[select2]
    
    iron_info['final_count'] = len(fastspec2) #Quality cuts  
    print(len(fastspec2), "galaxies in final catalog")
    
    
    #split into NGC and SGC
    select = (redshifts['RA'] < 304) * (redshifts['RA'] > 83)
    select2 = (redshifts2['RA'] < 304) * (redshifts2['RA'] > 83)

    iron_info['final_count_ngc'] = len(fastspec2[select2])
    print(len(fastspec2[select2]),'galaxies in NGC')
    iron_info['final_count_sgc'] = len(fastspec2[~select2])
    print(len(fastspec2[~select2]),'galaxies in SGC')
    
    print(len(fastspec),'galaxies in final catalog with magnitude and redshift buffer')
    iron_info['final_count_buffer'] = len(fastspec)
    
    #count how many galaxies are in smoothed mask
    #The pixel IDs for every (ra, dec) position
    pxid=hp.ang2pix(nside, redshifts2['RA'], redshifts2['DEC'], nest=True,lonlat=True)
    #Select the galaxies that fall within the mask
    select3 = np.isin(pxid, smoothed_mask[smoothed_mask["DONE"]==1]["HPXPIXEL"])
    iron_info['smoothed_final_count'] = len(redshifts2[select3])
    iron_info['smoothed_final_count_ngc'] = len(redshifts2[select3*select2])
    iron_info['smoothed_final_count_sgc'] =  len(redshifts2[select3*~select2])
    pxid=hp.ang2pix(nside, redshifts['RA'], redshifts['DEC'], nest=True,lonlat=True)
    select3 = np.isin(pxid, smoothed_mask[smoothed_mask["DONE"]==1]["HPXPIXEL"])
    
    #save catalog metadata
    pickle.dump(iron_info, iron_info_file)
    
    #save catalog
    out=Table([redshifts['TARGETID'],
               redshifts['RA'],
               redshifts['DEC'],
               redshifts['Z'],
               fastspec['ABSMAG01_SDSS_R'],
               fastspec['ABSMAG01_IVAR_SDSS_R'],
              ],
               
               names=['targetID','ra','dec','redshift',
                        'rabsmag','rabsmag_inv'])
    #unsmoothed mask catalog
    out[select].write('iron_ngc.fits', format='fits', overwrite=True) 
    out[~select].write('iron_sgc.fits', format='fits', overwrite=True) 
    
    #smoothed mask catalog
    out[select3*select].write('iron_smoothed_ngc.fits', format='fits', overwrite=True) 
    out[select3*~select].write('iron_smoothed_sgc.fits', format='fits', overwrite=True) 

