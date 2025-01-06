import os
import numpy as np
import healpy as hp
import fitsio
from astropy.table import Table
import pickle
import kcorrect.kcorrect
#import kcorrect.response
#kcorrect.response.all_responses()

# Path to Kibo database
# This path has been changed before, so make sure it is up to date when running the code for the first time in a while
#iron_spec_path = "/global/cfs/cdirs/desi/public/dr1/vac/dr1/fastspecfit/iron/v2.1/catalogs/fastspec-iron.fits"
iron_spec_path = '/global/cfs/cdirs/desi/spectro/redux/kibo/zcatalog/v1/ztile-main-bright-cumulative.fits'
# Mask of 100% Survey Completeness
mask_file = "./mask/kibo_mask.fits"
smoothed_mask_file = "./mask/kibo_mask_smoothed.fits"


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
with fitsio.FITS(iron_spec_path) as iron, open(f'kibo_info.pickle', 'wb') as iron_info_file:
    
    print("reading data")
    redshifts = iron[1]['TARGETID','Z','ZWARN','DELTACHI2','SPECTYPE','TARGET_RA','TARGET_DEC','BGS_TARGET', 
                        'FLUX_R', 'FLUX_IVAR_R', 'FLUX_G', 'FLUX_IVAR_G', 'FLUX_Z', 'FLUX_IVAR_Z', 'FLUX_W1', 'FLUX_IVAR_W1', 'FLUX_W2', 'FLUX_IVAR_W2',
                        'FIBERFLUX_R'
                       ][:]
    
    print('Data loaded')
        
    iron_info['file_count'] = len(redshifts) #target observations read in
    print(len(redshifts), "target observations read in")
    
    # Select BGS Bright galaxies
    select = np.where(
               (redshifts['SPECTYPE']=='GALAXY') 
               & (redshifts['BGS_TARGET'] & 2**1 != 0) #BGS bright
             )
    
    redshifts=redshifts[select]
    iron_info['bgs_bright_count'] = len(redshifts) #BGS Bright galaxies
    print(len(redshifts), "BGS Bright galaxies")
    
    #setup column for absmag_r
    colnames = list(redshifts.dtype.names)
    colnames=list(map(lambda x: x.replace('FIBERFLUX_R','ABSMAG01_SDSS_R'),colnames))
    colnames=list(map(lambda x: x.replace('TARGET_RA','RA'),colnames))
    colnames=list(map(lambda x: x.replace('TARGET_DEC','DEC'),colnames))
    redshifts.dtype.names = colnames
    
    #needed for abs mag calculation to work
    redshifts=redshifts[(redshifts['Z']>zmin)*(redshifts['Z']<zmax_buffer)*(redshifts['FLUX_R']>0)]

    #k correct flux
    responses = ['sdss_g0', 'sdss_r0', 'sdss_z0', 'wise_w1', 'wise_w2']
    kc = kcorrect.kcorrect.Kcorrect(responses=responses)
    
    
    maggies = np.array([redshifts['FLUX_G'], redshifts['FLUX_R'], redshifts['FLUX_Z'], redshifts['FLUX_W1'], redshifts['FLUX_W2']]).T * 1e-9
    
    ivar = np.array([redshifts['FLUX_IVAR_G'], redshifts['FLUX_IVAR_R'], redshifts['FLUX_IVAR_Z'], redshifts['FLUX_IVAR_W1'], redshifts['FLUX_IVAR_W2']]).T * 1e-9

    coeffs = kc.fit_coeffs(redshift=redshifts['Z'], maggies=maggies, ivar=ivar)

    #calculate absolute magnitudes
    absmag = kc.absmag(redshift=redshifts['Z'], maggies=maggies, ivar=ivar, coeffs=coeffs, band_shift=0.1)
    
    print('Absolute magnitude calculated')
    
    
    #K-corrected/E-corrected r-band absolute magnitude
    redshifts['ABSMAG01_SDSS_R'] = absmag[:,1] + 0.97*(redshifts['Z']-0.1)

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
    
    redshifts=redshifts[select]
    
    iron_info['zlim_count'] = len(redshifts2) #BGS Bright galaxies
    print(len(redshifts2), "galaxies in redshift limits")
    
    # Cut on (un-smoothed) mask
    #The pixel IDs for every (ra, dec) position
    pxid=hp.ang2pix(nside, redshifts['RA'], redshifts['DEC'], nest=True,lonlat=True)
    #Select the galaxies that fall within the mask
    select = np.isin(pxid, mask[mask["DONE"]==1]["HPXPIXEL"])
    redshifts=redshifts[select]
    
    #The pixel IDs for every (ra, dec) position
    pxid=hp.ang2pix(nside, redshifts2['RA'], redshifts2['DEC'], nest=True,lonlat=True)
    #Select the galaxies that fall within the mask
    select2 = np.isin(pxid, mask[mask["DONE"]==1]["HPXPIXEL"])
    redshifts2=redshifts2[select2]
    
    
    iron_info['masked_count'] = len(redshifts2) #targets within angular mask
    print(len(redshifts2), "galaxies within angular mask")
    
    #Absolute magnitude cut
    select = np.where(redshifts['ABSMAG01_SDSS_R'] <= mag_lim_buffer) 
    select2 = np.where(redshifts2['ABSMAG01_SDSS_R'] <= mag_lim) 
    
    # Preserve two versions of the catalog, with and without the buffer values
    redshifts=redshifts[select]
    redshifts2=redshifts2[select2]
    
    iron_info['vlim_count'] = len(redshifts2) #Volume limited Catalog
    print(len(redshifts2), "galaxies within magnitude cut")
    
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
    redshifts2=redshifts2[select2]
    
    iron_info['final_count'] = len(redshifts2) #Quality cuts  
    print(len(redshifts2), "galaxies in final catalog")
    
    
    #split into NGC and SGC
    select = (redshifts['RA'] < 304) * (redshifts['RA'] > 83)
    select2 = (redshifts2['RA'] < 304) * (redshifts2['RA'] > 83)

    iron_info['final_count_ngc'] = len(redshifts2[select2])
    print(len(redshifts2[select2]),'galaxies in NGC')
    iron_info['final_count_sgc'] = len(redshifts2[~select2])
    print(len(redshifts2[~select2]),'galaxies in SGC')
    
    print(len(redshifts),'galaxies in final catalog with magnitude and redshift buffer')
    iron_info['final_count_buffer'] = len(redshifts)
    
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
               redshifts['ABSMAG01_SDSS_R'],
              ],
               
               names=['targetID','ra','dec','redshift',
                        'rabsmag'])
    
    ##unsmoothed mask catalog
    ##out[select].write('iron_ngc.fits', format='fits', overwrite=True) 
    ##out[~select].write('iron_sgc.fits', format='fits', overwrite=True) 
    
    #smoothed mask catalog
    out[select3*select].write('kibo_smoothed_ngc.fits', format='fits', overwrite=True) 
    out[select3*~select].write('kibo_smoothed_sgc.fits', format='fits', overwrite=True) 


#for debugging

"""

import os
import numpy as np
import healpy as hp
import fitsio
from astropy.table import Table
import pickle
import kcorrect.kcorrect
#import kcorrect.response
#kcorrect.response.all_responses()

# Path to Kibo database
# This path has been changed before, so make sure it is up to date when running the code for the first time in a while
#iron_spec_path = "/global/cfs/cdirs/desi/public/dr1/vac/dr1/fastspecfit/iron/v2.1/catalogs/fastspec-iron.fits"
iron_spec_path = '/global/cfs/cdirs/desi/spectro/redux/kibo/zcatalog/v1/ztile-main-bright-cumulative.fits'

#Open the Iron VAC
with fitsio.FITS(iron_spec_path) as iron:
    
    print("reading data")
    redshifts = iron[1]['TARGETID','Z','ZWARN','DELTACHI2','SPECTYPE','TARGET_RA','TARGET_DEC','BGS_TARGET', 
                        'FLUX_R', 'FLUX_IVAR_R', 'FLUX_G', 'FLUX_IVAR_G', 'FLUX_Z', 'FLUX_IVAR_Z', 'FLUX_W1', 'FLUX_IVAR_W1', 'FLUX_W2', 'FLUX_IVAR_W2',
                        'FIBERFLUX_R'
                       ][:]
    
redshifts=redshifts[(redshifts['Z']>0)*(redshifts['Z']<.25)*(redshifts['FLUX_R']>0)]

responses = ['sdss_g0', 'sdss_r0', 'sdss_z0', 'wise_w1', 'wise_w2']
kc = kcorrect.kcorrect.Kcorrect(responses=responses)


red = redshifts[80963: 80964]


maggies = np.array([red['FLUX_G'], red['FLUX_R'], red['FLUX_Z'], red['FLUX_W1'], red['FLUX_W2']]).T * 1e-9

ivar = np.array([red['FLUX_IVAR_G'], red['FLUX_IVAR_R'], red['FLUX_IVAR_Z'], red['FLUX_IVAR_W1'], red['FLUX_IVAR_W2']]).T * 1e-9

print(maggies)
print(ivar)
coeffs = kc.fit_coeffs(redshift=red['Z'], maggies=maggies, ivar=ivar)
print(coeffs)
#calculate absolute magnitudes
absmag = kc.absmag(redshift=red['Z'], maggies=maggies, ivar=ivar, coeffs=coeffs, band_shift=0.1)
print(absmag)

"""