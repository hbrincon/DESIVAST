import numpy as np
from astropy.table import Table, vstack
from astropy.io import fits
import json
import time

import sys
sys.path.insert(1, '/global/homes/h/hrincon/python_tools')
import VoidVolume as vol
import Env as env
from VoidCatalog import VoidFinderCatalog, VoidFinderCatalogStacked, V2Catalog, V2CatalogStacked
import VoidCatalog as vc
from multiprocessing import Process, Manager


masks_path = '/global/homes/h/hrincon/DESIVAST/galaxy_catalog/mask/alt_masks/masks.fits'
hr4_path = '/global/homes/h/hrincon/DESIVAST/mocks/HorizonRun4/'

manager = Manager()
gal_frac = manager.dict()


def mock_path(idx):
    return hr4_path+f'DR7m_{idx}.fits'
def vf_path(idx):
    return hr4_path+f'HR4_{idx}_VoidFinder_Output.fits' 
def vv_path(idx):
    return hr4_path+f'HR4_{idx}_V2_VIDE_Output.fits'
def vr_path(idx):
    return hr4_path+f'HR4_{idx}_V2_REVOLVER_Output.fits'


def get_overlap_vf(idx):

    vfc = VoidFinderCatalog(vf_path(idx))
    vfc.add_galaxies(mock_path(idx), redshift_name='z')
    
    vfc.info['MAGLIM'] = -20.0
    vfc.galaxies['rabsmag'] = -21.0
    vfc.calculate_vflag()
    masks = fits.open(masks_path)

    tmp = vfc.galaxy_membership(masks['SDSS'])

    frac = 100 * tmp[0] / tmp[1]
    
    gal_frac[f"VF ({idx})"] = frac
    
def get_overlap_vv(idx):
    
    v2v = V2Catalog(vv_path(idx))
    v2v.add_galaxies(mock_path(idx), redshift_name='z') 
    v2v.info['MAGLIM'] = -20.0
    v2v.galaxies['rabsmag'] = -21.0
    #v2v.add_mask(vfc_ref)
    
    masks = fits.open(masks_path)
    
    tmp = v2v.galaxy_membership(masks['SDSS'])

    frac = 100 * tmp[0] / tmp[1]
    
    gal_frac[f"VV ({idx})"] = frac

    
def get_overlap_vr(idx):

    v2r = V2Catalog(vr_path(idx))
    v2r.add_galaxies(mock_path(idx), redshift_name='z') 
    v2r.info['MAGLIM'] = -20.0
    v2r.galaxies['rabsmag'] = -21.0
    #v2r.add_mask(vfc_ref)
    
    masks = fits.open(masks_path)
    
    tmp = v2r.galaxy_membership(masks['SDSS'])

    frac = 100 * tmp[0] / tmp[1]
    
    gal_frac[f"VR ({idx})"] = frac
    
    
def get_overlap_vf_vv(idx):

    vfc = VoidFinderCatalog(vf_path(idx))
    vfc.add_galaxies(mock_path(idx), redshift_name='z') 
    vfc.info['MAGLIM'] = -20.0
    vfc.galaxies['rabsmag'] = -21.0

    v2v = V2Catalog(vv_path(idx))
    v2v.add_galaxies(mock_path(idx), redshift_name='z') 
    v2v.info['MAGLIM'] = -20.0
    v2v.galaxies['rabsmag'] = -21.0
    
    masks = fits.open(masks_path)
    
    v2v_vfc = vc.combined_galaxy_membership(v2v, vfc, masks['SDSS'])
    frac = 100 * v2v_vfc[0]/v2v_vfc[1]
    
    gal_frac[f"VF/VV ({idx})"] = frac
    
def get_overlap_vf_vr(idx):

    vfc = VoidFinderCatalog(vf_path(idx))
    vfc.add_galaxies(mock_path(idx), redshift_name='z') 
    vfc.info['MAGLIM'] = -20.0
    vfc.galaxies['rabsmag'] = -21.0
    
    masks = fits.open(masks_path)

    v2r = V2Catalog(vr_path(idx))
    v2r.add_galaxies(mock_path(idx), redshift_name='z') 
    v2r.info['MAGLIM'] = -20.0
    v2r.galaxies['rabsmag'] = -21.0
    
    v2r_vfc = vc.combined_galaxy_membership(v2r, vfc, masks['SDSS'])
    frac = 100 * v2r_vfc[0]/v2r_vfc[1]
    
    gal_frac[f"VF/VR ({idx})"] = frac
    
def get_overlap_vv_vr(idx):

    v2v = V2Catalog(vv_path(idx))
    v2v.add_galaxies(mock_path(idx), redshift_name='z') 
    v2v.info['MAGLIM'] = -20.0
    v2v.galaxies['rabsmag'] = -21.0

    v2r = V2Catalog(vr_path(idx))
    v2r.add_galaxies(mock_path(idx), redshift_name='z') 
    v2r.info['MAGLIM'] = -20.0
    v2r.galaxies['rabsmag'] = -21.0
    
    masks = fits.open(masks_path)
    
    v2v_v2r = vc.combined_galaxy_membership(v2v, v2r, masks['SDSS'])
    frac = 100 * v2v_v2r[0]/ v2v_v2r[1]
    
    gal_frac[f"VV/VR ({idx})"] = frac
    

#quaternary indexing
idx='000'
idx_list = [idx]
for i in range(24):
    idx = idx[0]+idx[1]+str(int(idx[2]) + 1)
    if int(idx[2]) == 4:
        idx = idx[0] +  str(int(idx[1]) + 1) + '0'

        if int(idx[1]) == 4:
            idx = str(int(idx[0]) + 1) + '00'
    idx_list.append(idx)
    
def save_json():
    with open('gal_frac.json', 'w', encoding='utf-8') as f:
        json.dump(gal_frac.copy(), f, ensure_ascii=False, indent=4)


def run_mock(idx):
    get_overlap_vf(idx) 
    get_overlap_vv(idx)
    get_overlap_vr(idx)
    get_overlap_vf_vv(idx)
    get_overlap_vf_vr(idx)
    get_overlap_vv_vr(idx)
    
print('Launching processes')
processes = []

for idx in idx_list:
    processes.append(Process(target = run_mock, args=[idx]))
    
for p in processes:
    p.start()
    
for p in processes:
    p.join()

gal_frac = gal_frac.copy()

with open('gal_frac.json', 'w', encoding='utf-8') as f:
    json.dump(gal_frac, f, ensure_ascii=False, indent=4)