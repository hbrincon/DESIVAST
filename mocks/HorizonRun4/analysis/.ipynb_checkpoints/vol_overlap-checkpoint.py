import numpy as np
from astropy.table import Table, vstack
from astropy.io import fits
import json
import time
import os

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

if os.path.isfile('void_frac.json'):
    with open('void_frac.json') as f:
        void_frac = json.load(f)
    
    void_frac = manager.dict(void_frac)
else:
    void_frac = manager.dict()


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
    masks = fits.open(masks_path)

    frac = vfc.get_single_overlap( masks['SDSS'])[1] * 100
    
    void_frac[f"VF ({idx})"] = frac
    
def get_overlap_vv(idx):
    
    v2v = V2Catalog(vv_path(idx))
    v2v.add_galaxies(mock_path(idx), redshift_name='z') 
    v2v.info['MAGLIM'] = -20.0
    v2v.galaxies['rabsmag'] = -21.0
    #v2v.add_mask(vfc_ref)
    
    masks = fits.open(masks_path)
    
    frac = v2v.get_single_overlap( masks['SDSS'])[1] * 100
    
    void_frac[f"VV ({idx})"] = frac

    
def get_overlap_vr(idx):

    v2r = V2Catalog(vr_path(idx))
    v2r.add_galaxies(mock_path(idx), redshift_name='z') 
    v2r.info['MAGLIM'] = -20.0
    v2r.galaxies['rabsmag'] = -21.0
    #v2r.add_mask(vfc_ref)
    
    masks = fits.open(masks_path)
    
    frac = v2r.get_single_overlap( masks['SDSS'])[1] * 100
    
    void_frac[f"VR ({idx})"] = frac
    
    
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
    
    v2v_vf = vc.get_overlap(v2v, vfc, 
                            masks['SDSS'])
    
    frac = {}
    frac['Shared'] = v2v_vf[1]/v2v_vf[0]
    frac['Cat 1'] = (v2v_vf[1]+v2v_vf[4])/v2v_vf[0]
    frac['Cat 2'] = (v2v_vf[1]+v2v_vf[3])/v2v_vf[0]
    frac['Num. points'] = int(v2v_vf[0])
    
    void_frac[f"VF/VV ({idx})"] = frac
    
def get_overlap_vf_vr(idx):

    vfc = VoidFinderCatalog(vf_path(idx))
    vfc.add_galaxies(mock_path(idx), redshift_name='z') 
    vfc.info['MAGLIM'] = -20.0
    vfc.galaxies['rabsmag'] = -21.0
    
    v2r = V2Catalog(vr_path(idx))
    v2r.add_galaxies(mock_path(idx), redshift_name='z') 
    v2r.info['MAGLIM'] = -20.0
    v2r.galaxies['rabsmag'] = -21.0
    
    masks = fits.open(masks_path)

    
    v2r_vf = vc.get_overlap(v2r, vfc, 
                            masks['SDSS'])
    
    frac = {}
    frac['Shared'] = v2r_vf[1]/v2r_vf[0]
    frac['Cat 1'] = (v2r_vf[1]+v2r_vf[4])/v2r_vf[0]
    frac['Cat 2'] = (v2r_vf[1]+v2r_vf[3])/v2r_vf[0]
    frac['Num. points'] = int(v2r_vf[0])
    
    void_frac[f"VF/VR ({idx})"] = frac
    
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
    
    v2v_v2r = vc.get_overlap(v2v, v2r, 
                        masks['SDSS'])
        
    frac = {}
    frac['Shared'] = v2v_v2r[1]/v2v_v2r[0]
    frac['Cat 1'] = (v2v_v2r[1]+v2v_v2r[4])/v2v_v2r[0]
    frac['Cat 2'] = (v2v_v2r[1]+v2v_v2r[3])/v2v_v2r[0]
    frac['Num. points'] = int(v2v_v2r[0])
    
    void_frac[f"VV/VR ({idx})"] = frac
    


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
    
processes = []

for idx in idx_list:
    # 20 GB per process * 6 processes = 120 GB 
    # (512 GB total on node)

    if f"VF ({idx})" not in void_frac.keys():
        processes.append(Process(target = get_overlap_vf, args = [idx]))
    if f"VV ({idx})" not in void_frac.keys():
        processes.append(Process(target = get_overlap_vv, args = [idx]))
    if f"VR ({idx})" not in void_frac.keys():
        processes.append(Process(target = get_overlap_vr, args = [idx]))
    if f"VF/VV ({idx})" not in void_frac.keys():
        processes.append(Process(target = get_overlap_vf_vv, args = [idx]))
    if f"VF/VR ({idx})" not in void_frac.keys():
        processes.append(Process(target = get_overlap_vf_vr, args = [idx]))
    if f"VV/VR ({idx})" not in void_frac.keys():
        processes.append(Process(target = get_overlap_vv_vr, args = [idx]))

def save_json():
    with open('void_frac.json', 'w', encoding='utf-8') as f:
        json.dump(void_frac.copy(), f, ensure_ascii=False, indent=4)
                
print('Launching processes')

i = 0

while i*12 < len(processes):

    for p in processes[i*12:i*12+12]:
        p.start()

    for p in processes[i*12:i*12+12]:
        p.join()
        
    i+=1
    save_json()
    print(i)
        

print(f'Output: {len(void_frac)} element dictionary')