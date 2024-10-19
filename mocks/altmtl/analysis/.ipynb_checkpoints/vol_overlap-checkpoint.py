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

manager = Manager()

if os.path.isfile('void_frac.json'):
    with open('void_frac.json') as f:
        void_frac = json.load(f)
    
    void_frac = manager.dict(void_frac)
else:
    void_frac = manager.dict()

#mock galaxies
def ngc_gals (idx):
    return f'/global/homes/h/hrincon/DESIVAST/mocks/altmtl/altmtl{idx}_ngc.fits'

def sgc_gals (idx):
    return f'/global/homes/h/hrincon/DESIVAST/mocks/altmtl/altmtl{idx}_sgc.fits'

#VoidFinder
def ngc_vf (idx):
    return f'/global/homes/h/hrincon/DESIVAST/mocks/altmtl/altmtl{idx}_ngc_VoidFinder_Output.fits'

def sgc_vf (idx):
    return f'/global/homes/h/hrincon/DESIVAST/mocks/altmtl/altmtl{idx}_sgc_VoidFinder_Output.fits'

#VIDE
def ngc_vv (idx):
    return f'/global/homes/h/hrincon/DESIVAST/mocks/altmtl/altmtl{idx}_ngc_V2_VIDE_Output.fits'

def sgc_vv (idx):
    return f'/global/homes/h/hrincon/DESIVAST/mocks/altmtl/altmtl{idx}_sgc_V2_VIDE_Output.fits'

#REVOLVER
def ngc_vr (idx):
    return f'/global/homes/h/hrincon/DESIVAST/mocks/altmtl/altmtl{idx}_ngc_V2_REVOLVER_Output.fits'

def sgc_vr (idx):
    return f'/global/homes/h/hrincon/DESIVAST/mocks/altmtl/altmtl{idx}_sgc_V2_REVOLVER_Output.fits'


def get_overlap_vf(idx):

    
    voids_ngc_path = ngc_vf(idx)
    voids_sgc_path = sgc_vf(idx)
    vfc = VoidFinderCatalogStacked(['NGC','SGC'],[voids_ngc_path,voids_sgc_path])
    vfc.add_galaxies([ngc_gals(idx),sgc_gals(idx)]) 
    masks = fits.open(masks_path)

    
    tmp_ngc = vfc['NGC'].get_single_overlap( masks['NGCFID'])
    tmp_sgc = vfc['SGC'].get_single_overlap( masks['SGCFID'])

    frac = 100 * (tmp_ngc[0]*tmp_ngc[1]+tmp_sgc[0]*tmp_sgc[1]) / (tmp_ngc[0]+tmp_sgc[0])
    
    void_frac[f"VF ({idx})"] = frac
    
def get_overlap_vv(idx):
    
    
    voids_ngc_path = ngc_vv(idx)
    voids_sgc_path = sgc_vv(idx)

    v2v = V2CatalogStacked(['NGC','SGC'],[voids_ngc_path,voids_sgc_path])
    v2v.add_galaxies([ngc_gals(idx),sgc_gals(idx)]) 
    masks = fits.open(masks_path)
    
    tmp_ngc = v2v['NGC'].get_single_overlap( masks['NGCFID'])
    tmp_sgc = v2v['SGC'].get_single_overlap( masks['SGCFID'])

    frac = 100 * (tmp_ngc[0]*tmp_ngc[1]+tmp_sgc[0]*tmp_sgc[1]) / (tmp_ngc[0]+tmp_sgc[0])
    
    void_frac[f"VV ({idx})"] = frac

    
def get_overlap_vr(idx):

    voids_ngc_path = ngc_vr(idx)
    voids_sgc_path = sgc_vr(idx)

    v2r = V2CatalogStacked(['NGC','SGC'],[voids_ngc_path,voids_sgc_path])
    v2r.add_galaxies([ngc_gals(idx),sgc_gals(idx)]) 
    masks = fits.open(masks_path)
    
    tmp_ngc = v2r['NGC'].get_single_overlap( masks['NGCFID'])
    tmp_sgc = v2r['SGC'].get_single_overlap( masks['SGCFID'])

    frac = 100 * (tmp_ngc[0]*tmp_ngc[1]+tmp_sgc[0]*tmp_sgc[1]) / (tmp_ngc[0]+tmp_sgc[0])
    
    void_frac[f"VR ({idx})"] = frac
    
    
def get_overlap_vf_vv(idx):
    
    voids_ngc_path = ngc_vf(idx)
    voids_sgc_path = sgc_vf(idx)
    
    vfc = VoidFinderCatalogStacked(['NGC','SGC'],[voids_ngc_path,voids_sgc_path])
    vfc.add_galaxies([ngc_gals(idx),sgc_gals(idx)]) 
    
    voids_ngc_path = ngc_vv(idx)
    voids_sgc_path = sgc_vv(idx)
    
    v2v = V2CatalogStacked(['NGC','SGC'],[voids_ngc_path,voids_sgc_path])
    v2v.add_galaxies([ngc_gals(idx),sgc_gals(idx)]) 
    
    masks = fits.open(masks_path)
    
    v2v_vf_ngc = vc.get_overlap(v2v['NGC'], vfc['NGC'], 
                            masks['NGCFID'])
    v2v_vf_sgc = vc.get_overlap(v2v['SGC'], vfc['SGC'], 
                                masks['SGCFID'])
    res = vc.combine_overlaps([v2v_vf_ngc, v2v_vf_sgc], do_print = False, do_return = True)
    
    frac = {}
    frac['Shared'] = res[0]
    frac['Cat 1'] = res[1]
    frac['Cat 2'] = res[2]
    frac['Num. points'] = int(res[3])
    
    void_frac[f"VF/VV ({idx})"] = frac
    
def get_overlap_vf_vr(idx):
    
    voids_ngc_path = ngc_vf(idx)
    voids_sgc_path = sgc_vf(idx)
    
    vfc = VoidFinderCatalogStacked(['NGC','SGC'],[voids_ngc_path,voids_sgc_path])
    vfc.add_galaxies([ngc_gals(idx),sgc_gals(idx)]) 

    voids_ngc_path = ngc_vr(idx)
    voids_sgc_path = sgc_vr(idx)
    
    masks = fits.open(masks_path)

    v2r = V2CatalogStacked(['NGC','SGC'],[voids_ngc_path,voids_sgc_path])
    v2r.add_galaxies([ngc_gals(idx),sgc_gals(idx)]) 
    v2r_vf_ngc = vc.get_overlap(v2r['NGC'], vfc['NGC'], 
                                masks['NGCFID'])
    v2r_vf_sgc = vc.get_overlap(v2r['SGC'], vfc['SGC'], 
                                masks['SGCFID'])
    res = vc.combine_overlaps([v2r_vf_ngc, v2r_vf_sgc], do_print = False, do_return = True)
    
    frac = {}
    frac['Shared'] = res[0]
    frac['Cat 1'] = res[1]
    frac['Cat 2'] = res[2]
    frac['Num. points'] = int(res[3])
    
    void_frac[f"VF/VR ({idx})"] = frac
    
def get_overlap_vv_vr(idx):

    voids_ngc_path = ngc_vv(idx)
    voids_sgc_path = sgc_vv(idx)
    
    v2v = V2CatalogStacked(['NGC','SGC'],[voids_ngc_path,voids_sgc_path])
    v2v.add_galaxies([ngc_gals(idx),sgc_gals(idx)]) 
    
    voids_ngc_path = ngc_vr(idx)
    voids_sgc_path = sgc_vr(idx)

    v2r = V2CatalogStacked(['NGC','SGC'],[voids_ngc_path,voids_sgc_path])
    v2r.add_galaxies([ngc_gals(idx),sgc_gals(idx)])   
    
    masks = fits.open(masks_path)
    
    v2v_v2r_ngc = vc.get_overlap(v2v['NGC'], v2r['NGC'], 
                                masks['NGCFID'])
    v2v_v2r_sgc = vc.get_overlap(v2v['SGC'], v2r['SGC'], 
                                masks['SGCFID'])
    res = vc.combine_overlaps([v2v_v2r_ngc, v2v_v2r_sgc], do_print = False, do_return = True)
    
    frac = {}
    frac['Shared'] = res[0]
    frac['Cat 1'] = res[1]
    frac['Cat 2'] = res[2]
    frac['Num. points'] = int(res[3])
    
    void_frac[f"VV/VR ({idx})"] = frac
    

processes = []

for idx in range(25):
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