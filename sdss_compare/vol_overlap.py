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

# IMPORTANT NOTE: to run this code, I have manually set the edge_buffer default paramter in python_tools.VoidOverlap.is_edge_point 
# to equal 10 Mpc/h, wheras it usually defaults to 30 Mpc/h. This change must be manually re-made every time this notebook is ran


masks_path = '../galaxy_catalog/mask/alt_masks/masks.fits'

iron_path = '../galaxy_catalog/iron_smoothed_ngc.fits'

nsa_path = ''#'/global/homes/h/hrincon/sdss_compare/nsa_v1_0_1_main.txt'

nsak1_path = '../galaxy_catalog/SDSS_galaxies/nsa_k1.fits'

desi_VF_path = '../voids/data/minus20.0/DESIVAST_NGC_VoidFinder_Output.fits'

desi_VV_path = '../voids/data/minus20.0/DESIVAST_NGC_V2_VIDE_Output.fits'

desi_VR_path = '../voids/data/minus20.0/DESIVAST_NGC_V2_REVOLVER_Output.fits'

sdss_VF_path = ''#'/global/homes/h/hrincon/DESIVAST/alt_mag_cuts/minus20.0/SDSSDR7_VoidFinder_Output.fits'

sdss_VV_path = ''#'/global/homes/h/hrincon/DESIVAST/alt_mag_cuts/minus20.0/SDSSDR7_V2_VIDE_Output.fits'

sdss_VR_path = ''#'/global/homes/h/hrincon/DESIVAST/alt_mag_cuts/minus20.0/SDSSDR7_V2_REVOLVER_Output.fits'

sdssk1_VF_path = '../voids/data/minus20.0/SDSSK1_VoidFinder_Output.fits'

sdssk1_VV_path = '../voids/data/minus20.0/SDSSK1_V2_VIDE_Output.fits'

sdssk1_VR_path = '../voids/data/minus20.0/SDSSK1_V2_REVOLVER_Output.fits'


"""with open('void_frac.json') as f:
    void_frac = json.load(f)
    
manager = Manager()
void_frac = manager.dict(void_frac)"""

manager = Manager()
void_frac = manager.dict()

"""
def get_overlap_vf():

    desi = VoidFinderCatalog(desi_VF_path)
    desi.add_galaxies(iron_path) 
    
    sdss = VoidFinderCatalog(sdss_VF_path)
    sdss.add_galaxies(nsa_path) 
    
    masks = fits.open(masks_path)
    
    overlap = vc.get_overlap(desi, sdss, 
                            masks['COMP'])
        
    frac = {}
    frac['Shared'] = overlap[1]/overlap[0]
    frac['Cat 1'] = overlap[4]/overlap[0]
    frac['Cat 2'] = overlap[3]/overlap[0]
    frac['Num. points'] = int(overlap[0])
    
    void_frac["VF"] = frac
    
def get_overlap_vv():
    
    desi = V2Catalog(desi_VV_path)
    desi.add_galaxies(iron_path) 
    
    sdss = V2Catalog(sdss_VV_path)
    sdss.add_galaxies(nsa_path) 
    
    masks = fits.open(masks_path)
    
    overlap = vc.get_overlap(desi, sdss, 
                            masks['COMP'])
        
    frac = {}
    frac['Shared'] = overlap[1]/overlap[0]
    frac['Cat 1'] = overlap[4]/overlap[0]
    frac['Cat 2'] = overlap[3]/overlap[0]
    frac['Num. points'] = int(overlap[0])
    
    void_frac["VV"] = frac

    
def get_overlap_vr():

    desi = V2Catalog(desi_VR_path)
    desi.add_galaxies(iron_path) 
    
    sdss = V2Catalog(sdss_VR_path)
    sdss.add_galaxies(nsa_path) 
    
    masks = fits.open(masks_path)
    
    overlap = vc.get_overlap(desi, sdss, 
                            masks['COMP'])
        
    frac = {}
    frac['Shared'] = overlap[1]/overlap[0]
    frac['Cat 1'] = overlap[4]/overlap[0]
    frac['Cat 2'] = overlap[3]/overlap[0]
    frac['Num. points'] = int(overlap[0])
    
    void_frac["VR"] = frac
"""    
def get_overlap_k1_vf():

    desi = VoidFinderCatalog(desi_VF_path)
    desi.add_galaxies(iron_path) 
    
    sdss = VoidFinderCatalog(sdssk1_VF_path)
    sdss.add_galaxies(nsak1_path) 
    
    masks = fits.open(masks_path)
    
    overlap = vc.get_overlap(desi, sdss, 
                            masks['COMP'])
        
    frac = {}
    frac['Shared'] = overlap[1]/overlap[0]
    frac['Cat 1'] = overlap[4]/overlap[0]
    frac['Cat 2'] = overlap[3]/overlap[0]
    frac['Num. points'] = int(overlap[0])
    
    void_frac["VF K1"] = frac
    
def get_overlap_k1_vv():
    
    desi = V2Catalog(desi_VV_path)
    desi.add_galaxies(iron_path) 
    
    sdss = V2Catalog(sdssk1_VV_path)
    sdss.add_galaxies(nsak1_path) 
    
    masks = fits.open(masks_path)
    
    overlap = vc.get_overlap(desi, sdss, 
                            masks['COMP'])
        
    frac = {}
    frac['Shared'] = overlap[1]/overlap[0]
    frac['Cat 1'] = overlap[4]/overlap[0] #orders are purposelfully reveresed because vo.report returns catalogs in wrong order
    frac['Cat 2'] = overlap[3]/overlap[0]
    frac['Num. points'] = int(overlap[0])
    
    void_frac["VV K1"] = frac

    
def get_overlap_k1_vr():

    desi = V2Catalog(desi_VR_path)
    desi.add_galaxies(iron_path) 
    
    sdss = V2Catalog(sdssk1_VR_path)
    sdss.add_galaxies(nsak1_path) 
    
    masks = fits.open(masks_path)
    
    overlap = vc.get_overlap(desi, sdss, 
                            masks['COMP'])
        
    frac = {}
    frac['Shared'] = overlap[1]/overlap[0]
    frac['Cat 1'] = overlap[4]/overlap[0] #orders are purposelfully reveresed because vo.report returns catalogs in wrong order
    frac['Cat 2'] = overlap[3]/overlap[0]
    frac['Num. points'] = int(overlap[0])
    
    void_frac["VR K1"] = frac
      


processes = []

#processes.append(Process(target = get_overlap_vf))
#processes.append(Process(target = get_overlap_vv))
#processes.append(Process(target = get_overlap_vr))
processes.append(Process(target = get_overlap_k1_vf))
processes.append(Process(target = get_overlap_k1_vv))
processes.append(Process(target = get_overlap_k1_vr))

"""for  in range(25):
    # 20 GB per process * 6 processes = 120 GB 
    # (512 GB total on node)

    if f"VF ({})" not in void_frac.keys():
        processes.append(Process(target = get_overlap_vf, args = []))
    if f"VV ({})" not in void_frac.keys():
        processes.append(Process(target = get_overlap_vv, args = []))
    if f"VR ({})" not in void_frac.keys():
        processes.append(Process(target = get_overlap_vr, args = []))
    if f"VF/VV ({})" not in void_frac.keys():
        processes.append(Process(target = get_overlap_vf_vv, args = []))
    if f"VF/VR ({})" not in void_frac.keys():
        processes.append(Process(target = get_overlap_vf_vr, args = []))
    if f"VV/VR ({})" not in void_frac.keys():
        processes.append(Process(target = get_overlap_vv_vr, args = []))"""

def save_json():
    with open('void_frac.json', 'w', encoding='utf-8') as f:
        json.dump(void_frac.copy(), f, ensure_ascii=False, indent=4)
                
print('Launching processes')


for p in processes:
    p.start()

for p in processes:
    p.join()

save_json()

        

