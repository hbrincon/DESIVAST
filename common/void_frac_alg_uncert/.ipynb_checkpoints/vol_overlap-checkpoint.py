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
from vast.voidfinder.voidfinder_functions import xyz_to_radecz


masks_path = '/global/homes/h/hrincon/DESIVAST/galaxy_catalog/mask/alt_masks/masks.fits'

iron_ngc_path = '/global/homes/h/hrincon/DESIVAST/galaxy_catalog/iron_ngc.fits'
iron_sgc_path = '/global/homes/h/hrincon/DESIVAST/galaxy_catalog/iron_sgc.fits'
nsa_path = '/global/homes/h/hrincon/sdss_compare/nsa_v1_0_1_main.txt'

manager = Manager()

if os.path.isfile('void_frac.json'):
    with open('void_frac.json') as f:
        void_frac = json.load(f)
    
    void_frac = manager.dict(void_frac)
else:
    void_frac = manager.dict()


def get_overlap_vf(delta):
    
    voids_ngc_path = f'/global/homes/h/hrincon/DESIVAST/VoidFinder/DESIVAST_NGC_VoidFinder_Output.fits'
    voids_sgc_path = f'/global/homes/h/hrincon/DESIVAST/VoidFinder/DESIVAST_SGC_VoidFinder_Output.fits'

    vfc = VoidFinderCatalogStacked(['NGC','SGC'],[voids_ngc_path,voids_sgc_path])
    vfc.add_galaxies([iron_ngc_path,iron_sgc_path]) 
    masks = fits.open(masks_path)
    
    
    tmp_ngc = vfc['NGC'].get_single_overlap( masks['NGCFID'])
    tmp_sgc = vfc['SGC'].get_single_overlap( masks['SGCFID'])

    frac = 100 * (tmp_ngc[0]*tmp_ngc[1]+tmp_sgc[0]*tmp_sgc[1]) / (tmp_ngc[0]+tmp_sgc[0])
    
    void_frac[f"VF ({delta})"] = frac
    
def get_overlap_vv(delta):
    
    
    voids_ngc_path = f'/global/homes/h/hrincon/DESIVAST/V2/DESIVAST_NGC_V2_VIDE_Output.fits'
    voids_sgc_path = f'/global/homes/h/hrincon/DESIVAST/V2/DESIVAST_SGC_V2_VIDE_Output.fits'

    v2v = V2CatalogStacked(['NGC','SGC'],[voids_ngc_path,voids_sgc_path])
    v2v.add_galaxies([iron_ngc_path,iron_sgc_path]) 
    masks = fits.open(masks_path)
    
    tmp_ngc = v2v['NGC'].get_single_overlap( masks['NGCFID'])
    tmp_sgc = v2v['SGC'].get_single_overlap( masks['SGCFID'])

    frac = 100 * (tmp_ngc[0]*tmp_ngc[1]+tmp_sgc[0]*tmp_sgc[1]) / (tmp_ngc[0]+tmp_sgc[0])
    
    void_frac[f"VV ({delta})"] = frac

    
def get_overlap_vr(delta):
    

    voids_ngc_path = f'/global/homes/h/hrincon/DESIVAST/V2/DESIVAST_NGC_V2_REVOLVER_Output.fits'
    voids_sgc_path = f'/global/homes/h/hrincon/DESIVAST/V2/DESIVAST_SGC_V2_REVOLVER_Output.fits'

    v2r = V2CatalogStacked(['NGC','SGC'],[voids_ngc_path,voids_sgc_path])
    v2r.add_galaxies([iron_ngc_path,iron_sgc_path]) 
    masks = fits.open(masks_path)
    
    tmp_ngc = v2r['NGC'].get_single_overlap( masks['NGCFID'])
    tmp_sgc = v2r['SGC'].get_single_overlap( masks['SGCFID'])

    frac = 100 * (tmp_ngc[0]*tmp_ngc[1]+tmp_sgc[0]*tmp_sgc[1]) / (tmp_ngc[0]+tmp_sgc[0])
    
    void_frac[f"VR ({delta})"] = frac
    
"""    
def get_overlap_vf_vv():
    

    voids_ngc_path = f'/global/homes/h/hrincon/DESIVAST/VoidFinder/DESIVAST_NGC_VoidFinder_Output.fits'
    voids_sgc_path = f'/global/homes/h/hrincon/DESIVAST/VoidFinder/DESIVAST_SGC_VoidFinder_Output.fits'

    vfc = VoidFinderCatalogStacked(['NGC','SGC'],[voids_ngc_path,voids_sgc_path])
    vfc.add_galaxies([iron_ngc_path,iron_sgc_path]) 
    
    voids_ngc_path = f'/global/homes/h/hrincon/DESIVAST/V2/DESIVAST_NGC_V2_VIDE_Output.fits'
    voids_sgc_path = f'/global/homes/h/hrincon/DESIVAST/V2/DESIVAST_SGC_V2_VIDE_Output.fits'

    v2v = V2CatalogStacked(['NGC','SGC'],[voids_ngc_path,voids_sgc_path])
    v2v.add_galaxies([iron_ngc_path,iron_sgc_path]) 
    
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
    
    desi_void_frac[f"VF/VV"] = frac
    
def get_overlap_vf_vr():
    

    voids_ngc_path = f'/global/homes/h/hrincon/DESIVAST/VoidFinder/DESIVAST_NGC_VoidFinder_Output.fits'
    voids_sgc_path = f'/global/homes/h/hrincon/DESIVAST/VoidFinder/DESIVAST_SGC_VoidFinder_Output.fits'

    vfc = VoidFinderCatalogStacked(['NGC','SGC'],[voids_ngc_path,voids_sgc_path])
    vfc.add_galaxies([iron_ngc_path,iron_sgc_path]) 

    voids_ngc_path = f'/global/homes/h/hrincon/DESIVAST/V2/DESIVAST_NGC_V2_REVOLVER_Output.fits'
    voids_sgc_path = f'/global/homes/h/hrincon/DESIVAST/V2/DESIVAST_SGC_V2_REVOLVER_Output.fits'
    
    masks = fits.open(masks_path)

    v2r = V2CatalogStacked(['NGC','SGC'],[voids_ngc_path,voids_sgc_path])
    v2r.add_galaxies([iron_ngc_path,iron_sgc_path]) 
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
    
    desi_void_frac[f"VF/VR"] = frac
    
def get_overlap_vv_vr():
    

    voids_ngc_path = f'/global/homes/h/hrincon/DESIVAST/V2/DESIVAST_NGC_V2_VIDE_Output.fits'
    voids_sgc_path = f'/global/homes/h/hrincon/DESIVAST/V2/DESIVAST_SGC_V2_VIDE_Output.fits'

    v2v = V2CatalogStacked(['NGC','SGC'],[voids_ngc_path,voids_sgc_path])
    v2v.add_galaxies([iron_ngc_path,iron_sgc_path]) 
    
    voids_ngc_path = f'/global/homes/h/hrincon/DESIVAST/V2/DESIVAST_NGC_V2_REVOLVER_Output.fits'
    voids_sgc_path = f'/global/homes/h/hrincon/DESIVAST/V2/DESIVAST_SGC_V2_REVOLVER_Output.fits'

    v2r = V2CatalogStacked(['NGC','SGC'],[voids_ngc_path,voids_sgc_path])
    v2r.add_galaxies([iron_ngc_path,iron_sgc_path])   
    
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
    
    desi_void_frac[f"VV/VR"] = frac
    
    
    
    
    
    
    
    
def sdss_get_overlap_vf():

    
    voids_sdss = f'/global/homes/h/hrincon/DESIVAST/alt_mag_cuts/minus20.0/SDSSDR7_VoidFinder_Output.fits'
    vfc = VoidFinderCatalog(voids_sdss)
    vfc.add_galaxies(nsa_path)
    masks = fits.open(masks_path)

    
    frac = vfc.get_single_overlap( masks['SDSS'])[1] * 100
    
    sdss_void_frac[f"VF"] = frac
    
def sdss_get_overlap_vv():
    
    voids_sdss = f'/global/homes/h/hrincon/DESIVAST/alt_mag_cuts/minus20.0/SDSSDR7_V2_VIDE_Output.fits'
    v2v = V2Catalog(voids_sdss)
    v2v.add_galaxies(nsa_path)  
    masks = fits.open(masks_path)
    
    frac = v2v.get_single_overlap( masks['SDSS'])[1] * 100
    
    sdss_void_frac[f"VV"] = frac

    
def sdss_get_overlap_vr():

    voids_sdss = f'/global/homes/h/hrincon/DESIVAST/alt_mag_cuts/minus20.0/SDSSDR7_V2_REVOLVER_Output.fits'
    v2r = V2Catalog(voids_sdss)
    v2r.add_galaxies(nsa_path)
    masks = fits.open(masks_path)
    
    frac = v2r.get_single_overlap( masks['SDSS']) [1] * 100
    
    sdss_void_frac[f"VR"] = frac
    
    
def sdss_get_overlap_vf_vv():
    
    voids_sdss = f'/global/homes/h/hrincon/DESIVAST/alt_mag_cuts/minus20.0/SDSSDR7_VoidFinder_Output.fits'
    vfc = VoidFinderCatalog(voids_sdss)
    vfc.add_galaxies(nsa_path)
    
    voids_sdss = f'/global/homes/h/hrincon/DESIVAST/alt_mag_cuts/minus20.0/SDSSDR7_V2_VIDE_Output.fits'
    v2v = V2Catalog(voids_sdss)
    v2v.add_galaxies(nsa_path) 
    
    masks = fits.open(masks_path)
    
    v2v_vf = vc.get_overlap(v2v, vfc, 
                            masks['SDSS'])
    
    frac = {}
    frac['Shared'] = v2v_vf[1]/v2v_vf[0]
    frac['Cat 1'] = (v2v_vf[1]+v2v_vf[4])/v2v_vf[0]
    frac['Cat 2'] = (v2v_vf[1]+v2v_vf[3])/v2v_vf[0]
    frac['Num. points'] = int(v2v_vf[0])
    
    sdss_void_frac[f"VF/VV"] = frac
    
def sdss_get_overlap_vf_vr():

    voids_sdss = f'/global/homes/h/hrincon/DESIVAST/alt_mag_cuts/minus20.0/SDSSDR7_VoidFinder_Output.fits'
    vfc = VoidFinderCatalog(voids_sdss)
    vfc.add_galaxies(nsa_path)  

    voids_sdss = f'/global/homes/h/hrincon/DESIVAST/alt_mag_cuts/minus20.0/SDSSDR7_V2_REVOLVER_Output.fits'
    v2r = V2Catalog(voids_sdss)
    v2r.add_galaxies(nsa_path)  
    
    masks = fits.open(masks_path)
    
    v2r_vf = vc.get_overlap(v2r, vfc, 
                            masks['SDSS'])
    
    frac = {}
    frac['Shared'] = v2r_vf[1]/v2r_vf[0]
    frac['Cat 1'] = (v2r_vf[1]+v2r_vf[4])/v2r_vf[0]
    frac['Cat 2'] = (v2r_vf[1]+v2r_vf[3])/v2r_vf[0]
    frac['Num. points'] = int(v2r_vf[0])
    
    sdss_void_frac[f"VF/VR"] = frac
    
def sdss_get_overlap_vv_vr():

    voids_sdss = f'/global/homes/h/hrincon/DESIVAST/alt_mag_cuts/minus20.0/SDSSDR7_V2_VIDE_Output.fits'
    v2v = V2Catalog(voids_sdss)
    v2v.add_galaxies(nsa_path)  
    
    voids_sdss = f'/global/homes/h/hrincon/DESIVAST/alt_mag_cuts/minus20.0/SDSSDR7_V2_REVOLVER_Output.fits'
    v2r = V2Catalog(voids_sdss)
    v2r.add_galaxies(nsa_path)  
    
    masks = fits.open(masks_path)
    
    v2v_v2r = vc.get_overlap(v2v, v2r, 
                        masks['SDSS'])
        
    frac = {}
    frac['Shared'] = v2v_v2r[1]/v2v_v2r[0]
    frac['Cat 1'] = (v2v_v2r[1]+v2v_v2r[4])/v2v_v2r[0]
    frac['Cat 2'] = (v2v_v2r[1]+v2v_v2r[3])/v2v_v2r[0]
    frac['Num. points'] = int(v2v_v2r[0])
    
    sdss_void_frac[f"VV/VR"] = frac
    
"""  

def save_json():
    with open('void_frac.json', 'w', encoding='utf-8') as f:
        json.dump(void_frac.copy(), f, ensure_ascii=False, indent=4)

print('Launching processes')
processes = []

for delta in range(25):
    if f"VF ({delta})" not in void_frac.keys():
        processes.append(Process(target = get_overlap_vf, args = [delta]))
    if f"VV ({delta})" not in void_frac.keys():
        processes.append(Process(target = get_overlap_vv, args = [delta]))
    if f"VR ({delta})" not in void_frac.keys():
        processes.append(Process(target = get_overlap_vr, args = [delta]))

"""processes.append(Process(target = get_overlap_vf_vv))
processes.append(Process(target = get_overlap_vf_vr))
processes.append(Process(target = get_overlap_vv_vr))
processes.append(Process(target = sdss_get_overlap_vf))
processes.append(Process(target = sdss_get_overlap_vv))
processes.append(Process(target = sdss_get_overlap_vr))
processes.append(Process(target = sdss_get_overlap_vf_vv))
processes.append(Process(target = sdss_get_overlap_vf_vr))
processes.append(Process(target = sdss_get_overlap_vv_vr))"""
    

i = 0

while i*12 < len(processes):

    for p in processes[i*12:i*12+12]:
        p.start()

    for p in processes[i*12:i*12+12]:
        p.join()
        
    i+=1
    save_json()
    print(i)
        