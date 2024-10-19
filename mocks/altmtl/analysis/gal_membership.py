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
manager = Manager()
gal_frac = manager.dict()

def mock_ngc_path(idx):
    return f'/global/homes/h/hrincon/DESIVAST/mocks/altmtl/altmtl{idx}_ngc.fits'
def mock_sgc_path(idx):
    return f'/global/homes/h/hrincon/DESIVAST/mocks/altmtl/altmtl{idx}_sgc.fits'
def vf_ngc_path(idx):
    return f'/global/homes/h/hrincon/DESIVAST/mocks/altmtl/altmtl{idx}_ngc_VoidFinder_Output.fits'
def vf_sgc_path(idx):
    return f'/global/homes/h/hrincon/DESIVAST/mocks/altmtl/altmtl{idx}_sgc_VoidFinder_Output.fits'
def vv_ngc_path(idx):
    return f'/global/homes/h/hrincon/DESIVAST/mocks/altmtl/altmtl{idx}_ngc_V2_VIDE_Output.fits'
def vv_sgc_path(idx):
    return f'/global/homes/h/hrincon/DESIVAST/mocks/altmtl/altmtl{idx}_sgc_V2_VIDE_Output.fits'
def vr_ngc_path(idx):
    return f'/global/homes/h/hrincon/DESIVAST/mocks/altmtl/altmtl{idx}_ngc_V2_REVOLVER_Output.fits'
def vr_sgc_path(idx):
    return f'/global/homes/h/hrincon/DESIVAST/mocks/altmtl/altmtl{idx}_sgc_V2_REVOLVER_Output.fits'


def get_overlap_vf(idx):

    vfc = VoidFinderCatalogStacked(['NGC','SGC'],[vf_ngc_path(idx),vf_sgc_path(idx)])
    vfc['NGC'].info['MAGLIM'] = -20.0
    vfc['SGC'].info['MAGLIM'] = -20.0
    vfc.add_galaxies([mock_ngc_path(idx),mock_sgc_path(idx)]) 
    vfc.calculate_vflag()
    masks = fits.open(masks_path)

    tmp_ngc = vfc['NGC'].galaxy_membership(masks['NGCFID'])
    tmp_sgc = vfc['SGC'].galaxy_membership(masks['SGCFID'])    

    frac = 100 * (tmp_ngc[0]+tmp_sgc[0]) / (tmp_ngc[1]+tmp_sgc[1])
    
    gal_frac[f"VF ({idx})"] = frac
    
def get_overlap_vv(idx):
    
    v2v = V2CatalogStacked(['NGC','SGC'],[vv_ngc_path(idx),vv_sgc_path(idx)])
    v2v.add_galaxies([mock_ngc_path(idx),mock_sgc_path(idx)]) 
    #v2v.add_mask(vfc_ref)
    
    masks = fits.open(masks_path)
    
    tmp_ngc = v2v['NGC'].galaxy_membership(masks['NGCFID'])
    tmp_sgc = v2v['SGC'].galaxy_membership(masks['SGCFID'])    

    frac = 100 * (tmp_ngc[0]+tmp_sgc[0]) / (tmp_ngc[1]+tmp_sgc[1])   
    
    gal_frac[f"VV ({idx})"] = frac

    
def get_overlap_vr(idx):

    v2r = V2CatalogStacked(['NGC','SGC'],[vr_ngc_path(idx),vr_sgc_path(idx)])
    v2r.add_galaxies([mock_ngc_path(idx),mock_sgc_path(idx)]) 
    #v2r.add_mask(vfc_ref)
    
    masks = fits.open(masks_path)
    
    tmp_ngc = v2r['NGC'].galaxy_membership(masks['NGCFID'])
    tmp_sgc = v2r['SGC'].galaxy_membership(masks['SGCFID'])    

    frac = 100 * (tmp_ngc[0]+tmp_sgc[0]) / (tmp_ngc[1]+tmp_sgc[1])
    
    gal_frac[f"VR ({idx})"] = frac
    
    
def get_overlap_vf_vv(idx):

    vfc = VoidFinderCatalogStacked(['NGC','SGC'],[vf_ngc_path(idx),vf_sgc_path(idx)])
    vfc['NGC'].info['MAGLIM'] = -20.0
    vfc['SGC'].info['MAGLIM'] = -20.0
    vfc.add_galaxies([mock_ngc_path(idx),mock_sgc_path(idx)]) 

    v2v = V2CatalogStacked(['NGC','SGC'],[vv_ngc_path(idx),vv_sgc_path(idx)])
    v2v.add_galaxies([mock_ngc_path(idx),mock_sgc_path(idx)]) 
    #v2v.add_mask(vfc_ref)
    
    masks = fits.open(masks_path)
    
    v2v_vfc_ngc = vc.combined_galaxy_membership(v2v['NGC'], vfc['NGC'], masks['NGCFID'])
    v2v_vfc_sgc = vc.combined_galaxy_membership(v2v['SGC'], vfc['SGC'], masks['SGCFID'])

    frac = 100 * (v2v_vfc_ngc[0]+v2v_vfc_sgc[0]) / (v2v_vfc_ngc[1]+v2v_vfc_sgc[1])
    
    gal_frac[f"VF/VV ({idx})"] = frac
    
def get_overlap_vf_vr(idx):

    vfc = VoidFinderCatalogStacked(['NGC','SGC'],[vf_ngc_path(idx),vf_sgc_path(idx)])
    vfc['NGC'].info['MAGLIM'] = -20.0
    vfc['SGC'].info['MAGLIM'] = -20.0
    vfc.add_galaxies([mock_ngc_path(idx),mock_sgc_path(idx)]) 
    
    masks = fits.open(masks_path)

    v2r = V2CatalogStacked(['NGC','SGC'],[vr_ngc_path(idx),vr_sgc_path(idx)])
    v2r.add_galaxies([mock_ngc_path(idx),mock_sgc_path(idx)]) 
    #v2r.add_mask(vfc_ref)
    
    v2r_vfc_ngc = vc.combined_galaxy_membership(v2r['NGC'], vfc['NGC'], masks['NGCFID'])
    v2r_vfc_sgc = vc.combined_galaxy_membership(v2r['SGC'], vfc['SGC'], masks['SGCFID'])

    frac = 100 * (v2r_vfc_ngc[0]+v2r_vfc_sgc[0]) / (v2r_vfc_ngc[1]+v2r_vfc_sgc[1])
    
    gal_frac[f"VF/VR ({idx})"] = frac
    
def get_overlap_vv_vr(idx):

    v2v = V2CatalogStacked(['NGC','SGC'],[vv_ngc_path(idx),vv_sgc_path(idx)])
    v2v.add_galaxies([mock_ngc_path(idx),mock_sgc_path(idx)]) 

    v2r = V2CatalogStacked(['NGC','SGC'],[vr_ngc_path(idx),vr_sgc_path(idx)])
    v2r.add_galaxies([mock_ngc_path(idx),mock_sgc_path(idx)])   
    
    masks = fits.open(masks_path)
    
    v2v_v2r_ngc = vc.combined_galaxy_membership(v2v['NGC'], v2r['NGC'], masks['NGCFID'])
    v2v_v2r_sgc = vc.combined_galaxy_membership(v2v['SGC'], v2r['SGC'], masks['SGCFID'])

    frac = 100 * (v2v_v2r_ngc[0]+v2v_v2r_sgc[0]) / (v2v_v2r_ngc[1]+v2v_v2r_sgc[1])

    
    gal_frac[f"VV/VR ({idx})"] = frac
    


print('Launching processes')
processes = []
for idx in range(25):
    processes.append(Process(target = get_overlap_vf,args=[idx]))
    
for p in processes:
    p.start()
    
for p in processes:
    p.join()

for idx in range(25):
    #get_overlap_vf(idx) # this is done in parallel now
    get_overlap_vv(idx)
    get_overlap_vr(idx)
    get_overlap_vf_vv(idx)
    get_overlap_vf_vr(idx)
    get_overlap_vv_vr(idx)

gal_frac = gal_frac.copy()

with open('gal_frac.json', 'w', encoding='utf-8') as f:
    json.dump(gal_frac, f, ensure_ascii=False, indent=4)