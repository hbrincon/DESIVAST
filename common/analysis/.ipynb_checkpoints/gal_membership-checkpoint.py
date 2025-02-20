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

manager = Manager()

if os.path.isfile('gal_frac.json'):
    with open('gal_frac.json') as f:
        gal_frac = manager.dict(json.load(f))
else:
    gal_frac = manager.dict()
    
  
       
def get_single_overlap(survey_name, voids_path, galaxies, mask, edge_buffer):
    
    if 'VIDE' in voids_path:
        label = 'VV'
        catalog = V2Catalog(voids_path, edge_buffer=edge_buffer)
    elif 'REVOLVER' in voids_path:
        label = 'VR'
        catalog = V2Catalog(voids_path, edge_buffer=edge_buffer)
    elif 'VoidFinder' in voids_path:
        label = 'VF'
        catalog = VoidFinderCatalog(voids_path, edge_buffer=edge_buffer)

    catalog.add_galaxies(galaxies)
    
    overlap = catalog.galaxy_membership(mask)
    
    frac = {}
    frac['Num. void'] = int(overlap[0])
    frac['Num. total'] = int(overlap[1])
    
    mag = catalog.info['MAGLIM']
        
    gal_frac[f"{survey_name} : {label} ({str(mag)}) "] = frac
    
    
def get_overlap(survey_name, voids_path_1, voids_path_2, galaxies, mask, edge_buffer):
    
    if 'VIDE' in voids_path_1:
        label_1 = 'VV'
        catalog_1 = V2Catalog(voids_path_1, edge_buffer=edge_buffer)
    elif 'REVOLVER' in voids_path_1:
        label_1 = 'VR'
        catalog_1 = V2Catalog(voids_path_1, edge_buffer=edge_buffer)
    elif 'VoidFinder' in voids_path_1:
        label_1 = 'VF'
        catalog_1 = VoidFinderCatalog(voids_path_1, edge_buffer=edge_buffer)
        
    if 'VIDE' in voids_path_2:
        label_2 = 'VV'
        catalog_2 = V2Catalog(voids_path_2, edge_buffer=edge_buffer)
    elif 'REVOLVER' in voids_path_2:
        label_2 = 'VR'
        catalog_2 = V2Catalog(voids_path_2, edge_buffer=edge_buffer)
    elif 'VoidFinder' in voids_path_2:
        label_2 = 'VF'
        catalog_2 = VoidFinderCatalog(voids_path_2, edge_buffer=edge_buffer)
    
    assert catalog_1.info['MAGLIM'] == catalog_2.info['MAGLIM']
    
    catalog_1.add_galaxies(galaxies)
    catalog_2.add_galaxies(galaxies)
    
    
    overlap = vc.combined_galaxy_membership(catalog_1, catalog_2, mask)
    
    frac = {}
    frac['Num. void'] = int(overlap[0])
    frac['Num. total'] = int(overlap[1])
    
    mag = catalog_1.info['MAGLIM']
    
    gal_frac[f"{survey_name} : {label_1}/{label_2} ({mag}) "] = frac
    

print('Calculating Void Properity Fractions')

masks = fits.open('../../galaxy_catalog/mask/alt_masks/masks.fits')
sdss_galaxies = '../../galaxy_catalog/SDSS_galaxies/nsa_k1.fits'
ngc_galaxies = '../../galaxy_catalog/iron_smoothed_ngc.fits'
sgc_galaxies = '../../galaxy_catalog/iron_smoothed_sgc.fits'

def catalog_path(survey_name, algorithm, magnitude):
    return 

print('Launching processes')
processes = []

edge_buffer=30#Mpc/h

# SDSS
for magnitude in [-19.94, -20.0, -20.11]:
    void_finder = f'../../voids/data/minus{str(np.abs(magnitude))}/SDSSK1_VoidFinder_Output.fits'
    vide = f'../../voids/data/minus{str(np.abs(magnitude))}/SDSSK1_V2_VIDE_Output.fits'
    revolver = f'../../voids/data/minus{str(np.abs(magnitude))}/SDSSK1_V2_REVOLVER_Output.fits'
    
    processes.append(Process(target = get_single_overlap,args=['SDSS', void_finder, sdss_galaxies, masks['SDSS'], edge_buffer]))
    processes.append(Process(target = get_single_overlap,args=['SDSS', vide, sdss_galaxies, masks['SDSS'], edge_buffer]))
    processes.append(Process(target = get_single_overlap,args=['SDSS', revolver, sdss_galaxies, masks['SDSS'], edge_buffer]))
    processes.append(Process(target = get_overlap,args=['SDSS', void_finder, vide, sdss_galaxies, masks['SDSS'], edge_buffer]))
    processes.append(Process(target = get_overlap,args=['SDSS', void_finder, revolver, sdss_galaxies, masks['SDSS'], edge_buffer]))
    processes.append(Process(target = get_overlap,args=['SDSS', vide, revolver, sdss_galaxies, masks['SDSS'], edge_buffer]))

#DESI
for magnitude in [19.89, -20.0, -20.06]:
    
    void_finder = f'../../voids/data/minus{str(np.abs(magnitude))}/DESIVAST_NGC_VoidFinder_Output.fits'
    vide = f'../../voids/data/minus{str(np.abs(magnitude))}/DESIVAST_NGC_V2_VIDE_Output.fits'
    revolver = f'../../voids/data/minus{str(np.abs(magnitude))}/DESIVAST_NGC_V2_REVOLVER_Output.fits'
    
    processes.append(Process(target = get_single_overlap,args=['NGC', void_finder, ngc_galaxies, masks['NGC'], edge_buffer]))
    processes.append(Process(target = get_single_overlap,args=['NGC', vide, ngc_galaxies, masks['NGC'], edge_buffer]))
    processes.append(Process(target = get_single_overlap,args=['NGC', revolver, ngc_galaxies, masks['NGC'], edge_buffer]))
    processes.append(Process(target = get_overlap,args=['NGC', void_finder, vide, ngc_galaxies, masks['NGC'], edge_buffer]))
    processes.append(Process(target = get_overlap,args=['NGC', void_finder, revolver, ngc_galaxies, masks['NGC'], edge_buffer]))
    processes.append(Process(target = get_overlap,args=['NGC', vide, revolver, ngc_galaxies, masks['NGC'], edge_buffer]))

    
    void_finder = f'../../voids/data/minus{str(np.abs(magnitude))}/DESIVAST_SGC_VoidFinder_Output.fits'
    vide = f'../../voids/data/minus{str(np.abs(magnitude))}/DESIVAST_SGC_V2_VIDE_Output.fits'
    revolver = f'../../voids/data/minus{str(np.abs(magnitude))}/DESIVAST_SGC_V2_REVOLVER_Output.fits'
    
    processes.append(Process(target = get_single_overlap,args=['SGC', void_finder, sgc_galaxies, masks['SGC'], edge_buffer]))
    processes.append(Process(target = get_single_overlap,args=['SGC', vide, sgc_galaxies, masks['SGC'], edge_buffer]))
    processes.append(Process(target = get_single_overlap,args=['SGC', revolver, sgc_galaxies, masks['SGC'], edge_buffer]))
    processes.append(Process(target = get_overlap,args=['SGC', void_finder, vide, sgc_galaxies, masks['SGC'], edge_buffer]))
    processes.append(Process(target = get_overlap,args=['SGC', void_finder, revolver, sgc_galaxies, masks['SGC'], edge_buffer]))
    processes.append(Process(target = get_overlap,args=['SGC', vide, revolver, sgc_galaxies, masks['SGC'], edge_buffer]))


for p in processes:
    p.start()

for p in processes:
    p.join()


gal_frac = gal_frac.copy()

with open('gal_frac.json', 'w', encoding='utf-8') as f:
    json.dump(gal_frac, f, ensure_ascii=False, indent=4)








