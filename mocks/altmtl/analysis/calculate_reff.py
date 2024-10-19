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

#mock galaxies
def ngc_gals (idx):
    return f'/global/homes/h/hrincon/DESIVAST/mocks/altmtl/altmtl{idx}_ngc.fits'

def sgc_gals (idx):
    return f'/global/homes/h/hrincon/DESIVAST/mocks/altmtl/altmtl{idx}_sgc.fits'

def sdss_gals (idx):
    idx = idx_list[idx] #convert idx from decimal to base 4
    return f'/global/homes/h/hrincon/DESIVAST/mocks/HorizonRun4/DR7m_{idx}.fits'

#VoidFinder
def ngc_vf (idx):
    return f'/global/homes/h/hrincon/DESIVAST/mocks/altmtl/altmtl{idx}_ngc_VoidFinder_Output.fits'

def sgc_vf (idx):
    return f'/global/homes/h/hrincon/DESIVAST/mocks/altmtl/altmtl{idx}_sgc_VoidFinder_Output.fits'

def sdss_vf(idx):
    idx = idx_list[idx] #convert idx from decimal to base 4
    return f'/global/homes/h/hrincon/DESIVAST/mocks/HorizonRun4/HR4_{idx}_VoidFinder_Output.fits'

def get_overlap_vf(idx):

    
    voids_ngc_path = ngc_vf(idx)
    voids_sgc_path = sgc_vf(idx)
    vfc = VoidFinderCatalogStacked(['NGC','SGC'],[voids_ngc_path,voids_sgc_path])
    vfc.add_galaxies([ngc_gals(idx),sgc_gals(idx)]) 

    vfc.calculate_r_eff()

    
def sdss_get_overlap_vf(idx):

    
    voids_path = sdss_vf(idx)
    vfc = VoidFinderCatalog(voids_path)
    vfc.add_galaxies(sdss_gals(idx), redshift_name='z') 

    vfc.calculate_r_eff()
    
processes = []

for idx in range(25):

    processes.append(Process(target = get_overlap_vf, args = [idx]))
    processes.append(Process(target = sdss_get_overlap_vf, args = [idx]))
    
       
print('Launching processes')

for p in processes:
    p.start()

for p in processes:
    p.join()
        