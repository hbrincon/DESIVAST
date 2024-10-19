#!/usr/bin/env python
"""
Driver program for Voronoi Voids (V^2) void finder module using the ZOBOV
algorithm.
"""

import configparser, os

from vast.vsquared.zobov import Zobov

from multiprocessing import Process


hr4_path = '/global/homes/h/hrincon/DESIVAST/mocks/HorizonRun4/'
config_file_template = hr4_path+'ini/vsquared.ini'


#create list oif quaternary indices
idx='000'
idx_list = [idx]
for i in range(24):
    idx = idx[0]+idx[1]+str(int(idx[2]) + 1)
    if int(idx[2]) == 4:
        idx = idx[0] +  str(int(idx[1]) + 1) + '0'

        if int(idx[1]) == 4:
            idx = str(int(idx[0]) + 1) + '00'
    idx_list.append(idx)
    
    
#voifinding
def find_voids(idx):

    parser = configparser.SafeConfigParser()
    parser.read(config_file_template)

    parser.set('Paths', 'Input Catalog', hr4_path + f'DR7m_{idx}.fits')
    parser.set('Paths', 'Survey Name', f'HR4_{idx}_')
    
    config_file = hr4_path+f'ini/vsquared_{idx}.ini'

    # Writing our configuration file to 'example.ini'
    with open(config_file, 'w') as configfile:
        parser.write(configfile)


    newZobov = Zobov(config_file,
                     start=0,
                     end=3,
                     save_intermediate=False,
                     visualize=True,
                     periodic=False,
                     xyz=False,
                     capitalize_colnames=True)

    #VIDE
    newZobov.sortVoids(method=0)
    newZobov.saveVoids()
    newZobov.saveZones()
    #newZobov.preViz()

    #REVOLVER
    newZobov.sortVoids(method=4)
    newZobov.saveVoids()
    newZobov.saveZones()
    #newZobov.preViz()


print('Launching processes')
processes = []
for idx in idx_list:
    processes.append(Process(target = find_voids, args=[idx]))

for p in processes:
    p.start()
    
for p in processes:
    p.join()
