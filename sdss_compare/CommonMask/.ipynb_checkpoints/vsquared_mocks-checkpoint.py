#!/usr/bin/env python
"""
Driver program for Voronoi Voids (V^2) void finder module using the ZOBOV
algorithm.
"""

import configparser, os

from vast.vsquared.zobov import Zobov


config_file = 'vsquared.ini'


for i in range (25):
    
    print('--------------------')
    print('altmtl', i)
    print('--------------------')
    
    out_dir = './altmtl/'
    
    parser = configparser.SafeConfigParser()
    parser.read(config_file)
    
    parser.set('Paths','input catalog', f'./altmtl/altmtl{i}_ngc.fits')
    parser.set('Paths','survey name', f'altmtl{i}')
    parser.set('Paths','output directory', out_dir)

    # Writing our configuration file to 'example.ini'
    with open(out_dir + config_file, 'w') as configfile:
        parser.write(configfile)

    newZobov = Zobov(out_dir + config_file,
                     start=0,
                     end=3,
                     save_intermediate=False,
                     visualize=True,
                     #visualize=False,
                     periodic=False,
                     xyz=False,
                     capitalize_colnames=True)

    #VIDE
    newZobov.sortVoids(method=0)
    newZobov.saveVoids()
    newZobov.saveZones()
    newZobov.preViz()

    #REVOLVER
    newZobov.sortVoids(method=4)
    newZobov.saveVoids()
    newZobov.saveZones()
    #newZobov.preViz()

idx='000'
idx_list = [idx]
for i in range(24):
    idx = idx[0]+idx[1]+str(int(idx[2]) + 1)
    if int(idx[2]) == 4:
        idx = idx[0] +  str(int(idx[1]) + 1) + '0'

        if int(idx[1]) == 4:
            idx = str(int(idx[0]) + 1) + '00'
    idx_list.append(idx)
    
for i in idx_list:
    
    print('--------------------')
    print('HR4', i)
    print('--------------------')
    
    out_dir = './HR4/'
    
    parser = configparser.SafeConfigParser()
    parser.read(config_file)
    
    parser.set('Paths','input catalog', f'./HR4/DR7m_{i}.fits')
    parser.set('Paths','survey name', f'HR4_{i}')
    parser.set('Paths','output directory', out_dir)
    parser.set('Settings','rabsmag_min', 'None')

    # Writing our configuration file to 'example.ini'
    with open(out_dir + config_file, 'w') as configfile:
        parser.write(configfile)

    newZobov = Zobov(out_dir + config_file,
                     start=0,
                     end=3,
                     save_intermediate=False,
                     visualize=True,
                     #visualize=False,
                     periodic=False,
                     xyz=False,
                     capitalize_colnames=True)

    #VIDE
    newZobov.sortVoids(method=0)
    newZobov.saveVoids()
    newZobov.saveZones()
    newZobov.preViz()

    #REVOLVER
    newZobov.sortVoids(method=4)
    newZobov.saveVoids()
    newZobov.saveZones()
    #newZobov.preViz()