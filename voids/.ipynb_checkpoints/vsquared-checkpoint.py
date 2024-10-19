#!/usr/bin/env python
"""
Driver program for Voronoi Voids (V^2) void finder module using the ZOBOV
algorithm.
"""

import configparser, os

from vast.vsquared.zobov import Zobov


config_file = 'vsquaredSDSS.ini'


for mag in (-19.94, -20.0, -20.11):
    
    print('------------------------')
    print('SDSS:', mag)
    print('------------------------')
    
    parser = configparser.SafeConfigParser()
    parser.read(config_file)
    
    parser.set('Settings', 'rabsmag_min', str(mag))
    parser.set('Paths','Output Directory', f'./data/minus{abs(mag)}/')

    # Writing our configuration file to 'example.ini'
    with open(config_file, 'w') as configfile:
        parser.write(configfile)
        
    if mag == -20.0:
        visualize = True
    else:
        visualize = False

    newZobov = Zobov(config_file,
                     start=0,
                     end=3,
                     save_intermediate=False,
                     visualize=visualize,
                     periodic=False,
                     xyz=False,
                     capitalize_colnames=True)

    #VIDE
    newZobov.sortVoids(method=0)
    newZobov.saveVoids()
    newZobov.saveZones()
    if visualize:
        newZobov.preViz()

    #REVOLVER
    newZobov.sortVoids(method=4)
    newZobov.saveVoids()
    newZobov.saveZones()
    if visualize:
        newZobov.preViz()


config_file ='vsquaredDESI.ini' 

for mag in (-19.89, -20.0, -20.06):
    
    for cap in ('ngc','sgc'):
        
        print('------------------------')
        print('DESI', cap, ':', mag)
        print('------------------------')
        
        if mag == -20.0:
            visualize = True
        else:
            visualize = False
        
        parser = configparser.SafeConfigParser()
        parser.read(config_file)

        parser.set('Paths', 'Input Catalog', f'../galaxy_catalog/iron_smoothed_{cap}.fits')
        parser.set('Paths', 'Survey Name', f'DESIVAST_{str.upper(cap)}')

        parser.set('Settings', 'rabsmag_min', str(mag))
        parser.set('Paths','Output Directory', f'./data/minus{abs(mag)}/')

        # Writing our configuration file to 'example.ini'
        with open(config_file, 'w') as configfile:
            parser.write(configfile)

        newZobov = Zobov(config_file,
                         start=0,
                         end=3,
                         save_intermediate=False,
                         visualize=visualize,
                         periodic=False,
                         xyz=False,
                         capitalize_colnames=True)

        #VIDE
        newZobov.sortVoids(method=0)
        newZobov.saveVoids()
        newZobov.saveZones()
        if visualize:
            newZobov.preViz()

        #REVOLVER
        newZobov.sortVoids(method=4)
        newZobov.saveVoids()
        newZobov.saveZones()
        if visualize:
            newZobov.preViz()
