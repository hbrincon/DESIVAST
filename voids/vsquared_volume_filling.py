#!/usr/bin/env python
"""
Driver program for Voronoi Voids (V^2) void finder module using the ZOBOV
algorithm.
"""

import configparser, os

from vast.vsquared.zobov import Zobov


config_file ='vsquaredDESI_volume_filling.ini' 


for cap in ('ngc','sgc'):

    print('------------------------')
    print('DESI', cap)
    print('------------------------')

    #Output visualization data
    visualize = True

    parser = configparser.SafeConfigParser()
    parser.read(config_file)

    parser.set('Paths', 'Input Catalog', f'../galaxy_catalog/iron_smoothed_{cap}.fits')
    parser.set('Paths', 'Survey Name', f'DESIVAST_{str.upper(cap)}')

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

    #ZOBOV
    newZobov.sortVoids(method=1)
    newZobov.saveVoids()
    newZobov.saveZones()
    if visualize:
        newZobov.preViz()
