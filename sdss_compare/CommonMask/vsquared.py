#!/usr/bin/env python
"""
Driver program for Voronoi Voids (V^2) void finder module using the ZOBOV
algorithm.
"""

import configparser, os

from vast.vsquared.zobov import Zobov


config_file = 'vsquared.ini'


for survey in ('DESI','SDSS'):
    
    print('--------------------')
    print(survey)
    print('--------------------')
    
    parser = configparser.SafeConfigParser()
    parser.read(config_file)
    
    parser.set('Paths','input catalog', f'./{str.lower(survey)}_fiducial.fits')
    parser.set('Paths','survey name', survey)

    # Writing our configuration file to 'example.ini'
    with open(config_file, 'w') as configfile:
        parser.write(configfile)

    newZobov = Zobov(config_file,
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
    
    mags = None
    if survey == 'DESI': mags = (-19.89, -20.06) 
    if survey == 'SDSS': mags = (-19.94, -20.11)
    
    for mag in mags:


        print('--------------------')
        print(survey, mag)
        print('--------------------')
        
        out_directory = f'./minus{abs(mag)}/'

        parser = configparser.SafeConfigParser()
        parser.read(config_file)

        parser.set('Paths','input catalog', f'./{str.lower(survey)}_fiducial.fits')
        parser.set('Paths','survey name', survey)
        parser.set('Settings', 'rabsmag_min', str(mag))
        parser.set('Paths','Output Directory', out_directory)

        # Writing our configuration file to 'example.ini'
        with open(out_directory + config_file, 'w') as configfile:
            parser.write(configfile)

        newZobov = Zobov(out_directory + config_file,
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