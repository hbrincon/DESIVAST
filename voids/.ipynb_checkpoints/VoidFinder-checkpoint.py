################################################################################
# VoidFinder - Hoyle & Vogeley (2002)
#
# This is where VoidFinder is run on the DESI BGS Bright
# galaxy catalog.
################################################################################




################################################################################
# IMPORT MODULES
#
# If you have control over your python environment, voidfinder can be installed
# as a normal python package via 'python setup.py install', in which case the 
# below import of 'sys' and 'sys.path.insert(0, '/abspath/to/VoidFinder/python'
# is unnecessary.  If you aren't able to install the voidfinder package,
# you can use the sys.path.insert to add it to the list of available packages
# in your python environment.
#
# Alternately, "python setup.py develop" will 'install' some symlinks which
# point back to the current directory and you can run off the same voidfinder
# repository that you're working on as if it was installed
#-------------------------------------------------------------------------------
#import sys
#sys.path.insert(1, 'local/path/VAST/VoidFinder/vast/voidfinder/')

from vast.voidfinder import find_voids, filter_galaxies

from vast.voidfinder.multizmask import generate_mask
from vast.voidfinder.preprocessing import file_preprocess

from vast.voidfinder.constants import c

import pickle
import numpy as np

import os
import sys
sys.path.insert(1, '/global/homes/h/hrincon/python_tools')
import VoidVolume as vol
import VoidOverlap as vo
import VoidCatalog as vc
import VoidSlicePlots as vsp

from vast.voidfinder.postprocessing import mknum

from multiprocessing import Process, Manager

################################################################################

# The first step is set up for 1 process, while the second step is parellizied
find_void_catalogs = False
calc_void_properties = True


################################################################################
# USER INPUTS
#-------------------------------------------------------------------------------
# Number of CPUs available for analysis.
# A value of None will use one less than all available CPUs.
num_cpus = 1

#-------------------------------------------------------------------------------
# File name details
#-------------------------------------------------------------------------------
# File header

survey_names = ['SDSSK1_']


# Change these directory paths to where your data is stored, and where you want 
# the output to be saved.
in_directory = '../galaxy_catalog/SDSS_galaxies/'


# Input file name
# File format: RA, dec, redshift, comoving distance, absolute magnitude
galaxies_filenames = ['nsa_k1.fits']

sdss_galaxies = in_directory + galaxies_filenames[0]
sdss_voids = []
#-------------------------------------------------------------------------------


#-------------------------------------------------------------------------------
# Survey parameters
#-------------------------------------------------------------------------------
# Redshift limits
# Note: These can be set to None, in which case VoidFinder will use the limits 
# of the galaxy catalog.
min_z = 0.
max_z = 0.114


# Cosmology (uncomment and change values to change cosmology)
# Need to also uncomment relevent inputs in function calls below
Omega_M = 0.315
#h = 1


# Uncomment if you do NOT want to use comoving distances
# Need to also uncomment relevent inputs in function calls below
dist_metric = 'comoving'
#-------------------------------------------------------------------------------


#-------------------------------------------------------------------------------
# Galaxy pruning details
#-------------------------------------------------------------------------------
# Uncomment if you do NOT want to remove galaxies with Mr > -20
# Need to also uncomment relevent input in function calls below
#mag_cut = False
magnitude_limits = [-19.94, -20.0, -20.11]

# Uncomment if you do NOT want to remove isolated galaxies
# Need to also uncomment relevent input in function calls below
#rm_isolated = False
#-------------------------------------------------------------------------------
################################################################################

for magnitude_limit in magnitude_limits:
    #print(np.abs(magnitude_limit))
    out_directory = f'./data/minus{np.abs(magnitude_limit)}/'  

    for survey_name, galaxies_filename in zip (survey_names, galaxies_filenames):
        sdss_voids.append(out_directory+survey_name+'VoidFinder_Output.fits')
        
        if find_void_catalogs:
            print('------------------------')
            print(survey_name, magnitude_limit)
            print('------------------------')

            ################################################################################
            # PREPROCESS DATA
            #-------------------------------------------------------------------------------
            galaxy_data_table, dist_limits = file_preprocess(
                galaxies_filename, 
                survey_name,
                in_directory, 
                out_directory, 
                #mag_cut=mag_cut,
                dist_metric=dist_metric,
                min_z=min_z, 
                max_z=max_z,
                Omega_M=Omega_M,
                #h=h,
            )

            print("Dist limits: ", dist_limits)
            print("Number of galaxies read in: ", len(galaxy_data_table))
            ################################################################################




            ################################################################################
            # GENERATE MASK
            #-------------------------------------------------------------------------------
            mask, mask_resolution = generate_mask(galaxy_data_table, 
                                                  max_z, 
                                                  survey_name,
                                                  out_directory,
                                                  dist_metric=dist_metric, 
                                                  smooth_mask=True,
                                                  #h=h,
                                                  )

            # The mask is automatically saved in the survey_name+'VoidFinder_Output.fits'
            # file and can be read in for future use
            ################################################################################




            ################################################################################
            # FILTER GALAXIES
            #-------------------------------------------------------------------------------
            # If you are rerunning the code, you can comment out the mask generation step 
            # above and just load it here instead. Use something in the vein of the below 
            # (untested) code:
            # with fits.open(out_directory+survey_name+'_VoidFinder_Output.fits') as output:
            #   mask = output['MASK'].data
            #   mask_resolution  = output['MASK'].header['MSKRES']
            #   dist_limits = (output['PRIMARY'].header['DLIML'], output['PRIMARY'].header['DLIMU'])

            wall_coords_xyz, field_coords_xyz = filter_galaxies(galaxy_data_table,
                                                                survey_name,
                                                                out_directory,
                                                                dist_limits=dist_limits,
                                                                magnitude_limit=magnitude_limit,
                                                                #mag_cut_flag=mag_cut,
                                                                #rm_isolated_flag=rm_isolated,
                                                                dist_metric=dist_metric,
                                                                capitalize_colnames=True,
                                                                #h=h,
                                                                )
            print("Number of filtred galaxies: ", len(wall_coords_xyz)+len(field_coords_xyz))
            del galaxy_data_table
            # The galaxies are automatically saved in the survey_name+'VoidFinder_Output.fits'
            # file and can be read in for future use, o long as the write_table parameter of
            # filter_galaxies is True (it is True by default)
            ################################################################################


            coords_min = np.min(np.concatenate([wall_coords_xyz, field_coords_xyz]), axis=0)


            ################################################################################
            # FIND VOIDS
            #-------------------------------------------------------------------------------
            # Again, if you are running the code and have not changed any of the above steps 
            # from a previous run, you can comment out most of the above function calls and 
            # load all the details in here to start over. Use something in the vein of the 
            # below (untested) code:
            # with fits.open(out_directory+survey_name+'_VoidFinder_Output.fits') as output:
            #   wall_coords_xyz = np.array([output['WALL'].data['X'],
            #                               output['WALL'].data['Y'],
            #                               output['WALL'].data['Z']]).T



            find_voids(wall_coords_xyz, 
                       survey_name,
                       out_directory,
                       mask_type='ra_dec_z',
                       mask=mask, 
                       mask_resolution=mask_resolution,
                       dist_limits=dist_limits,
                       grid_origin=coords_min,
                       #save_after=50000,
                       #use_start_checkpoint=True,
                       verbose=1,
                       num_cpus=num_cpus,
                       capitalize_colnames=True)

        


################################################################################
# USER INPUTS
#-------------------------------------------------------------------------------
# Number of CPUs available for analysis.
# A value of None will use one less than all available CPUs.
num_cpus = 1

#-------------------------------------------------------------------------------
# File name details
#-------------------------------------------------------------------------------
# File header

survey_names = ['DESIVAST_NGC_','DESIVAST_SGC_']


# Change these directory paths to where your data is stored, and where you want 
# the output to be saved.
in_directory = '../galaxy_catalog/'


# Input file name
# File format: RA, dec, redshift, comoving distance, absolute magnitude
galaxies_filenames = ['iron_smoothed_ngc.fits', 'iron_smoothed_sgc.fits']

ngc_galaxies = in_directory + galaxies_filenames[0]
sgc_galaxies = in_directory + galaxies_filenames[1]
ngc_voids = []
sgc_voids = []
#-------------------------------------------------------------------------------


#-------------------------------------------------------------------------------
# Survey parameters
#-------------------------------------------------------------------------------
# Redshift limits
# Note: These can be set to None, in which case VoidFinder will use the limits 
# of the galaxy catalog.
min_z = 0.
max_z = 0.24


# Cosmology (uncomment and change values to change cosmology)
# Need to also uncomment relevent inputs in function calls below
Omega_M = 0.315
#h = 1


# Uncomment if you do NOT want to use comoving distances
# Need to also uncomment relevent inputs in function calls below
dist_metric = 'comoving'
#-------------------------------------------------------------------------------


#-------------------------------------------------------------------------------
# Galaxy pruning details
#-------------------------------------------------------------------------------
# Uncomment if you do NOT want to remove galaxies with Mr > -20
# Need to also uncomment relevent input in function calls below
#mag_cut = False
magnitude_limits = [-19.89, -20.0, -20.06]


# Uncomment if you do NOT want to remove isolated galaxies
# Need to also uncomment relevent input in function calls below
#rm_isolated = False
#-------------------------------------------------------------------------------
################################################################################

for magnitude_limit in magnitude_limits:
    
    out_directory = f'./data/minus{np.abs(magnitude_limit)}/'
    

    for survey_name, galaxies_filename, voids_list in zip (survey_names, galaxies_filenames, (ngc_voids, sgc_voids)):
        voids_list.append(out_directory+survey_name+'VoidFinder_Output.fits')
        
        if find_void_catalogs:
        
            print('------------------------')
            print(survey_name)
            print('------------------------')

            ################################################################################
            # PREPROCESS DATA
            #-------------------------------------------------------------------------------
            galaxy_data_table, dist_limits = file_preprocess(
                galaxies_filename, 
                survey_name,
                in_directory, 
                out_directory, 
                #mag_cut=mag_cut,
                dist_metric=dist_metric,
                min_z=min_z, 
                max_z=max_z,
                Omega_M=Omega_M,
                #h=h,
            )

            print("Dist limits: ", dist_limits)
            print("Number of galaxies read in: ", len(galaxy_data_table))
            ################################################################################




            ################################################################################
            # GENERATE MASK
            #-------------------------------------------------------------------------------
            mask, mask_resolution = generate_mask(galaxy_data_table, 
                                                  max_z, 
                                                  survey_name,
                                                  out_directory,
                                                  dist_metric=dist_metric, 
                                                  smooth_mask=True,
                                                  #h=h,
                                                  )

            # The mask is automatically saved in the survey_name+'VoidFinder_Output.fits'
            # file and can be read in for future use
            ################################################################################




            ################################################################################
            # FILTER GALAXIES
            #-------------------------------------------------------------------------------
            # If you are rerunning the code, you can comment out the mask generation step 
            # above and just load it here instead. Use something in the vein of the below 
            # (untested) code:
            # with fits.open(out_directory+survey_name+'_VoidFinder_Output.fits') as output:
            #   mask = output['MASK'].data
            #   mask_resolution  = output['MASK'].header['MSKRES']
            #   dist_limits = (output['PRIMARY'].header['DLIML'], output['PRIMARY'].header['DLIMU'])

            wall_coords_xyz, field_coords_xyz = filter_galaxies(galaxy_data_table,
                                                                survey_name,
                                                                out_directory,
                                                                dist_limits=dist_limits,
                                                                magnitude_limit=magnitude_limit,
                                                                #mag_cut_flag=mag_cut,
                                                                #rm_isolated_flag=rm_isolated,
                                                                dist_metric=dist_metric,
                                                                capitalize_colnames=True,
                                                                #h=h,
                                                                )
            print("Number of filtered galaxies: ", len(wall_coords_xyz)+len(field_coords_xyz))
            del galaxy_data_table
            # The galaxies are automatically saved in the survey_name+'VoidFinder_Output.fits'
            # file and can be read in for future use, o long as the write_table parameter of
            # filter_galaxies is True (it is True by default)
            ################################################################################


            coords_min = np.min(np.concatenate([wall_coords_xyz, field_coords_xyz]), axis=0)


            ################################################################################
            # FIND VOIDS
            #-------------------------------------------------------------------------------
            # Again, if you are running the code and have not changed any of the above steps 
            # from a previous run, you can comment out most of the above function calls and 
            # load all the details in here to start over. Use something in the vein of the 
            # below (untested) code:
            # with fits.open(out_directory+survey_name+'_VoidFinder_Output.fits') as output:
            #   wall_coords_xyz = np.array([output['WALL'].data['X'],
            #                               output['WALL'].data['Y'],
            #                               output['WALL'].data['Z']]).T



            find_voids(wall_coords_xyz, 
                       survey_name,
                       out_directory,
                       mask_type='ra_dec_z',
                       mask=mask, 
                       mask_resolution=mask_resolution,
                       dist_limits=dist_limits,
                       grid_origin=coords_min,
                       #save_after=50000,
                       #use_start_checkpoint=True,
                       verbose=1,
                       num_cpus=num_cpus,
                       capitalize_colnames=True)

print('Void-finding complete')            
            
def prepare_catalog (voids_path, galaxies):
    if os.path.exists(voids_path):
        print("exists:", voids_path)
        vfc = vc.VoidFinderCatalog(voids_path)
        vfc.add_galaxies(galaxies) 
        print('r_eff:',voids_path)
        vfc.calculate_r_eff()
        print('vflag_eff:',voids_path)
        vfc.calculate_vflag(overwrite=True) #Restarting the VoidFinding won't clear the previous vflag output
        print('complete:',voids_path)
    else:
        print("failure:", voids_path)

if calc_void_properties:
    print('Calculating Void Properites')
    print('Launching processes')
    processes = []

    for voids_path in sdss_voids:
        processes.append(Process(target = prepare_catalog,args=[voids_path, sdss_galaxies]))

    for voids_path in ngc_voids:
        processes.append(Process(target = prepare_catalog,args=[voids_path, ngc_galaxies]))

    for voids_path in sgc_voids:
        processes.append(Process(target = prepare_catalog,args=[voids_path, sgc_galaxies]))

    for p in processes:
        p.start()

    for p in processes:
        p.join()
