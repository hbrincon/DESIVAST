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

import numpy as np

from vast.voidfinder import find_voids, wall_field_separation

from vast.voidfinder.preprocessing import load_data_to_Table

import numpy as np
################################################################################

# load mock
# (2000 Mpc/h)^3 box from the fiducial AbacusSummit_base_c000_ph000 sim
#columns of interest include
#rabsmag
#gr
#halo mass
#central flag
#velocity info
#halo ID
"""
t1 = Table.read('/global/cfs/cdirs/desi/cosmosim/SecondGenMocks/AbacusSummit/CubicBox/BGS/v0.1/z0.200/AbacusSummit_base_c000_ph000/BGS_box_ph000.fits')
sel = (t1['x'] < 500) * (t1['y'] < 500) * (t1['z'] < 500)

t1 = t1[sel]
t1.write('bgs_cubic.fits')
"""

################################################################################
# USER INPUTS
#-------------------------------------------------------------------------------

# Coordinate limits of the simulation
xyz_limits = np.array([[0.,0.,0.],[500.,500.,500.]])

# Size of a single grid cell
hole_grid_edge_length = 5.0

# Number of CPUs available for analysis.
# A value of None will use one less than all available CPUs.
num_cpus = None

#-------------------------------------------------------------------------------
# File name details
#-------------------------------------------------------------------------------
# File header

survey_name = 'bgs_cubic'

# Change these directory paths to where your data is stored, and where you want 
# the output to be saved.
in_directory = './'
out_directory = './'


# Input file name
# File format: RA, dec, redshift, comoving distance, absolute magnitude
galaxies_filename = 'bgs_cubic.fits'
#-------------------------------------------------------------------------------



#-------------------------------------------------------------------------------
# Galaxy pruning details
#-------------------------------------------------------------------------------
# Uncomment if you do NOT want to remove galaxies with Mr > -20
# Need to also uncomment relevent input in function calls below
#mag_cut = False
magnitude_limit = -20.0


# Uncomment if you do NOT want to remove isolated galaxies
# Need to also uncomment relevent input in function calls below
#rm_isolated = False
#-------------------------------------------------------------------------------
################################################################################



################################################################################
# Read in data
#-------------------------------------------------------------------------------
# Read in the simulated data
coords_xyz = load_data_to_Table(galaxies_filename)

#magnitude cut
coords_xyz['R_MAG_ABS'].name = 'rabsmag'
coords_xyz = coords_xyz[coords_xyz['rabsmag'] <= magnitude_limit]

#-------------------------------------------------------------------------------
# Restructure the data for the find_voids function
#-------------------------------------------------------------------------------
x = coords_xyz['x']
y = coords_xyz['y']
z = coords_xyz['z']

num_gal = x.shape[0]

coords_xyz = np.concatenate((x.reshape(num_gal,1),
                             y.reshape(num_gal,1),
                             z.reshape(num_gal,1)), axis=1)

#-------------------------------------------------------------------------------
# Remove isolated galaxies
#-------------------------------------------------------------------------------
wall_coords_xyz, field_coords_xyz = wall_field_separation(coords_xyz,
                                                          survey_name=survey_name, 
                                                          out_directory=out_directory,
                                                          write_galaxies=True)




################################################################################
# Find voids
#-------------------------------------------------------------------------------
find_voids(wall_coords_xyz,
           survey_name,
           out_directory,
           mask_type='xyz',
           xyz_limits=xyz_limits,
           #save_after=50000,
           #use_start_checkpoint=True,
           num_cpus=num_cpus,
           batch_size=10000)
