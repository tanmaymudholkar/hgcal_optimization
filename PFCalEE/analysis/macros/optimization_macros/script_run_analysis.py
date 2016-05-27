#!/usr/bin/env python

from __future__ import print_function, division
from os import getcwd, system

version_number = 30
version_name = "version30_test_latest_full"
datadir = "/afs/cern.ch/user/t/tmudholk/private/mount/eos_mount/cms/store/user/tmudholk/hexagonal_geometry/version_30/githexaV02-01-01/gamma"
current_directory = getcwd()
outputdir = current_directory
threshold_adc = 5
digi_or_raw_switch = 3
use_modified_weighting_scheme = "false"
multiple_runs_present = "true"
get_hits_distribution = "true"
get_shower_profile = "true"
get_xy_positions_gen_particle = "true"
get_xy_positions_by_layer = "true"
max_events = 0

root_executable_string = "root -b -q \"%s/run_analysis.C++(%i, \\\"%s\\\", \\\"%s\\\", \\\"%s\\\", %.1f, %i, %s, %s, %s, %s, %s, %s, %i)\""%(current_directory, version_number, version_name, datadir, outputdir, threshold_adc, digi_or_raw_switch, use_modified_weighting_scheme, multiple_runs_present, get_hits_distribution, get_shower_profile, get_xy_positions_gen_particle, get_xy_positions_by_layer, max_events)

print ("Calling: %s"%(root_executable_string))
system(root_executable_string)
