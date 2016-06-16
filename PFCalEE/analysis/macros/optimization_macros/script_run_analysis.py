#!/usr/bin/env python

from __future__ import print_function, division
from os import getcwd, system

version_number = 30
version_name = "version30_suppression_study"
datadir = "/afs/cern.ch/user/t/tmudholk/private/mount/eos_mount/cms/store/user/tmudholk/hexagonal_geometry/version_30/githexaV02-01-01/gamma"
current_directory = getcwd()
outputdir = current_directory
# thresholds_adc = [5, 10, 20, 35, 50]
thresholds_adc = [75, 100]
digi_or_raw_switch = 3
use_modified_weighting_scheme = "false"
multiple_runs_present = "true"
get_hits_distribution = "true"
get_shower_profile = "true"
get_xy_positions_gen_particle = "true"
get_xy_positions_by_layer = "true"
max_events = 0

system ("mv ~/.bash_logout ~/.bash_logout_temp")

system("if mountpoint -q \"/afs/cern.ch/user/t/tmudholk/private/mount/eos_mount\" ; then echo \"already mounted, no need to mount again\"; else /afs/cern.ch/project/eos/installation/0.3.84-aquamarine/bin/eos.select -b fuse mount /afs/cern.ch/user/t/tmudholk/private/mount/eos_mount; fi")

for threshold_adc in thresholds_adc:
    root_executable_string = "root -b -q \"%s/run_analysis.C++(%i, \\\"%s\\\", \\\"%s\\\", \\\"%s\\\", %.1f, %i, %s, %s, %s, %s, %s, %s, %i)\""%(current_directory, version_number, version_name, datadir, outputdir, threshold_adc, digi_or_raw_switch, use_modified_weighting_scheme, multiple_runs_present, get_hits_distribution, get_shower_profile, get_xy_positions_gen_particle, get_xy_positions_by_layer, max_events)
    
    print ("Calling: %s"%(root_executable_string))
    system(root_executable_string)

system("mv ~/.bash_logout_temp ~/.bash_logout")
system("source ~/.bash_logout")
system("cd %s"%(current_directory))    
system("echo \"Analysis completed for version_number = %i, version_name = %s, datadir = %s, outputdir = %s, thresholds_adc = %s, digi_or_raw_switch = %i, use_modified_weighting_scheme = %s, multiple_runs_present = %s, get_hits_distribution = %s, get_shower_profile = %s, get_xy_positions_gen_particle = %s, get_xy_positions_by_layer = %s, max_events = %i\" > run_analysis_completed.txt"%(version_number, version_name, datadir, outputdir, thresholds_adc, digi_or_raw_switch, use_modified_weighting_scheme, multiple_runs_present, get_hits_distribution, get_shower_profile, get_xy_positions_gen_particle, get_xy_positions_by_layer, max_events))
