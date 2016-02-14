#!/usr/bin/env python

import os,sys
import optparse
import commands
import math
import numpy as np
#import matplotlib.pyplot as plt

# usage = 'usage: %prog [options]'
# parser = optparse.OptionParser(usage)
# parser.add_option('-p', '--path-of-g4', dest='path', help='path of g4.log', default='')
# (opt, args) = parser.parse_args()

twcu_t = 39.6
tw_t = 41.6
x0wcu = 5.1223
x0w = 3.50418
rwcu_t = twcu_t/x0wcu
rw_t = tw_t/x0w

x0air = 303921
x0pcb = 187.31
x0si = 93.6607
x0cu = 14.3558
x0cfm = 3494.61
lcueven = 1.0
lcfm = 1.0
lair = 3.0
lpcb = 2.0
lsi = 0.3
lcuodd = 6.0;

rw_min = lcueven/x0cu + lcfm/x0cfm + lair/x0air + lpcb/x0pcb + lsi/x0si
rwcu_min = lcuodd/x0cu + lsi/x0si + lpcb/x0pcb + lair/x0air

for index_f in range(0,11):
    rwcu_f = index_f*0.1
    rw_f = -rw_min + rwcu_min + rwcu_f
    twcu_f = 0.5*rwcu_f*x0wcu
    tw_f = rw_f*x0w
    for index_m in range(0,11):
        rwcu_m = index_m*0.1
        rw_m = -rw_min + rwcu_min + rwcu_m
        twcu_m = 0.5*rwcu_m*x0wcu
        tw_m = rw_m*x0w
        rwcu_b = 0.25*rwcu_t - rwcu_f - rwcu_m
        rw_b = 0.25*rw_t - rw_f - rw_m
        twcu_b = 0.5*rwcu_b*x0wcu
        tw_b = rw_b*x0w
        if (tw_f < 0 or tw_m < 0 or tw_b < 0 or twcu_f < 0 or twcu_m < 0 or twcu_b < 0):
            print "this g4 does not exist"
            continue
        else:
            path_to_g4='workfolder_results/syst_opt/git_V05-02-04/version_45/model_2/rwcuf%.1f/rwcum%.1f/gamma/BOFF/et_5/eta_2.100'%(rwcu_f,rwcu_m)
            os.system('ls %s/g4.log'%(path_to_g4))
            os.system('cp %s/g4.log g4copy.log'%(path_to_g4))
            os.system('./scriptportion_calculate_layer_weights.sh')
            layer_numbers_raw=np.loadtxt('layer_numbers_for_weight')
            layer_numbers=np.array(layer_numbers_raw)
            layer_types_raw=np.genfromtxt('layer_types_for_weight',dtype=None)
            layer_types=np.array(layer_types_raw,dtype=None)
            dEdXs_raw=np.loadtxt('dEdXs_for_weight')
            dEdXs=np.array(dEdXs_raw)
            thicknesses_raw=np.loadtxt('thicknesses_for_weight')
            thicknesses=np.array(thicknesses_raw)
            # print layer_numbers,layer_types,dEdXs,thicknesses
            # print layer_numbers.shape[0],layer_types.shape[0],dEdXs.shape[0],thicknesses.shape[0]
            if((layer_types.shape[0]-layer_numbers.shape[0])*(layer_types.shape[0]-layer_numbers.shape[0])+(dEdXs.shape[0]-layer_numbers.shape[0])*(dEdXs.shape[0]-layer_numbers.shape[0])+(thicknesses.shape[0]-layer_numbers.shape[0])*(thicknesses.shape[0]-layer_numbers.shape[0]) > 0):
                print "something went wrong; not all arrays are of the same length"
                sys.exit()
            total_dEdXs_times_thickness=[0]*int(1+np.amax(layer_numbers))
            # print total_dEdXs_times_thickness
            for iterator_over_sublayers in range(layer_numbers.shape[0]):
                total_dEdXs_times_thickness[int(layer_numbers[iterator_over_sublayers])] += dEdXs[iterator_over_sublayers]*thicknesses[iterator_over_sublayers]

            # print total_dEdXs_times_thickness
            normalizer=total_dEdXs_times_thickness[0]
            for iterator_over_layers in range(len(total_dEdXs_times_thickness)):
                total_dEdXs_times_thickness[iterator_over_layers]=total_dEdXs_times_thickness[iterator_over_layers]/normalizer
            
            # print total_dEdXs_times_thickness

            np.savetxt('layer_weights_rwcuf%.1f_rwcum%.1f.dat'%(rwcu_f,rwcu_m),total_dEdXs_times_thickness,fmt='%.3f')
