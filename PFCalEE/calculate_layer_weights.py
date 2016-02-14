#!/usr/bin/env python

import os,sys
import optparse
import commands
import math
import numpy as np
#import matplotlib.pyplot as plt

usage = 'usage: %prog [options]'
parser = optparse.OptionParser(usage)
parser.add_option('-p', '--path-of-g4', dest='path', help='path of g4.log', default='')
(opt, args) = parser.parse_args()

os.system('cp %s/g4.log g4copy.log'%(opt.path))
os.system('./scriptportion_calculate_layer_weights.sh')

layer_numbers_raw=np.loadtxt('layer_numbers_for_weight')
layer_numbers=np.array(layer_numbers_raw)
layer_types_raw=np.genfromtxt('layer_types_for_weight',dtype=None)
layer_types=np.array(layer_types_raw,dtype=None)
dEdXs_raw=np.loadtxt('dEdXs_for_weight')
dEdXs=np.array(dEdXs_raw)
thicknesses_raw=np.loadtxt('thicknesses_for_weight')
thicknesses=np.array(thicknesses_raw)

#print layer_numbers,layer_types,dEdXs,thicknesses
#print layer_numbers.shape[0],layer_types.shape[0],dEdXs.shape[0],thicknesses.shape[0]

if((layer_types.shape[0]-layer_numbers.shape[0])*(layer_types.shape[0]-layer_numbers.shape[0])+(dEdXs.shape[0]-layer_numbers.shape[0])*(dEdXs.shape[0]-layer_numbers.shape[0])+(thicknesses.shape[0]-layer_numbers.shape[0])*(thicknesses.shape[0]-layer_numbers.shape[0]) > 0):
    print "something went wrong; not all arrays are of the same length"
    sys.exit()

total_dEdXs_times_thickness=[0]*int(1+np.amax(layer_numbers))
print total_dEdXs_times_thickness
for iterator_over_sublayers in range(layer_numbers.shape[0]):
    total_dEdXs_times_thickness[int(layer_numbers[iterator_over_sublayers])] = total_dEdXs_times_thickness[int(layer_numbers[iterator_over_sublayers])]+dEdXs[iterator_over_sublayers]*thicknesses[iterator_over_sublayers]

print total_dEdXs_times_thickness

normalizer=total_dEdXs_times_thickness[0]

for iterator_over_layers in range(len(total_dEdXs_times_thickness)):
    total_dEdXs_times_thickness[iterator_over_layers]=total_dEdXs_times_thickness[iterator_over_layers]/normalizer

print total_dEdXs_times_thickness

np.savetxt('layer_weights.dat',total_dEdXs_times_thickness,fmt='%.3f')
