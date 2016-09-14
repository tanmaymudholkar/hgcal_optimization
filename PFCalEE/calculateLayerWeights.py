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
parser.add_option('-n', '--version-number', dest='version_number', help='version number', default=0,type=int)
(opt, args) = parser.parse_args()

os.system('cp %s/g4.log g4copy.log'%(opt.path))
os.system('./calculateLayerWeightsHelper.sh')

layer_numbers_raw=np.loadtxt('layer_numbers_for_weight')
layer_numbers=np.array(layer_numbers_raw)
layer_types_raw=np.genfromtxt('layer_types_for_weight',dtype=None)
layer_types=np.array(layer_types_raw,dtype=None)
dEdXs_raw=np.loadtxt('dEdXs_for_weight')
dEdXs=np.array(dEdXs_raw)
X0s_raw=np.loadtxt('X0s_for_weight')
X0s=np.array(X0s_raw)
thicknesses_raw=np.loadtxt('thicknesses_for_weight')
thicknesses=np.array(thicknesses_raw)

#print layer_numbers,layer_types,dEdXs,thicknesses
#print layer_numbers.shape[0],layer_types.shape[0],dEdXs.shape[0],thicknesses.shape[0]

if((layer_types.shape[0]-layer_numbers.shape[0])*(layer_types.shape[0]-layer_numbers.shape[0])+(dEdXs.shape[0]-layer_numbers.shape[0])*(dEdXs.shape[0]-layer_numbers.shape[0])+(thicknesses.shape[0]-layer_numbers.shape[0])*(thicknesses.shape[0]-layer_numbers.shape[0])+(X0s.shape[0]-layer_numbers.shape[0])*(X0s.shape[0]-layer_numbers.shape[0]) > 0):
    print "something went wrong; not all arrays are of the same length"
    sys.exit()


print "Calculating dE/dx-based weights..."
print
total_dEdXs_times_thickness=[0]*int(1+np.amax(layer_numbers))
print "Initializing... total_dEdXs_times_thickness:"
print total_dEdXs_times_thickness
print
for iterator_over_sublayers in range(layer_numbers.shape[0]):
    total_dEdXs_times_thickness[int(layer_numbers[iterator_over_sublayers])] = total_dEdXs_times_thickness[int(layer_numbers[iterator_over_sublayers])]+dEdXs[iterator_over_sublayers]*thicknesses[iterator_over_sublayers]

print "Calculating unnormalized sums... total_dEdXs_times_thickness:"
print total_dEdXs_times_thickness
print

normalizer=sum(total_dEdXs_times_thickness)/len(total_dEdXs_times_thickness)

print "Normalizing... total_dEdXs_times_thickness:"
for iterator_over_layers in range(len(total_dEdXs_times_thickness)):
    total_dEdXs_times_thickness[iterator_over_layers]=total_dEdXs_times_thickness[iterator_over_layers]/normalizer
print total_dEdXs_times_thickness
print

np.savetxt('layer_weights_version%i.dat'%(opt.version_number),total_dEdXs_times_thickness,fmt='%.3f')

print "Calculating averaged dE/dx-based weights..."
print
total_dEdXs_times_thickness_averaged=[0]*len(total_dEdXs_times_thickness)
print total_dEdXs_times_thickness_averaged

print "Averaging... dEdXs_times_thickness:"
print
for iterator_over_layers in range(len(total_dEdXs_times_thickness)-1):
    total_dEdXs_times_thickness_averaged[iterator_over_layers]=0.5*(total_dEdXs_times_thickness[iterator_over_layers]+total_dEdXs_times_thickness[iterator_over_layers+1])
total_dEdXs_times_thickness_averaged[-1]=total_dEdXs_times_thickness[-1]
print total_dEdXs_times_thickness_averaged

print "Normalizing... total_dEdXs_times_thickness:"
normalizer=sum(total_dEdXs_times_thickness_averaged)/len(total_dEdXs_times_thickness_averaged)
for iterator_over_layers in range(len(total_dEdXs_times_thickness)):
    total_dEdXs_times_thickness_averaged[iterator_over_layers]=total_dEdXs_times_thickness_averaged[iterator_over_layers]/normalizer
print total_dEdXs_times_thickness_averaged

np.savetxt('layer_weights_version%i_averaged.dat'%(opt.version_number),total_dEdXs_times_thickness_averaged,fmt='%.3f')



print "Calculating X0-based weights..."
print
total_inv_X0s_times_thickness=[0]*int(1+np.amax(layer_numbers))
print "Initializing... total_inv_X0_times_thickness:"
print total_inv_X0s_times_thickness
print
for iterator_over_sublayers in range(layer_numbers.shape[0]):
    total_inv_X0s_times_thickness[int(layer_numbers[iterator_over_sublayers])] = total_inv_X0s_times_thickness[int(layer_numbers[iterator_over_sublayers])]+(thicknesses[iterator_over_sublayers]/X0s[iterator_over_sublayers])

print "Calculating unnormalized sums... total_inv_X0s_times_thickness:"
print total_inv_X0s_times_thickness
print

normalizer=sum(total_inv_X0s_times_thickness)/len(total_inv_X0s_times_thickness)
print "Normalizing... total_inv_X0s_times_thickness:"
for iterator_over_layers in range(len(total_inv_X0s_times_thickness)):
    total_inv_X0s_times_thickness[iterator_over_layers]=total_inv_X0s_times_thickness[iterator_over_layers]/normalizer
print total_inv_X0s_times_thickness
print

np.savetxt('layer_weights_x0_version%i.dat'%(opt.version_number),total_inv_X0s_times_thickness,fmt='%.3f')

print "Calculating averaged X0-based weights..."
print
total_inv_X0s_times_thickness_averaged=[0]*len(total_inv_X0s_times_thickness)
print total_inv_X0s_times_thickness_averaged

print "Averaging... inv_X0s_times_thickness:"
print
for iterator_over_layers in range(len(total_inv_X0s_times_thickness)-1):
    total_inv_X0s_times_thickness_averaged[iterator_over_layers]=0.5*(total_inv_X0s_times_thickness[iterator_over_layers]+total_inv_X0s_times_thickness[iterator_over_layers+1])
total_inv_X0s_times_thickness_averaged[-1]=total_inv_X0s_times_thickness[-1]
print total_inv_X0s_times_thickness_averaged

print "Normalizing... total_inv_X0s_times_thickness:"
normalizer=sum(total_inv_X0s_times_thickness_averaged)/len(total_inv_X0s_times_thickness_averaged)
for iterator_over_layers in range(len(total_inv_X0s_times_thickness)):
    total_inv_X0s_times_thickness_averaged[iterator_over_layers]=total_inv_X0s_times_thickness_averaged[iterator_over_layers]/normalizer
print total_inv_X0s_times_thickness_averaged

np.savetxt('layer_weights_x0_version%i_averaged.dat'%(opt.version_number),total_inv_X0s_times_thickness_averaged,fmt='%.3f')
