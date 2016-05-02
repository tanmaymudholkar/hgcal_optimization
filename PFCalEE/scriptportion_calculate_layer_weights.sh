#!/bin/bash

awk '/ layer / {print $4;}' g4copy.log > layer_numbers_for_weight
awk '/ layer / {print $2;}' g4copy.log > layer_types_for_weight
awk '/ layer / {print $5;}' g4copy.log | sed 's:dEdx=\(.*\):\1:' > dEdXs_for_weight
awk '/ layer / {print $6;}' g4copy.log | sed 's:X0=\(.*\):\1:' > X0s_for_weight
awk '/ layer / {print $9;}' g4copy.log | sed 's:w=\(.*\)mm:\1:' > thicknesses_for_weight
awk '/ layer / {print $2"    "$9;}' g4copy.log | sed 's:w=\(.*\)mm:\1:' > geometry_of_detector
