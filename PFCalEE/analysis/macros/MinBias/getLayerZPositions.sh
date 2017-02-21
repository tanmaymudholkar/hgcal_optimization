#!/bin/bash

# cat /afs/cern.ch/user/t/tmudholk/public/research/hgcal_minbias/PFCal/PFCalEE/LSFJOB_807527118/STDOUT | egrep "(Si|Scintillator) layer"
awk '/(Si|Scintillator) layer/ {print $4, "    ", $8}' /afs/cern.ch/user/t/tmudholk/public/research/hgcal_minbias/PFCal/PFCalEE/LSFJOB_807527118/STDOUT | sed "s|zpos=\([^m]*\)mm|\1|" > getLayerZPositionsTempOutput.dat
./getLayerZPositionsHelper.py > layerZPositions.dat
rm getLayerZPositionsTempOutput.dat
