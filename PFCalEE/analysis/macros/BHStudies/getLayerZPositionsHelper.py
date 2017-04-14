#!/usr/bin/env python

from __future__ import print_function, division

inputFileName = "getLayerZPositionsTempOutput.dat"

inputFile = open(inputFileName)
layerNumbers = []
numberOfActiveLayers = []
zPositions = []

runningLayerNumber = -1
for line in inputFile:
    splitline = line.rstrip().split()
    layerNumber = int(splitline[0])
    zPosition = float(splitline[1])
    if (runningLayerNumber != layerNumber):
        layerNumbers.append(layerNumber)
        numberOfActiveLayers.append(1)
        zPositions.append(zPosition)
        runningLayerNumber = layerNumber
    else:
        numberOfActiveLayers[layerNumber] += 1
        zPositions[layerNumber] += zPosition

for layerNumber in layerNumbers:
    avgZPosition = zPositions[layerNumber]/numberOfActiveLayers[layerNumber]
    print("%.2f"%(avgZPosition))
