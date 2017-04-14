#!/bin/bash

minR=900
maxR=2700
nRadialBins=16
totalNumberOfCellsDir=totalNumberOfCellsInRadialBinsInfo
occupiedNumberOfCellsDir=occupiedNumberOfCellsInRadialBinsInfo
outputDir=occupancy2DPlots
# outputFileNamePrefix=occupiedCellsTrees_Minbias_14TeV_10Events_radial
threshold=0.5

./compileGetOccupancyMapsFromFiles.sh &&\
    ./getOccupancyMapsFromFiles.out ${minR} ${maxR} ${nRadialBins} ${totalNumberOfCellsDir} ${occupiedNumberOfCellsDir} ${outputDir} ${threshold}
