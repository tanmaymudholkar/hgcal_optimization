#!/bin/bash

minR=900
maxR=2700
nRadialBins=16
dataDir=root://eoscms//eos/cms/store/cmst3/group/hgcal/BHStudies_Standalone/Minbias_14TeV_10Events/
# inputFileNamePrefix=Digi_Pu200_IC3_version33_model2_Minbias_14TeV_10Events_
inputFileNamePrefix=Digi_Pu200_IC3_bhnoise0.20mips_version33_model2_Minbias_14TeV_10Events_
# inputRadialRangesFileName=radialBinBoundaries.dat
# totalNumberOfCellsInRadialBinsDir=totalNumberOfCellsInRadialBinsInfo
outputDir=occupiedNumberOfCellsInRadialBinsInfo_withnoise
# outputFileNamePrefix=occupiedCellsTrees_Minbias_14TeV_10Events_radial
threshold=0.5
hardcodedMaxDataCounter=0

./compileCompareOccupancyProfiles.sh &&\
    ./compareOccupancyProfiles.out ${minR} ${maxR} ${nRadialBins} ${dataDir} ${inputFileNamePrefix} ${outputDir} ${threshold} ${hardcodedMaxDataCounter}
