#!/bin/bash

minR=900
maxR=2700
nRadialBins=16
dataDir=root://eoscms//eos/cms/store/cmst3/group/hgcal/BHStudies_Standalone/Minbias_14TeV_10Events/
# inputFileNamePrefix=Digi_Pu200_IC3_version33_model2_Minbias_14TeV_10Events_
# inputFileNamePrefix=Digi_Pu200_IC3_bhnoise0.20mips_version33_model2_Minbias_14TeV_10Events_
inputFileNamePrefix=DigivariableNoise_Pu200_IC3_version33_model2_Minbias_14TeV_10Events_
# outputDir=occupiedNumberOfCellsInRadialBinsInfo_base
# outputDir=occupiedNumberOfCellsInRadialBinsInfo_constantNoise
outputDir=occupiedNumberOfCellsInRadialBinsInfo_variableNoise
thresholdsList=(0.5 5.0)
hardcodedMaxDataCounter=0

./compileCountOccupiedRadialBins.sh &&\
    for threshold in ${thresholdsList[@]}; do
        ./countOccupiedRadialBins.out ${minR} ${maxR} ${nRadialBins} ${dataDir} ${inputFileNamePrefix} ${outputDir} ${threshold} ${hardcodedMaxDataCounter}
    done
