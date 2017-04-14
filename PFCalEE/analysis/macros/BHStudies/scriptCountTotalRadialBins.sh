#!/bin/bash

minR=900
maxR=2700
nRadialBins=16
inputFileName=root://eoscms//eos/cms/store/cmst3/group/hgcal/BHStudies_Standalone/Minbias_14TeV_10Events/Digi_Pu200_IC3_version33_model2_Minbias_14TeV_10Events_0001.root
outputDir=totalNumberOfCellsInRadialBinsInfo

./compileCountTotalRadialBins.sh &&\
    ./countTotalRadialBins.out ${minR} ${maxR} ${nRadialBins} ${inputFileName} ${outputDir}
