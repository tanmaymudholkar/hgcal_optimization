#!/bin/bash

./plotRechitDistributions.py --inputDigiFileNamePrefix=Digi_Pu200_IC3_version33_model2_Minbias_14TeV_10Events_ --outputFileName=recHitDistributions_zeroNoise --resetOutputDir
./plotRechitDistributions.py --inputDigiFileNamePrefix=DigivariableNoise_Pu200_IC3_version33_model2_Minbias_14TeV_10Events_ --outputFileName=recHitDistributions_variableNoise
./plotRechitDistributions.py --inputDigiFileNamePrefix=Digi_Pu200_IC3_bhnoise0.20mips_version33_model2_Minbias_14TeV_10Events_ --outputFileName=recHitDistributions_constNoise