#!/usr/bin/env python

import os

# for pythiaIndex in range(0,1):
#     os.system("./submitDigi.py -q 8nh -t HGCalMinBias -r %i -v 33 -m 2 -d PythiaTest -n 20 -o /afs/cern.ch/work/t/tmudholk/public/simulation_results/minbias -e /store/cmst3/group/hgcal/HGCalMinbias -E /store/cmst3/group/hgcal/HGCalMinbias -N noXTalk_Test_noNoise_ICoff -S"%(pythiaIndex))

# hepMCInputFileNames = ['/afs/cern.ch/user/t/tmudholk/public/research/pythia/pythia8223/examples/HGCalHepMC_gg2bbbarTest.txt', '/afs/cern.ch/user/t/tmudholk/public/research/pythia/pythia8223/examples/HGCalHepMC_qqbar2bbbarTest.txt']
hepMCTypes = ['VBFHtoBB']
hepMCNEvtsList = [5]
versionNumber = 33
modelNumber = 2
versionName = "noiseOff"
inputSuffix = "modifiedGeometry"

for index in range(0, len(hepMCTypes)):
    # hepMCInputFileName = hepMCInputFileNames[index]
    hepMCType = hepMCTypes[index]
    hepMCNEvts = hepMCNEvtsList[index]
    os.system("./submitDigi.py -S -s 8nh -v " + str(versionNumber) + " -m " + str(modelNumber) + " -d  PythiaTest -n " + str(hepMCNEvts) + " -o /afs/cern.ch/work/t/tmudholk/public/bhStudies/HepMCTest_" + hepMCType + " -e /store/cmst3/group/hgcal/BHStudies_Standalone/HepMCTest_" + hepMCType + " -E /store/cmst3/group/hgcal/BHStudies_Standalone/HepMCTest_" + hepMCType + " -N " + versionName + " -x " + inputSuffix)
