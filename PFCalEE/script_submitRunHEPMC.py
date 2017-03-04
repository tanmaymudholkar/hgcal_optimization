#!/usr/bin/env python

import os

# hepMCInputFileNames = ['/afs/cern.ch/user/t/tmudholk/public/research/pythia/pythia8223/examples/HGCalHepMC_allSoftQCD.txt', '/afs/cern.ch/user/t/tmudholk/public/research/pythia/pythia8223/examples/HGCalHepMC_allHardQCD.txt', '/afs/cern.ch/user/t/tmudholk/public/research/pythia/pythia8223/examples/HGCalHepMC_allSoftAndHardQCD.txt']
hepMCInputFileNames = ['/afs/cern.ch/user/t/tmudholk/public/research/standalone/CMSSW_8_0_26_patch1/src/VBFHtoBB_HepMC.dat']
# hepMCTypes = ['softQCD', 'hardQCD', 'softAndHardQCD']
hepMCTypes = ['VBFHtoBB']
hepMCNEvtsList = [5]
versionNumber = 33
modelNumber = 2
suffix = "modifiedGeometry"

# for index in range(0,len(hep)):
#     os.system("./submitRunHEPMC.py -q 1nw -v 33 -m 2 -f /afs/cern.ch/work/t/tmudholk/public/pythia/Pythia140305_000%03d.dat -s 000%03d -n 5000 -o /afs/cern.ch/work/t/tmudholk/public/simulation_results/minbias -e /store/cmst3/group/hgcal/HGCalMinbias"%(pythiaIndex, pythiaIndex))

for index in range(0, len(hepMCInputFileNames)):
    hepMCInputFileName = hepMCInputFileNames[index]
    hepMCType = hepMCTypes[index]
    hepMCNEvts = hepMCNEvtsList[index]
    os.system("./submitRunHEPMC.py -S -q 1nd -v " + str(versionNumber) + " -m " + str(modelNumber) + " -f " + hepMCInputFileName + " -s " + suffix + " -n " + str(hepMCNEvts) + " -o /afs/cern.ch/work/t/tmudholk/public/bhStudies/HepMCTest_" + hepMCType + " -e /store/cmst3/group/hgcal/BHStudies_Standalone/HepMCTest_" + hepMCType)
