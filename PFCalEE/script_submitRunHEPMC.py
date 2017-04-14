#!/usr/bin/env python

import os

# hepMCInputFileNames = ['/afs/cern.ch/user/t/tmudholk/public/research/pythia/pythia8223/examples/HGCalHepMC_allSoftQCD.txt', '/afs/cern.ch/user/t/tmudholk/public/research/pythia/pythia8223/examples/HGCalHepMC_allHardQCD.txt', '/afs/cern.ch/user/t/tmudholk/public/research/pythia/pythia8223/examples/HGCalHepMC_allSoftAndHardQCD.txt']
# hepMCInputFileNames = ['/afs/cern.ch/user/t/tmudholk/public/research/hgcal_BHStudies/CMSSW_8_0_26_patch1/src/GenEvent_ASCII.dat']
# hepMCInputFileNames = ['/afs/cern.ch/user/t/tmudholk/public/research/hgcal_BHStudies/mcConversion_900pre6/CMSSW_9_0_0_pre6/src/Minbias_14TeV_9950Events_test.dat']
queue = "1nh"
hepMCTypePrefix = 'Minbias_14TeV_10Events'
hepMCInputFileNamePrefix = "/afs/cern.ch/work/t/tmudholk/public/Minbias_14TeV_HepMC/" + hepMCTypePrefix + "/" + hepMCTypePrefix
# hepMCNEvtsList = [20, 20, 20]
# hepMCNEvtsList = [9950]
nEvtsPerFile = 10
version = 33
model = 2
# suffix = "modifiedGeometry"
startIndex = 1
stopIndex = 110
outputPath = "/afs/cern.ch/work/t/tmudholk/public/bhStudies/" + hepMCTypePrefix
eosSavePath = "/store/cmst3/group/hgcal/BHStudies_Standalone/" + hepMCTypePrefix

# for pythiaIndex in range(123,124):

for index in range(startIndex, 1+stopIndex):
    hepMCInputFileName = hepMCInputFileNamePrefix + "_%04d"%(index) + ".dat"
    suffix = hepMCTypePrefix + "_%04d"%(index)
    # hepMCType = hepMCTypes[index]
    # hepMCNEvts = hepMCNEvtsList[index]
    # os.system("./submitRunHEPMC.py -S -q 8nh -v " + str(versionNumber) + " -m " + str(modelNumber) + " -f " + hepMCInputFileName + " -s " + suffix + " -n " + str(hepMCNEvts) + " -o /afs/cern.ch/work/t/tmudholk/public/bhStudies/HepMCTest_" + hepMCType + " -e /store/cmst3/group/hgcal/BHStudies_Standalone/HepMCTest_" + hepMCType)
    # os.system("./submitRunHEPMC.py -q " + queue + " -v " + str(versionNumber) + " -m " + str(modelNumber) + " -f " + hepMCInputFileName + " -x " + suffix + " -n " + str(nEvtsPerFile) + " -o " + outputPath + " -e " + eosSavePath)
    os.system("./submitRunHEPMC.py -q {queue} -v {version} -m {model} -f {hepMCInputFileName} -x {suffix} -n {nEvtsPerFile} -o {outputPath} -e {eosSavePath}".format(**locals()))
