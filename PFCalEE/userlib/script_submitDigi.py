#!/usr/bin/env python

from __future__ import print_function, division
import os

TESTRUN=False
# etlist = [100]
hepMCTypePrefix = 'Minbias_14TeV_10Events'
# hepMCInputFileNamePrefix = "/afs/cern.ch/work/t/tmudholk/public/Minbias_14TeV_HepMC/" + hepMCTypePrefix
PileupSourcePath = "root://eoscms//eos/cms/store/cmst3/group/hgcal/BHStudies_Standalone/Minbias_14TeV/"
queue = '1nd'
nEvts = 10
nAvgPU = 200
# startRunNumber = 1
# gittag = 'ModifiedGeometry'
version = 33
model = 2
# eta = 1.6
# datatype = 'pi+'
startIndex = 1
stopIndex = 110
outputPath = "/afs/cern.ch/work/t/tmudholk/public/bhStudies/" + hepMCTypePrefix
eosSavePath = "/store/cmst3/group/hgcal/BHStudies_Standalone/" + hepMCTypePrefix
BHNoisemips = 0
versionName = "noNoise_lowThreshold"

for index in range(startIndex, 1+stopIndex):
    # hepMCInputFileName = hepMCInputFileNamePrefix + "_%04d"%(index) + ".dat"
    suffix = hepMCTypePrefix + "_%04d"%(index)
    command = "./submitDigi.py -q {queue} -v {version} -m {model} -P {nAvgPU} -D {PileupSourcePath} -x {suffix} -n {nEvts} -o {outputPath} -e {eosSavePath} -E {eosSavePath}".format(**locals())
    if (BHNoisemips > 0): command += " -B {BHNoisemips}".format(**locals())
    if (len(versionName) > 0): command += " -N {versionName}".format(**locals())
    if TESTRUN: command += " -S"
    print ("About to execute: {command}".format(**locals()))
    # os.system("./submitDigi.py -S -q {queue} -v {version} -m {model} -P {nAvgPU} -D {PileupSourcePath} -x {suffix} -n {nEvts} -o {outputPath} -e {eosSavePath} -E {eosSavePath} -B {BHNoisemips}".format(**locals()))
    os.system(command)

# for et in etlist:
#     for runNumber in range(startRunNumber,startRunNumber+nRuns):
#         os.system("./submitDigi.py -S -q {queue} -t {gittag} -r {runNumber} -v {version} -m {model} -a {eta} -T {et} -d {datatype} -n {nEvts} -o {outputDirectory} -g".format(**locals()))


# for pythiaIndex in range(0,1):
#     os.system("./submitDigi.py -q 8nh -t HGCalMinBias -r %i -v 33 -m 2 -d PythiaTest -n 20 -o /afs/cern.ch/work/t/tmudholk/public/simulation_results/minbias -e /store/cmst3/group/hgcal/HGCalMinbias -E /store/cmst3/group/hgcal/HGCalMinbias -N noXTalk_Test_noNoise_ICoff -S"%(pythiaIndex))
