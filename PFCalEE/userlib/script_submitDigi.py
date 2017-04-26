#!/usr/bin/env python

from __future__ import print_function, division
import os, math

IS_TESTRUN = False
eosSave = False
etlist = [100]
typePrefix = 'pionGun'
datatype = 'pi+'
# hepMCInputFileNamePrefix = "/afs/cern.ch/work/t/tmudholk/public/Minbias_14TeV_HepMC/" + typePrefix
PileupSourcePath = "root://eoscms//eos/cms/store/cmst3/group/hgcal/BHStudies_Standalone/Minbias_14TeV/"
queue = '8nh'
nEvts = 10
# nAvgPU = 200
nAvgPU = 0
# startRunNumber = 1
gittag = 'SingleCell'
version = 33
model = 2
# eta = 1.6
# startIndex = 1
# stopIndex = 110
outputPath = "/afs/cern.ch/work/t/tmudholk/public/bhStudies/" + typePrefix
eosSavePath = "/store/cmst3/group/hgcal/BHStudies_Standalone/" + typePrefix
BHNoisemips = 0
versionName = "noNoise_lowThreshold_withExtraInfo"
forceInputDir = "/afs/cern.ch/work/t/tmudholk/public/bhStudies/pionGun/simulated"

# for index in range(startIndex, 1+stopIndex):
#     # hepMCInputFileName = hepMCInputFileNamePrefix + "_%04d"%(index) + ".dat"
#     suffix = typePrefix + "_%04d"%(index)
#     command = "./submitDigi.py -q {queue} -v {version} -m {model} -P {nAvgPU} -D {PileupSourcePath} -x {suffix} -n {nEvts} -o {outputPath}".format(**locals())
#     if (eosSave): command += " -e {eosSavePath} -E {eosSavePath}".format(**locals())
#     if (BHNoisemips > 0): command += " -B {BHNoisemips}".format(**locals())
#     if (len(versionName) > 0): command += " -N {versionName}".format(**locals())
#     if TESTRUN: command += " -S"
#     print ("About to execute: {command}".format(**locals()))
#     # os.system("./submitDigi.py -S -q {queue} -v {version} -m {model} -P {nAvgPU} -D {PileupSourcePath} -x {suffix} -n {nEvts} -o {outputPath} -e {eosSavePath} -E {eosSavePath} -B {BHNoisemips}".format(**locals()))
#     os.system(command)

# for et in etlist:
#     for runNumber in range(startRunNumber,startRunNumber+nRuns):
#         os.system("./submitDigi.py -S -q {queue} -t {gittag} -r {runNumber} -v {version} -m {model} -a {eta} -T {et} -d {datatype} -n {nEvts} -o {outputDirectory} -g".format(**locals()))

minEta = 1.5
maxEta = 2.1
deltaEta = 0.01
nEtas = 1 + int(math.floor((maxEta - minEta)/deltaEta + 0.5))


for etaCounter in range(0, nEtas):
    eta = minEta + deltaEta*etaCounter
    for et in etlist:
        # suffix = typePrefix
        # command = "./submitDigi.py -q {queue} -t {gittag} -v {version} -m {model} -a {eta} -T {et} -d {datatype} -x {suffix} -n {nEvts} -o {outputPath} -g".format(**locals())
        command = "./submitDigi.py -s {queue} -t {gittag} -v {version} -m {model} -a {eta} -T {et} -d {datatype} -n {nEvts} -o {outputPath} -g".format(**locals())
        if (nAvgPU > 0): command += " -P {nAvgPU} -D {PileupSourcePath}".format(**locals())
        if (eosSave): command += " -e {eosSavePath} -E {eosSavePath}".format(**locals())
        command += " -B {BHNoisemips}".format(**locals())
        if (len(versionName) > 0): command += " -N {versionName}".format(**locals())
        if (len(forceInputDir) > 0): command += " -F {forceInputDir}".format(**locals())
        if (IS_TESTRUN): command += ' -S'
        print ("About to execute: {command}".format(**locals()))
        os.system(command)
