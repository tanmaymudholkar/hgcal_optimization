#!/usr/bin/env python

from __future__ import print_function, division
import os

IS_TESTRUN=False

etlist = [100]
nEvts = 10
# nRuns = 10
# startRunNumber = 1
queue = '1nh'
gittag = 'SingleCell'
version = 33
model = 2
# eta = 1.6
nEtas = 151
minEta = 1.5
maxEta = 3.0
datatype = 'pi+'
outputDirectory = '/afs/cern.ch/work/t/tmudholk/public/bhStudies/pionGun'

deltaEta = (maxEta - minEta)/(nEtas-1)
for etaCounter in range(0, nEtas):
    eta = minEta + deltaEta*etaCounter
    for et in etlist:
        # for runNumber in range(startRunNumber,startRunNumber+nRuns):
            # commandToSubmit = "./submitProd.py -q {queue} -t {gittag} -r {runNumber} -v {version} -m {model} -a {eta} -T {et} -d {datatype} -n {nEvts} -o {outputDirectory} -g".format(**locals())
        commandToSubmit = "./submitProd.py -q {queue} -t {gittag} -v {version} -m {model} -a {eta} -T {et} -d {datatype} -n {nEvts} -o {outputDirectory} -g".format(**locals())
        if (IS_TESTRUN): commandToSubmit += ' -S'
        print ("About to execute: {commandToSubmit}".format(**locals()))
        os.system(commandToSubmit)
