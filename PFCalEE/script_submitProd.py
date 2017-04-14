#!/usr/bin/env python

from __future__ import print_function, division
import os

etlist = [100]
nEvts = 10
nRuns = 5
startRunNumber = 1
queue = '8nh'
gittag = 'ModifiedGeometry'
version = 33
model = 2
eta = 1.6
datatype = 'pi+'
outputDirectory = '/afs/cern.ch/work/t/tmudholk/public/bhStudies/pionGun'

for et in etlist:
    for runNumber in range(startRunNumber,startRunNumber+nRuns):
        os.system("./submitProd.py -S -q {queue} -t {gittag} -r {runNumber} -v {version} -m {model} -a {eta} -T {et} -d {datatype} -n {nEvts} -o {outputDirectory} -g".format(**locals()))
