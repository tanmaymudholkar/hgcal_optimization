#!/usr/bin/env python

import os

for pythiaIndex in range(0,1):
    os.system("./submitDigi.py -q 8nh -t HGCalMinBias -r %i -v 33 -m 2 -d PythiaTest -n 20 -o /afs/cern.ch/work/t/tmudholk/public/simulation_results/minbias -e /store/cmst3/group/hgcal/HGCalMinbias -E /store/cmst3/group/hgcal/HGCalMinbias -N noXTalk_Test_noNoise_ICoff -S"%(pythiaIndex))
