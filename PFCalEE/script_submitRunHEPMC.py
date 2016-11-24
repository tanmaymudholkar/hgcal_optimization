#!/usr/bin/env python

import os

for pythiaIndex in range(22,123):
    os.system("./submitRunHEPMC.py -q 1nw -v 33 -m 2 -f /afs/cern.ch/work/t/tmudholk/public/pythia/Pythia140305_000%03d.dat -s 000%03d -n 5000 -o /afs/cern.ch/work/t/tmudholk/public/simulation_results/minbias -e /store/cmst3/group/hgcal/HGCalMinbias"%(pythiaIndex, pythiaIndex))
