#!/usr/bin/env python

from __future__ import print_function, division

import os,sys,time
# import optparse
# import commands
# import math
# import random

# random.seed()

version = 100
eta = 2.0
# etlist=[3,5,7,10,20,30,40,50,60,70,80,90,100,125,150,175]

print ("WARNING: only integer energies supported!!!")

for run_number in range(1,511):
    time.sleep(1.01)
    print ("run number: %i"%(run_number))
    os.system("./submitProd.py -q 2nd -t V06e-04-06 -r %i -v %i -m 4 -a %.3f -b 0 -d gamma -n 200 -o /afs/cern.ch/work/t/tmudholk/public/simulation_results/cracks_study -e /store/cmst3/group/hgcal/HGCalCracks -g"%(run_number, version, eta))
