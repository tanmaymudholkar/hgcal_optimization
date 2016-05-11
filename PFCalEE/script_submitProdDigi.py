#!/usr/bin/env python

from __future__ import print_function, division

import os,sys,time
# import optparse
# import commands
# import math
# import random

# random.seed()

version = 30
etalist = [2.1]
etlist=[3,5,10,30,50,70,100,150]

print ("WARNING: only integer energies supported!!!")

for run_number in range(103,106):
    time.sleep(1.1)
    print ("run number: %i"%(run_number))
    for eta in etalist:
        print ("eta: %.1f"%(eta))
        for et in etlist:
            print ("et: %i"%(et))
            os.system("./submitProdDigi.py -s 8nh -q 1nd -t hexaV02-01-01 -r %i -v %i -m 2 -E %i -a %.3f -b 0 -d gamma -n 100 -T 5 -o /afs/cern.ch/work/t/tmudholk/public/simulation_results/hexagonal_geometry/version_%i -e /store/user/tmudholk/hexagonal_geometry/version_%i -g"%(run_number, version, et, eta, version, version))
    # os.system("./submitProd.py -s 8nh -q 1nd -t hexaV02-01-01 -v 33 -m 2 -a 2.1 -b 0 -d gamma -n 100 -g -o /afs/cern.ch/work/t/tmudholk/public/simulation_results/hexagonal_geometry/version_33 -r %i"%(run_number))
