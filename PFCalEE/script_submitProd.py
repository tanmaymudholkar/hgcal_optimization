#!/usr/bin/env python

import os,sys,time
# import optparse
# import commands
# import math
# import random

# random.seed()

for run_number in range(12,111):
    time.sleep(1)
    os.system("./submitProd.py -s 8nh -q 1nd -t hexaV02-01-01 -v 33 -m 2 -a 2.1 -b 0 -d gamma -n 100 -g -o /afs/cern.ch/work/t/tmudholk/public/simulation_results/hexagonal_geometry/version_33 -r %i -e /store/user/tmudholk/hexagonal_geometry/version_33"%(run_number))
    # os.system("./submitProd.py -s 8nh -q 1nd -t hexaV02-01-01 -v 33 -m 2 -a 2.1 -b 0 -d gamma -n 100 -g -o /afs/cern.ch/work/t/tmudholk/public/simulation_results/hexagonal_geometry/version_33 -r %i"%(run_number))
