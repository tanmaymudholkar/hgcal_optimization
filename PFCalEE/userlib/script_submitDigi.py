#!/usr/bin/env python

from __future__ import print_function

import os,sys,time
# import optparse
# import commands
# import math
# import random

# random.seed()

thresholds_adc = [5, 20, 50]
version = 30

for run_number in range(1,103):
    print ("Run number: %i"%(run_number))
    for threshold_adc in thresholds_adc:
        time.sleep(1.1)
        os.system("./submitDigi.py -s 8nh -q 1nd -t hexaV02-01-01 -v %i -m 2 -a 2.1 -b 0 -d gamma -n 100 -T %i -g -o /afs/cern.ch/work/t/tmudholk/public/simulation_results/hexagonal_geometry/version_30 -r %i -e /store/user/tmudholk/hexagonal_geometry/version_%i -E /store/user/tmudholk/hexagonal_geometry/version_%i"%(version, threshold_adc, run_number, version, version))
        # os.system("./submitDigi.py -s 8nh -q 8nh -t hexaV02-01-01 -v 33 -m 2 -a 2.1 -b 0 -d gamma -n 100 -g -o /afs/cern.ch/work/t/tmudholk/public/simulation_results/hexagonal_geometry/version_33 -r %i"%(run_number))
