#!/usr/bin/env python

import os,sys
# import optparse
# import commands
# import math
# import random

# random.seed()

twcu_t = 39.6
tw_t = 41.6
x0wcu = 5.1223
x0w = 3.50418
rwcu_t = twcu_t/x0wcu
rw_t = tw_t/x0w

x0air = 303921
x0pcb = 187.31
x0si = 93.6607
x0cu = 14.3558
x0cfm = 3494.61
lcueven = 1.0
lcfm = 1.0
lair = 3.0
lpcb = 2.0
lsi = 0.3
lcuodd = 6.0;

rw_min = lcueven/x0cu + lcfm/x0cfm + lair/x0air + lpcb/x0pcb + lsi/x0si
rwcu_min = lcuodd/x0cu + lsi/x0si + lpcb/x0pcb + lair/x0air

# print "rw min is ",rw_min
# print "rwcu min is ",rwcu_min

# for index_f in range(0,11):
# rwcu_f = index_f*0.1
rwcu_f = rwcu_t/8
rwcu_m = rwcu_t/8
rw_f = -rw_min + rwcu_min + rwcu_f
twcu_f = 0.5*rwcu_f*x0wcu
tw_f = rw_f*x0w
# for index_m in range(0,11):
# rwcu_m = index_m*0.1
rw_m = -rw_min + rwcu_min + rwcu_m
twcu_m = 0.5*rwcu_m*x0wcu
tw_m = rw_m*x0w
rwcu_b = 0.25*rwcu_t - rwcu_f - rwcu_m
rw_b = 0.25*rw_t - rw_f - rw_m
twcu_b = 0.5*rwcu_b*x0wcu
tw_b = rw_b*x0w

if (tw_f < 0 or tw_m < 0 or tw_b < 0 or twcu_f < 0 or twcu_m < 0 or twcu_b < 0):
    print "something has gone wrong here: you are not submitting the following:"
    os.system("echo \"./submitProdDigi.py -s 2nd -q 1nw -t V05-02-04 -v 45 -m 2 -a 2.1 -b 0 -d gamma -n 1000 -g -o /afs/cern.ch/work/t/tmudholk/public/simulation_results/syst_opt -w %.1f -W %.1f\""%(rwcu_f,rwcu_m))
else:
    # os.system("echo \"you would have submitted the following:\"")
    os.system("echo \"./submitProdDigi.py -s 2nd -q 1nw -t V05-02-04 -v 45 -m 2 -a 2.1 -b 0 -d gamma -n 1000 -g -o /afs/cern.ch/work/t/tmudholk/public/simulation_results/syst_opt -w %.6f -W %.6f\""%(rwcu_f,rwcu_m))
    # os.system("./submitProdDigi.py -s 2nd -q 1nw -t V05-02-04 -v 45 -m 2 -a 2.1 -b 0 -d gamma -n 1000 -g -o /afs/cern.ch/work/t/tmudholk/public/simulation_results/syst_opt -w %.1f -W %.1f"%(rwcu_f,rwcu_m))
