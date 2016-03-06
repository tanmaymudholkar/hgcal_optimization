#!/usr/bin/env python

from __future__ import print_function

import os

for index_f in range(0,11):
    rwcu_f = index_f*0.1
    for index_m in range(0,11):
        rwcu_m = index_m*0.1
        if(not(rwcu_f==1.0 and rwcu_m==1.0)):
            os.system("root -b -q \"/export/home/tmudholk/research/HGCstandalone/analysis/macros/plotE_resolution_high_stats_syst_opt.C++(45,\\\"high_stats_syst_opt\\\",\\\"/export/cmss2/tmudholk/HGCal/high_stats_syst_opt/et100\\\",\\\"/export/home/tmudholk/research/HGCstandalone/analysis/macros\\\",100.0,%.1f,%.1f)\""%(rwcu_f,rwcu_m))
