#!/usr/bin/env python

from __future__ import print_function

import os

for index_f in range(0,11):
    rwcu_f = index_f*0.1
    for index_m in range(0,11):
        rwcu_m = index_m*0.1
        if(not(rwcu_f==1.0 and rwcu_m==1.0)):
            # if (rwcu_f == 0.1 and rwcu_m == 0.2) :
                # os.system("root -b -q \"/home/tmudholk/research/hgcal_analysis/hgcal_optimization/analysis/macros/plotE_resolution_high_stats_syst_opt.C++(45,\\\"high_stats_syst_opt\\\",\\\"/home/tmudholk/data/high_stats_syst_opt/et150\\\",\\\"/home/tmudholk/research/hgcal_analysis/hgcal_optimization/analysis/macros\\\",150.0,%.1f,%.1f)\""%(rwcu_f,rwcu_m))
            os.system("root -b -q \"/home/tmudholk/research/hgcal_analysis/hgcal_optimization/analysis/macros/plotE_resolution_high_stats_syst_opt.C++(45,\\\"high_stats_syst_opt\\\",\\\"/home/tmudholk/data/high_stats_syst_opt\\\",\\\"/home/tmudholk/research/hgcal_analysis/hgcal_optimization/analysis/macros\\\",%.1f,%.1f)\""%(rwcu_f,rwcu_m))
