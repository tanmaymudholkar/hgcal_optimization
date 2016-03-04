#!/usr/bin/env python

from __future__ import print_function

import os

os.system("./script_generate_scripts_highstats_syst_opt.py")

for index_f in range(0,11):
    rwcu_f = index_f*0.1
    for index_m in range(0,11):
        rwcu_m = index_m*0.1
        if (not(rwcu_f==1.0 and rwcu_m==1.0)):
            scriptFile_path="/afs/cern.ch/user/t/tmudholk/public/research/hgcal_analysis/hgcal_optimization/analysis/macros/scripts/script_plotE_resolution_highstats_syst_opt_rwcuf_%.1f_rwcum_%.1f.sh"%(rwcu_f,rwcu_m)
            os.system("echo \"bsub -q 1nh %s\""%(scriptFile_path))
