#!/usr/bin/env python

from __future__ import print_function

import os

for index_f in range(0,11):
    rwcu_f = index_f*0.1
    for index_m in range(0,11):
        rwcu_m = index_m*0.1
        if(not(rwcu_f==1.0 and rwcu_m==1.0)):
            scriptFile_path="/afs/cern.ch/user/t/tmudholk/public/research/hgcal_analysis/hgcal_optimization/analysis/macros/scripts/script_plotE_resolution_highstats_syst_opt_rwcuf_%.1f_rwcum_%.1f.sh"%(rwcu_f,rwcu_m)
            scriptFile=open(scriptFile_path,'w')
            scriptFile.write("#!/bin/bash\n")
            scriptFile.write("source /afs/cern.ch/user/t/tmudholk/.bashrc\n")
            scriptFile.write("cp /afs/cern.ch/user/t/tmudholk/public/research/hgcal_analysis/hgcal_optimization/analysis/macros/plotE_resolution_high_stats_syst_opt.C /afs/cern.ch/user/t/tmudholk/public/research/hgcal_analysis/hgcal_optimization/analysis/macros/scripts/plotE_resolution_high_stats_syst_opt_rwcuf%i_rwcum%i_temp.C\n"%(index_f,index_m))
            scriptFile.write("sed \"s|void plotE_resolution_high_stats_syst_opt|void plotE_resolution_high_stats_syst_opt_rwcuf%i_rwcum%i|\" /afs/cern.ch/user/t/tmudholk/public/research/hgcal_analysis/hgcal_optimization/analysis/macros/scripts/plotE_resolution_high_stats_syst_opt_rwcuf%i_rwcum%i_temp.C > /afs/cern.ch/user/t/tmudholk/public/research/hgcal_analysis/hgcal_optimization/analysis/macros/scripts/plotE_resolution_high_stats_syst_opt_rwcuf%i_rwcum%i.C\n"%(index_f,index_m,index_f,index_m,index_f,index_m))
            scriptFile.write("root -b -q \"/afs/cern.ch/user/t/tmudholk/public/research/hgcal_analysis/hgcal_optimization/analysis/macros/scripts/plotE_resolution_high_stats_syst_opt_rwcuf%i_rwcum%i.C++(45,\\\"high_stats_syst_opt\\\",\\\"/afs/cern.ch/user/t/tmudholk/private/mount/eos_cern/cms/store/user/tmudholk/high_stats_syst_opt/gitV05-02-04/gamma\\\",\\\"/afs/cern.ch/user/t/tmudholk/public/research/hgcal_analysis/hgcal_optimization/analysis/macros\\\",150.0,%.1f,%.1f)\"\n"%(index_f,index_m,rwcu_f,rwcu_m))
            scriptFile.write("source /afs/cern.ch/user/t/tmudholk/.bash_logout\n")
            scriptFile.close()
            os.system("chmod u+rwx %s"%(scriptFile_path))
