#!/bin/bash
source /afs/cern.ch/user/t/tmudholk/.bashrc
cp /afs/cern.ch/user/t/tmudholk/public/research/hgcal_analysis/hgcal_optimization/analysis/macros/plotE_resolution_high_stats_syst_opt.C /afs/cern.ch/user/t/tmudholk/public/research/hgcal_analysis/hgcal_optimization/analysis/macros/scripts/plotE_resolution_high_stats_syst_opt_rwcuf9_rwcum4_temp.C
sed "s|void plotE_resolution_high_stats_syst_opt|void plotE_resolution_high_stats_syst_opt_rwcuf9_rwcum4|" /afs/cern.ch/user/t/tmudholk/public/research/hgcal_analysis/hgcal_optimization/analysis/macros/scripts/plotE_resolution_high_stats_syst_opt_rwcuf9_rwcum4_temp.C > /afs/cern.ch/user/t/tmudholk/public/research/hgcal_analysis/hgcal_optimization/analysis/macros/scripts/plotE_resolution_high_stats_syst_opt_rwcuf9_rwcum4.C
root -b -q "/afs/cern.ch/user/t/tmudholk/public/research/hgcal_analysis/hgcal_optimization/analysis/macros/scripts/plotE_resolution_high_stats_syst_opt_rwcuf9_rwcum4.C++(45,\"high_stats_syst_opt\",\"/afs/cern.ch/user/t/tmudholk/private/mount/eos_cern/cms/store/user/tmudholk/high_stats_syst_opt/gitV05-02-04/gamma\",\"/afs/cern.ch/user/t/tmudholk/public/research/hgcal_analysis/hgcal_optimization/analysis/macros\",150.0,0.9,0.4)"
source /afs/cern.ch/user/t/tmudholk/.bash_logout
