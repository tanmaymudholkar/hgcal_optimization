#!/usr/bin/env python

from __future__ import print_function

import os

consolidated_data_file_path="/afs/cern.ch/user/t/tmudholk/public/research/hgcal_analysis/hgcal_optimization/analysis/macros/resolutions/data_consolidated"
consolidated_data_file=open(consolidated_data_file_path,'w')
for index_f in range(0,11):
    rwcu_f = index_f*0.1
    for index_m in range(0,11):
        rwcu_m = index_m*0.1
        if (not(rwcu_f==1.0 and rwcu_m==1.0)):
            inputFile_path="/afs/cern.ch/user/t/tmudholk/public/research/hgcal_analysis/hgcal_optimization/analysis/macros/resolutions/data_resolution_et150_rwcuf%.1f_rwcum%.1f"%(rwcu_f,rwcu_m)
            if (os.path.exists(inputFile_path)):
                inputFile=open(inputFile_path,'r')
                for line in inputFile:
                    consolidated_data_file.write(line)
                inputFile.close()
            else:
                print ("%s does not exist"%(inputFile_path))

consolidated_data_file.close()
