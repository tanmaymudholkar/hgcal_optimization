#!/usr/bin/env python

from __future__ import print_function

import os

# energy=150.0

# macros_directory = os.getcwd()
# consolidated_data_file_path="%s/resolutions/data_high_stats_syst_opt_et%.0f"%(macros_directory,energy)
# consolidated_data_file=open(consolidated_data_file_path,'w')
# for index_f in range(0,11):
#     rwcu_f = index_f*0.1
#     for index_m in range(0,11):
#         rwcu_m = index_m*0.1
#         if (not(rwcu_f==1.0 and rwcu_m==1.0)):
#             inputFile_path="%s/resolutions/data_resolution_et%.0f_rwcuf%.1f_rwcum%.1f"%(macros_directory,energy,rwcu_f,rwcu_m)
#             if (os.path.exists(inputFile_path)):
#                 inputFile=open(inputFile_path,'r')
#                 for line in inputFile:
#                     consolidated_data_file.write(line)
#                 inputFile.close()
#             else:
#                 print ("%s does not exist"%(inputFile_path))

# consolidated_data_file.close()


macros_directory = os.getcwd()
consolidated_data_file_path="%s/resolutions/data_consolidated_high_stats_syst_opt"%(macros_directory)
consolidated_data_file=open(consolidated_data_file_path,'w')
for index_f in range(0,11):
    rwcu_f = index_f*0.1
    for index_m in range(0,11):
        rwcu_m = index_m*0.1
        if (not(rwcu_f==1.0 and rwcu_m==1.0)):
            inputFile_path="%s/resolutions/data_resolutions_high_stats_syst_opt_eta2.100_rwcuf%.1f_rwcum%.1f"%(macros_directory,rwcu_f,rwcu_m)
            if (os.path.exists(inputFile_path)):
                inputFile=open(inputFile_path,'r')
                for line in inputFile:
                    consolidated_data_file.write("%.1f    %.1f    %s"%(rwcu_f,rwcu_m,line))
                inputFile.close()
            else:
                print ("%s does not exist"%(inputFile_path))

consolidated_data_file.close()
