#!/bin/bash

cat plotE_resolution_drop_layers_modified_weights.C | sed "s|to_drop_array\[\] = .*;|to_drop_array[] = {25,26};|" > plotE_resolution_drop_layers_temp_modified_weights.C
mv plotE_resolution_drop_layers_temp_modified_weights.C plotE_resolution_drop_layers_modified_weights.C
root -b -q "plotE_resolution_drop_layers_modified_weights.C++(30,\"v30_modified_weights_dl_25_26\",\"\/home\/tmudholk\/data\/version30\",\"\/home\/tmudholk\/research\/hgcal_analysis\/hgcal_optimization\/analysis\/macros\")"

cat plotE_resolution_drop_layers_modified_weights.C | sed "s|to_drop_array\[\] = .*;|to_drop_array[] = {0,1};|" > plotE_resolution_drop_layers_temp_modified_weights.C
mv plotE_resolution_drop_layers_temp_modified_weights.C plotE_resolution_drop_layers_modified_weights.C
root -b -q "plotE_resolution_drop_layers_modified_weights.C++(30,\"v30_modified_weights_dl_0_1\",\"\/home\/tmudholk\/data\/version30\",\"\/home\/tmudholk\/research\/hgcal_analysis\/hgcal_optimization\/analysis\/macros\")"

cat plotE_resolution_drop_layers_modified_weights.C | sed "s|to_drop_array\[\] = .*;|to_drop_array[] = {0,26};|" > plotE_resolution_drop_layers_temp_modified_weights.C
mv plotE_resolution_drop_layers_temp_modified_weights.C plotE_resolution_drop_layers_modified_weights.C
root -b -q "plotE_resolution_drop_layers_modified_weights.C++(30,\"v30_modified_weights_dl_0_26\",\"\/home\/tmudholk\/data\/version30\",\"\/home\/tmudholk\/research\/hgcal_analysis\/hgcal_optimization\/analysis\/macros\")"

cat plotE_resolution_drop_layers_modified_weights.C | sed "s|to_drop_array\[\] = .*;|to_drop_array[] = {0,1,25,26};|" > plotE_resolution_drop_layers_temp_modified_weights.C
mv plotE_resolution_drop_layers_temp_modified_weights.C plotE_resolution_drop_layers_modified_weights.C
root -b -q "plotE_resolution_drop_layers_modified_weights.C++(30,\"v30_modified_weights_dl_0_1_25_26\",\"\/home\/tmudholk\/data\/version30\",\"\/home\/tmudholk\/research\/hgcal_analysis\/hgcal_optimization\/analysis\/macros\")"






cat plotE_resolution_drop_layers_x0_modified_weights.C | sed "s|to_drop_array\[\] = .*;|to_drop_array[] = {25,26};|" > plotE_resolution_drop_layers_x0_temp_modified_weights.C
mv plotE_resolution_drop_layers_x0_temp_modified_weights.C plotE_resolution_drop_layers_x0_modified_weights.C
root -b -q "plotE_resolution_drop_layers_x0_modified_weights.C++(30,\"v30_modified_weights_dl_x0_25_26\",\"\/home\/tmudholk\/data\/version30\",\"\/home\/tmudholk\/research\/hgcal_analysis\/hgcal_optimization\/analysis\/macros\")"

cat plotE_resolution_drop_layers_x0_modified_weights.C | sed "s|to_drop_array\[\] = .*;|to_drop_array[] = {0,1};|" > plotE_resolution_drop_layers_x0_temp_modified_weights.C
mv plotE_resolution_drop_layers_x0_temp_modified_weights.C plotE_resolution_drop_layers_x0_modified_weights.C
root -b -q "plotE_resolution_drop_layers_x0_modified_weights.C++(30,\"v30_modified_weights_dl_x0_0_1\",\"\/home\/tmudholk\/data\/version30\",\"\/home\/tmudholk\/research\/hgcal_analysis\/hgcal_optimization\/analysis\/macros\")"

cat plotE_resolution_drop_layers_x0_modified_weights.C | sed "s|to_drop_array\[\] = .*;|to_drop_array[] = {0,26};|" > plotE_resolution_drop_layers_x0_temp_modified_weights.C
mv plotE_resolution_drop_layers_x0_temp_modified_weights.C plotE_resolution_drop_layers_x0_modified_weights.C
root -b -q "plotE_resolution_drop_layers_x0_modified_weights.C++(30,\"v30_modified_weights_dl_x0_0_26\",\"\/home\/tmudholk\/data\/version30\",\"\/home\/tmudholk\/research\/hgcal_analysis\/hgcal_optimization\/analysis\/macros\")"

cat plotE_resolution_drop_layers_x0_modified_weights.C | sed "s|to_drop_array\[\] = .*;|to_drop_array[] = {0,1,25,26};|" > plotE_resolution_drop_layers_x0_temp_modified_weights.C
mv plotE_resolution_drop_layers_x0_temp_modified_weights.C plotE_resolution_drop_layers_x0_modified_weights.C
root -b -q "plotE_resolution_drop_layers_x0_modified_weights.C++(30,\"v30_modified_weights_dl_x0_0_1_25_26\",\"\/home\/tmudholk\/data\/version30\",\"\/home\/tmudholk\/research\/hgcal_analysis\/hgcal_optimization\/analysis\/macros\")"
