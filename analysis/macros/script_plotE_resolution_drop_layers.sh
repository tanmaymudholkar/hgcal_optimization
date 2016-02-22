#!/bin/bash

cat plotE_resolution_drop_layers.C | sed "s|to_drop_array\[\] = .*;|to_drop_array[] = {0,1};|" > plotE_resolution_drop_layers_temp.C
mv plotE_resolution_drop_layers_temp.C plotE_resolution_drop_layers.C
root -b -q "plotE_resolution_drop_layers.C++(30,\"v30_dl_0_1\",\"export_dir\/version30\")"

cat plotE_resolution_drop_layers.C | sed "s|to_drop_array\[\] = .*;|to_drop_array[] = {0,1,25,26};|" > plotE_resolution_drop_layers_temp.C
mv plotE_resolution_drop_layers_temp.C plotE_resolution_drop_layers.C
root -b -q "plotE_resolution_drop_layers.C++(30,\"v30_dl_0_1_25_26\",\"export_dir\/version30\")"

cat plotE_resolution_drop_layers.C | sed "s|to_drop_array\[\] = .*;|to_drop_array[] = {0,1,23,24,25,26};|" > plotE_resolution_drop_layers_temp.C
mv plotE_resolution_drop_layers_temp.C plotE_resolution_drop_layers.C
root -b -q "plotE_resolution_drop_layers.C++(30,\"v30_dl_0_1_23_24_25_26\",\"export_dir\/version30\")"

cat plotE_resolution_drop_layers.C | sed "s|to_drop_array\[\] = .*;|to_drop_array[] = {0,1,2,3,25,26};|" > plotE_resolution_drop_layers_temp.C
mv plotE_resolution_drop_layers_temp.C plotE_resolution_drop_layers.C
root -b -q "plotE_resolution_drop_layers.C++(30,\"v30_dl_0_1_2_3_25_26\",\"export_dir\/version30\")"

cat plotE_resolution_drop_layers.C | sed "s|to_drop_array\[\] = .*;|to_drop_array[] = {0,1,2,3,23,24,25,26};|" > plotE_resolution_drop_layers_temp.C
mv plotE_resolution_drop_layers_temp.C plotE_resolution_drop_layers.C
root -b -q "plotE_resolution_drop_layers.C++(30,\"v30_dl_0_1_2_3_23_24_25_26\",\"export_dir\/version30\")"

cat plotE_resolution_drop_layers.C | sed "s|to_drop_array\[\] = .*;|to_drop_array[] = {0,1,2,3,18,19,23,24,25,26};|" > plotE_resolution_drop_layers_temp.C
mv plotE_resolution_drop_layers_temp.C plotE_resolution_drop_layers.C
root -b -q "plotE_resolution_drop_layers.C++(30,\"v30_dl_0_1_2_3_18_19_23_24_25_26\",\"export_dir\/version30\")"

cat plotE_resolution_drop_layers.C | sed "s|to_drop_array\[\] = .*;|to_drop_array[] = {0,1,2,3,7,8,18,19,23,24,25,26};|" > plotE_resolution_drop_layers_temp.C
mv plotE_resolution_drop_layers_temp.C plotE_resolution_drop_layers.C
root -b -q "plotE_resolution_drop_layers.C++(30,\"v30_dl_0_1_2_3_7_8_18_19_23_24_25_26\",\"export_dir\/version30\")"

cat plotE_resolution_drop_layers.C | sed "s|to_drop_array\[\] = .*;|to_drop_array[] = {0,1,2,3,7,8,13,14,18,19,23,24,25,26};|" > plotE_resolution_drop_layers_temp.C
mv plotE_resolution_drop_layers_temp.C plotE_resolution_drop_layers.C
root -b -q "plotE_resolution_drop_layers.C++(30,\"v30_dl_0_1_2_3_7_8_13_14_18_19_23_24_25_26\",\"export_dir\/version30\")"





cat plotE_resolution_drop_layers_x0.C | sed "s|to_drop_array\[\] = .*;|to_drop_array[] = {25,26};|" > plotE_resolution_drop_layers_x0_temp.C
mv plotE_resolution_drop_layers_x0_temp.C plotE_resolution_drop_layers_x0.C
root -b -q "plotE_resolution_drop_layers_x0.C++(30,\"v30_dl_x0_25_26\",\"export_dir\/version30\")"

cat plotE_resolution_drop_layers_x0.C | sed "s|to_drop_array\[\] = .*;|to_drop_array[] = {0,1};|" > plotE_resolution_drop_layers_x0_temp.C
mv plotE_resolution_drop_layers_x0_temp.C plotE_resolution_drop_layers_x0.C
root -b -q "plotE_resolution_drop_layers_x0.C++(30,\"v30_dl_x0_0_1\",\"export_dir\/version30\")"

cat plotE_resolution_drop_layers_x0.C | sed "s|to_drop_array\[\] = .*;|to_drop_array[] = {0,1,25,26};|" > plotE_resolution_drop_layers_x0_temp.C
mv plotE_resolution_drop_layers_x0_temp.C plotE_resolution_drop_layers_x0.C
root -b -q "plotE_resolution_drop_layers_x0.C++(30,\"v30_dl_x0_0_1_25_26\",\"export_dir\/version30\")"

cat plotE_resolution_drop_layers_x0.C | sed "s|to_drop_array\[\] = .*;|to_drop_array[] = {0,1,23,24,25,26};|" > plotE_resolution_drop_layers_x0_temp.C
mv plotE_resolution_drop_layers_x0_temp.C plotE_resolution_drop_layers_x0.C
root -b -q "plotE_resolution_drop_layers_x0.C++(30,\"v30_dl_x0_0_1_23_24_25_26\",\"export_dir\/version30\")"

cat plotE_resolution_drop_layers_x0.C | sed "s|to_drop_array\[\] = .*;|to_drop_array[] = {0,1,2,3,25,26};|" > plotE_resolution_drop_layers_x0_temp.C
mv plotE_resolution_drop_layers_x0_temp.C plotE_resolution_drop_layers_x0.C
root -b -q "plotE_resolution_drop_layers_x0.C++(30,\"v30_dl_x0_0_1_2_3_25_26\",\"export_dir\/version30\")"

cat plotE_resolution_drop_layers_x0.C | sed "s|to_drop_array\[\] = .*;|to_drop_array[] = {0,1,2,3,23,24,25,26};|" > plotE_resolution_drop_layers_x0_temp.C
mv plotE_resolution_drop_layers_x0_temp.C plotE_resolution_drop_layers_x0.C
root -b -q "plotE_resolution_drop_layers_x0.C++(30,\"v30_dl_x0_0_1_2_3_23_24_25_26\",\"export_dir\/version30\")"

cat plotE_resolution_drop_layers_x0.C | sed "s|to_drop_array\[\] = .*;|to_drop_array[] = {0,1,2,3,18,19,23,24,25,26};|" > plotE_resolution_drop_layers_x0_temp.C
mv plotE_resolution_drop_layers_x0_temp.C plotE_resolution_drop_layers_x0.C
root -b -q "plotE_resolution_drop_layers_x0.C++(30,\"v30_dl_x0_0_1_2_3_18_19_23_24_25_26\",\"export_dir\/version30\")"

cat plotE_resolution_drop_layers_x0.C | sed "s|to_drop_array\[\] = .*;|to_drop_array[] = {0,1,2,3,7,8,18,19,23,24,25,26};|" > plotE_resolution_drop_layers_x0_temp.C
mv plotE_resolution_drop_layers_x0_temp.C plotE_resolution_drop_layers_x0.C
root -b -q "plotE_resolution_drop_layers_x0.C++(30,\"v30_dl_x0_0_1_2_3_7_8_18_19_23_24_25_26\",\"export_dir\/version30\")"

cat plotE_resolution_drop_layers_x0.C | sed "s|to_drop_array\[\] = .*;|to_drop_array[] = {0,1,2,3,7,8,13,14,18,19,23,24,25,26};|" > plotE_resolution_drop_layers_x0_temp.C
mv plotE_resolution_drop_layers_x0_temp.C plotE_resolution_drop_layers_x0.C
root -b -q "plotE_resolution_drop_layers_x0.C++(30,\"v30_dl_x0_0_1_2_3_7_8_13_14_18_19_23_24_25_26\",\"export_dir\/version30\")"

