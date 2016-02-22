#!/bin/bash

cat plotE_resolution_high_stats.C | sed "s|calibration_offset = .*;|calibration_offset = 11.3;|" > plotE_resolution_high_stats_temp.C
mv plotE_resolution_high_stats_temp.C plotE_resolution_high_stats.C
cat plotE_resolution_high_stats.C | sed "s|calibration_offset_error = .*;|calibration_offset_error = 1.0;|" > plotE_resolution_high_stats_temp.C
mv plotE_resolution_high_stats_temp.C plotE_resolution_high_stats.C
cat plotE_resolution_high_stats.C | sed "s|calibration_slope = .*;|calibration_slope = 46.600;|" > plotE_resolution_high_stats_temp.C
mv plotE_resolution_high_stats_temp.C plotE_resolution_high_stats.C
cat plotE_resolution_high_stats.C | sed "s|calibration_slope_error = .*;|calibration_slope_error = 0.008;|" > plotE_resolution_high_stats_temp.C
mv plotE_resolution_high_stats_temp.C plotE_resolution_high_stats.C

root -b -q "plotE_resolution_high_stats.C++(34,\"high_stats_v34\",\"export_dir\/high_stats_version_34\")"

cat plotE_resolution_high_stats.C | sed "s|calibration_offset = .*;|calibration_offset = 15.4;|" > plotE_resolution_high_stats_temp.C
mv plotE_resolution_high_stats_temp.C plotE_resolution_high_stats.C
cat plotE_resolution_high_stats.C | sed "s|calibration_offset_error = .*;|calibration_offset_error = 1.1;|" > plotE_resolution_high_stats_temp.C
mv plotE_resolution_high_stats_temp.C plotE_resolution_high_stats.C
cat plotE_resolution_high_stats.C | sed "s|calibration_slope = .*;|calibration_slope = 43.708;|" > plotE_resolution_high_stats_temp.C
mv plotE_resolution_high_stats_temp.C plotE_resolution_high_stats.C
cat plotE_resolution_high_stats.C | sed "s|calibration_slope_error = .*;|calibration_slope_error = 0.012;|" > plotE_resolution_high_stats_temp.C
mv plotE_resolution_high_stats_temp.C plotE_resolution_high_stats.C

root -b -q "plotE_resolution_high_stats.C++(44,\"high_stats_vflat\",\"export_dir\/high_stats_vflat\")"
