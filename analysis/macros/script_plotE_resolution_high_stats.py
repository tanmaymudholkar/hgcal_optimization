#!/usr/bin/env python

from __future__ import print_function

import os

# os.system("root -b -q \"/home/tmudholk/research/hgcal_analysis/hgcal_optimization/analysis/macros/plotE_resolution_high_stats.C++(34,\\\"high_stats_version34\\\",\\\"/home/tmudholk/data/high_stats_version34\\\",\\\"/home/tmudholk/research/hgcal_analysis/hgcal_optimization/analysis/macros\\\")\"")
# os.system("root -b -q \"/home/tmudholk/research/hgcal_analysis/hgcal_optimization/analysis/macros/plotE_resolution_high_stats.C++(44,\\\"high_stats_vflat\\\",\\\"/home/tmudholk/data/high_stats_vflat\\\",\\\"/home/tmudholk/research/hgcal_analysis/hgcal_optimization/analysis/macros\\\")\"")

os.system("root -b -q \"/home/tmudholk/research/hgcal_analysis/hgcal_optimization/analysis/macros/plotE_resolution_high_stats.C++(33,\\\"version33\\\",\\\"/home/tmudholk/data/hexagonal_geometry/version_33\\\",\\\"/home/tmudholk/research/hgcal_analysis/hgcal_optimization/analysis/macros\\\",0.5)\"")

os.system("root -b -q \"/home/tmudholk/research/hgcal_analysis/hgcal_optimization/analysis/macros/plotE_resolution_high_stats.C++(33,\\\"version33\\\",\\\"/home/tmudholk/data/hexagonal_geometry/version_33\\\",\\\"/home/tmudholk/research/hgcal_analysis/hgcal_optimization/analysis/macros\\\",2.0)\"")

os.system("root -b -q \"/home/tmudholk/research/hgcal_analysis/hgcal_optimization/analysis/macros/plotE_resolution_high_stats.C++(33,\\\"version33\\\",\\\"/home/tmudholk/data/hexagonal_geometry/version_33\\\",\\\"/home/tmudholk/research/hgcal_analysis/hgcal_optimization/analysis/macros\\\",5.0)\"")
