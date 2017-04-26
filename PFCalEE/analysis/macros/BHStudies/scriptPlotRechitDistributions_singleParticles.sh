#!/bin/bash

./plotRechitDistributions_singleParticles.py --inputDigiDirectory=/afs/cern.ch/work/t/tmudholk/public/bhStudies/muonGun/simulated/ --inputDigiFileNamePrefix=DiginoNoise_lowThreshold_IC3__version33_model2_BOFF_et100 --outputFileName=energyProfiles_muon --resetOutputDir --minHistEnergy=0.0 --maxHistEnergy=3.0 && ./plotRechitDistributions_singleParticles.py --inputDigiDirectory=/afs/cern.ch/work/t/tmudholk/public/bhStudies/pionGun/simulated/ --inputDigiFileNamePrefix=DiginoNoise_lowThreshold_IC3__version33_model2_BOFF_et100 --outputFileName=energyProfiles_pion --minHistEnergy=0.0 --maxHistEnergy=30.0 --nEnergyBins=300 --setLogX