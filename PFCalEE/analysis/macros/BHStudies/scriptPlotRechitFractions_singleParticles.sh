#!/bin/bash

./plotRechitFractions_singleParticles.py --inputDigiDirectory=/afs/cern.ch/work/t/tmudholk/public/bhStudies/pionGun/simulated/ --inputDigiFileNamePrefix=DiginoNoise_lowThreshold_withExtraInfo_IC3__version33_model2_BOFF_et100 --outputFileNamePrefix=pionGun_Scint --resetOutputDir --minHistEnergy=0.0 --maxHistEnergy=30.0 --nEnergyBins=300 --setLogX --forceNoTitle

./plotRechitFractions_singleParticles.py --inputDigiDirectory=/afs/cern.ch/work/t/tmudholk/public/bhStudies/pionGun/simulated/ --inputDigiFileNamePrefix=DiginoNoise_lowThreshold_withExtraInfo_IC3__version33_model2_BOFF_et100 --outputFileNamePrefix=pionGun_Si --minHistEnergy=0.0 --maxHistEnergy=30.0 --nEnergyBins=300 --setLogX --startLayer=0 --stopLayer=27 --forceNoTitle

./plotRechitFractions_singleParticles.py --inputDigiDirectory=/afs/cern.ch/work/t/tmudholk/public/bhStudies/muonGun/simulated/ --inputDigiFileNamePrefix=DiginoNoise_lowThreshold_withExtraInfo_IC3__version33_model2_BOFF_et100 --outputFileNamePrefix=muonGun_Scint --minHistEnergy=0.0 --maxHistEnergy=3.0 --boundaryRechitEnergyInMips=1.0 --forceNoTitle

./plotRechitFractions_singleParticles.py --inputDigiDirectory=/afs/cern.ch/work/t/tmudholk/public/bhStudies/muonGun/simulated/ --inputDigiFileNamePrefix=DiginoNoise_lowThreshold_withExtraInfo_IC3__version33_model2_BOFF_et100 --outputFileNamePrefix=muonGun_Si --minHistEnergy=0.0 --maxHistEnergy=3.0 --startLayer=0 --stopLayer=27 --boundaryRechitEnergyInMips=1.0 --forceNoTitle
