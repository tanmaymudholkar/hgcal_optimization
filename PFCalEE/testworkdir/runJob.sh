#!/bin/bash
source /afs/cern.ch/user/t/tmudholk/public/research/hgcal_minbias_modified_geometry/PFCal/PFCalEE/g4env.sh
cp /afs/cern.ch/work/t/tmudholk/public/bhStudies/HepMCTest_VBFHtoBB/version_33/PythiaTest/modifiedGeometry//g4steer.mac .
PFCalEE g4steer.mac 33 2
localdir=`pwd`
echo "--Local directory is " $localdir > g4.log
ls * >> g4.log
grep "alias eos=" /afs/cern.ch/project/eos/installation/cms/etc/setup.sh | sed "s/alias /export my/" > eosenv.sh
source eosenv.sh
$myeos mkdir -p /store/cmst3/group/hgcal/BHStudies_Standalone/HepMCTest_VBFHtoBB/PythiaTest
cmsStage -f PFcal.root /store/cmst3/group/hgcal/BHStudies_Standalone/HepMCTest_VBFHtoBB/PythiaTest/HGcal_version33_modifiedGeometry.root
if (( "$?" != "0" )); then
echo " --- Problem with copy of file PFcal.root to EOS. Keeping locally." >> g4.log
else
echo " --- File PFcal.root successfully copied to EOS: /store/cmst3/group/hgcal/BHStudies_Standalone/HepMCTest_VBFHtoBB/PythiaTest/HGcal_version33_modifiedGeometry.root" >> g4.log
rm PFcal.root
fi
cp * /afs/cern.ch/work/t/tmudholk/public/bhStudies/HepMCTest_VBFHtoBB/version_33/PythiaTest/modifiedGeometry//
echo "All done"
