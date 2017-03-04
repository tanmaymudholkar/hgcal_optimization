#!/bin/bash
source /afs/cern.ch/user/t/tmudholk/public/research/hgcal_minbias_modified_geometry/PFCal/PFCalEE/userlib/../g4env.sh
localdir=`pwd`
/afs/cern.ch/user/t/tmudholk/public/research/hgcal_minbias_modified_geometry/PFCal/PFCalEE/userlib/bin/digitizer 5 root://eoscms//eos/cms/store/cmst3/group/hgcal/BHStudies_Standalone/HepMCTest_VBFHtoBB/PythiaTest/HGcal_version33_modifiedGeometry.root $localdir/ 0-27:4,28-39:4,40-51:4 0-39:0,40-51:0 0-9:4.1,10-51:5 3 2 0 root://eoscms//eos/cms/store/cmst3/group/hgcal/HGCalMinbias/PythiaTest/ | tee /afs/cern.ch/work/t/tmudholk/public/bhStudies/HepMCTest_VBFHtoBB/version_33/PythiaTest//digitizernoiseOff_IC3.log
echo "--Local directory is " $localdir >> digijobnoiseOff_IC3.log
ls * >> digijobnoiseOff_IC3.log
grep "alias eos=" /afs/cern.ch/project/eos/installation/cms/etc/setup.sh | sed "s/alias /export my/" > eosenv.sh
source eosenv.sh
$myeos mkdir -p /store/cmst3/group/hgcal/BHStudies_Standalone/HepMCTest_VBFHtoBB/PythiaTest
cmsStage -f DigiPFcal.root /store/cmst3/group/hgcal/BHStudies_Standalone/HepMCTest_VBFHtoBB/PythiaTest/DiginoiseOff_IC3_version33_modifiedGeometry.root
if (( "$?" != "0" )); then
echo " --- Problem with copy of file DigiPFcal.root to EOS. Keeping locally." >> digijobnoiseOff_IC3.log
else
eossize=`$myeos ls -l /store/cmst3/group/hgcal/BHStudies_Standalone/HepMCTest_VBFHtoBB/PythiaTest/DiginoiseOff_IC3_version33_modifiedGeometry.root | awk '{print $5}'`
localsize=`ls -l DigiPFcal.root | awk '{print $5}'`
if (( "$eossize" != "$localsize" )); then
echo " --- Copy of digi file to eos failed. Localsize = $localsize, eossize = $eossize. Keeping locally..." >> digijobnoiseOff_IC3.log
else
echo " --- Size check done: Localsize = $localsize, eossize = $eossize" >> digijobnoiseOff_IC3.log
echo " --- File DigiPFcal.root successfully copied to EOS: /store/cmst3/group/hgcal/BHStudies_Standalone/HepMCTest_VBFHtoBB/PythiaTest/DiginoiseOff_IC3_version33_modifiedGeometry.root" >> digijobnoiseOff_IC3.log
rm DigiPFcal.root
fi
fi
echo "--deleting core files: too heavy!!"
rm core.*
cp * /afs/cern.ch/work/t/tmudholk/public/bhStudies/HepMCTest_VBFHtoBB/version_33/PythiaTest//
echo "All done"
