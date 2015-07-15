export USERBASE=`pwd`
ARCH=x86_64-slc6-gcc46-opt
export BOOSTSYS=/cvmfs/cms.cern.ch/slc6_amd64_gcc472/external/boost/1.47.0-cms
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$USERBASE/userlib/lib:$USERBASE/analysis/lib:$BOOSTSYS/lib
source /afs/cern.ch/sw/lcg/contrib/gcc/4.6/x86_64-slc6/setup.sh
cd /afs/cern.ch/sw/lcg/app/releases/ROOT/5.34.18/${ARCH}/root/
source bin/thisroot.sh
cd - &> /dev/null
