#!/bin/bash

g++ -g -Wall -W -O2 -std=c++0x `root-config --cflags --glibs` -I/include/ -isystem /cvmfs/cms.cern.ch/slc6_amd64_gcc472/external/boost/1.47.0-cms/include -I../../../userlib/include/ -I/afs/cern.ch/sw/lcg/external/fastjet/3.0.3/x86_64-slc6-gcc46-opt//include/ -I/afs/cern.ch/sw/lcg/external/HepMC/2.06.08/x86_64-slc6-gcc46-opt//include/ -fPIC -c getXYPositions.C -o getXYPositions.o &&

g++ -shared -Wall -W -o getXYPositionsLib.so getXYPositions.o -L/afs/cern.ch/sw/lcg/app/releases/ROOT/5.34.18/x86_64-slc6-gcc46-opt/root/lib -lGui -lCore -lCint -lRIO -lNet -lHist -lGraf -lGraf3d -lGpad -lTree -lRint -lPostscript -lMatrix -lPhysics -lMathCore -lThread -pthread -lm -ldl -rdynamic -lRooFit -lGenVector  -L../../../userlib/lib -lPFCalEEuserlib -Wl,-rpath,/afs/cern.ch/sw/lcg/external/fastjet/3.0.3/x86_64-slc6-gcc46-opt//lib -lm  -L/afs/cern.ch/sw/lcg/external/fastjet/3.0.3/x86_64-slc6-gcc46-opt//lib -lfastjettools -lfastjet -lfastjetplugins -lsiscone_spherical -lsiscone -L/cvmfs/cms.cern.ch/slc6_amd64_gcc472/external/boost/1.47.0-cms/lib -lboost_regex -lboost_program_options -lboost_filesystem &&

g++ -o getXYPositions.out -g -Wall -W -O2 -std=c++0x  -I/afs/cern.ch/sw/lcg/app/releases/ROOT/5.34.18/x86_64-slc6-gcc46-opt/root/include/ -I/include/ -isystem /cvmfs/cms.cern.ch/slc6_amd64_gcc472/external/boost/1.47.0-cms/include -I../../../userlib/include/ -I/afs/cern.ch/sw/lcg/external/fastjet/3.0.3/x86_64-slc6-gcc46-opt//include/ -I/afs/cern.ch/sw/lcg/external/HepMC/2.06.08/x86_64-slc6-gcc46-opt//include/ getXYPositions.C -L/afs/cern.ch/sw/lcg/app/releases/ROOT/5.34.18/x86_64-slc6-gcc46-opt/root/lib -lGui -lCore -lCint -lRIO -lNet -lHist -lGraf -lGraf3d -lGpad -lTree -lRint -lPostscript -lMatrix -lPhysics -lMathCore -lThread -pthread -lm -ldl -rdynamic -lRooFit -lGenVector  -L../../../userlib/lib -lPFCalEEuserlib -Wl,-rpath,/afs/cern.ch/sw/lcg/external/fastjet/3.0.3/x86_64-slc6-gcc46-opt//lib -lm  -L/afs/cern.ch/sw/lcg/external/fastjet/3.0.3/x86_64-slc6-gcc46-opt//lib -lfastjettools -lfastjet -lfastjetplugins -lsiscone_spherical -lsiscone -L/cvmfs/cms.cern.ch/slc6_amd64_gcc472/external/boost/1.47.0-cms/lib -lboost_regex -lboost_program_options -lboost_filesystem -L.-lgetXYPositionsLib

rm getXYPositions.o
rm getXYPositionsLib.so