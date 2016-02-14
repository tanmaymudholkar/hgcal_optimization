#!/usr/bin/env python

import os,sys
import optparse
import commands
import math
import random

random.seed()

usage = 'usage: %prog [options]'
parser = optparse.OptionParser(usage)
parser.add_option('-s', '--short-queue',    dest='squeue'             , help='short batch queue'            , default='1nd')
parser.add_option('-q', '--long-queue' ,    dest='lqueue'             , help='long batch queue'             , default='2nw')
parser.add_option('-t', '--git-tag'     ,    dest='gittag'             , help='git tag version'              , default='V00-00-00')
parser.add_option('-r', '--run'         ,    dest='run'                , help='stat run'                     , default=-1,      type=int)
parser.add_option('-v', '--version'     ,    dest='version'            , help='detector version'             , default=3,      type=int)
parser.add_option('-m', '--model'       ,    dest='model'              , help='detector model'               , default=3,      type=int)
parser.add_option('-a', '--eta'       ,    dest='eta'              , help='incidence eta'       , default=0,      type=float)
parser.add_option('-T', '--threshold'   ,    dest='threshold'          , help='noise threshold'     , default=5,      type=float)
parser.add_option('-b', '--Bfield'      ,    dest='Bfield'             , help='B field value in Tesla'       , default=0,      type=float)
parser.add_option('-d', '--datatype'    ,    dest='datatype'           , help='data type or particle to shoot', default='e-')
parser.add_option('-f', '--datafile'    ,    dest='datafile'           , help='full path to HepMC input file', default='data/example_MyPythia.dat')
parser.add_option('-n', '--nevts'       ,    dest='nevts'              , help='number of events to generate' , default=1000,    type=int)
parser.add_option('-o', '--out'         ,    dest='out'                , help='output directory'             , default=os.getcwd() )
parser.add_option('-w', '--wcu-radlengths-f'   ,    dest='rwcuf'          , help='total radiation lengths of wcu in front ecal'     , default=0,      type=float)
parser.add_option('-W', '--wcu-radlengths-m'   ,    dest='rwcum'          , help='total radiation lengths of wcu in middle ecal'     , default=0,      type=float)
# parser.add_option('-R', '--runnumber'   ,    dest='runno'              , help='run number'                   , default=1, type=int)
# parser.add_option('-e', '--eos'         ,    dest='eos'                , help='eos path to save root file to EOS',         default='')
# parser.add_option('-E', '--eosin'       ,    dest='eosin'              , help='eos path to read input root file from EOS',  default='')
parser.add_option('-g', '--gun'         ,    action="store_true",  dest='dogun'              , help='use particle gun.')
parser.add_option('-S', '--no-submit'   ,    action="store_true",  dest='nosubmit'           , help='Do not submit batch job.')
(opt, args) = parser.parse_args()

twcu_t = 39.6
tw_t = 41.6
x0wcu = 5.1223
x0w = 3.50418
rwcu_t = twcu_t/x0wcu
rw_t = tw_t/x0w

x0air = 303921
x0pcb = 187.31
x0si = 93.6607
x0cu = 14.3558
x0cfm = 3494.61
lcueven = 1.0
lcfm = 1.0
lair = 3.0
lpcb = 2.0
lsi = 0.3
lcuodd = 6.0;

rw_min = lcueven/x0cu + lcfm/x0cfm + lair/x0air + lpcb/x0pcb + lsi/x0si
rwcu_min = lcuodd/x0cu + lsi/x0si + lpcb/x0pcb + lair/x0air

rwcu_f=opt.rwcuf
rw_f = -rw_min + rwcu_min + rwcu_f
twcu_f = 0.5*rwcu_f*x0wcu
tw_f = rw_f*x0w
rwcu_m=opt.rwcum
rw_m = -rw_min + rwcu_min + rwcu_m
twcu_m = 0.5*rwcu_m*x0wcu
tw_m = rw_m*x0w
rwcu_b = 0.25*rwcu_t - rwcu_f - rwcu_m
rw_b = 0.25*rw_t - rw_f - rw_m
twcu_b = 0.5*rwcu_b*x0wcu
tw_b = rw_b*x0w

wthick='%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f'%(tw_f,tw_f,tw_f,tw_f,tw_m,tw_m,tw_m,tw_m,tw_b,tw_b,tw_b,tw_b)
wcuthick='%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f'%(twcu_f,twcu_f,twcu_f,twcu_f,twcu_m,twcu_m,twcu_m,twcu_m,twcu_b,twcu_b,twcu_b,twcu_b)
pbthick='1,1,1,1,1,2.1,2.1,2.1,2.1,2.1,4.4,4.4,4.4,4.4'

print "wthick is ",wthick
print "sumw is ",4*(tw_f+tw_m+tw_b)
print "wcuthick is ",wcuthick
print "sumwcu is ",8*(twcu_f+twcu_m+twcu_b)

enlist=[0]
if opt.dogun :
    # enlist=[3,7,20,35,50]
    # enlist=[70,100,125]
    enlist=[5,50]
    # enlist=[35]

# wthick='1.75,1.75,1.75,1.75,1.75,2.8,2.8,2.8,2.8,2.8,4.2,4.2,4.2,4.2,4.2'
# pbthick='1,1,1,1,1,2.1,2.1,2.1,2.1,2.1,4.4,4.4,4.4,4.4'
droplayers=''
label=''

INPATHPU="root://eoscms//eos/cms/store/cmst3/group/hgcal/Standalone/V12/MinBias/"

granularity='0-29:4,30-65:4'
noise='0-65:0.15'
threshold='0-65:5'

if (opt.version==8) :
    granularity='0-20:4,21-30:6'
    noise='0-30:0.14'
    threshold='0-30:2'
elif opt.version<20 :
    granularity='0-19:4,20-29:4'
    noise='0-29:0.14'
    threshold='0-29:5'
elif (opt.version==21 or opt.version==24):
    granularity='0-23:6,24-33:8'
    noise='0-33:0.14'
    threshold='0-33:2'
elif opt.version==22:
    granularity='0-9:8'
    noise='0-9:0.14'
    threshold='0-9:2'
elif opt.version==23:
    granularity='0-53:12'
    noise='0-53:0.14'
    threshold='0-53:2'
elif (opt.version==25 or opt.version==26):
    granularity='0-29:4,30-41:4,42-53:8'
    noise='0-41:0.14,42-53:0.2'
    threshold='0-53:5'
elif (opt.version==30):
    granularity='0-27:4'
    noise='0-27:0.14'
    threshold='0-27:5'
elif (opt.version==33):
    granularity='0-27:4,28-39:4,40-51:8'
    noise='0-39:0.14,40-51:0.2'
    threshold='0-51:5'
elif (opt.version==27 or opt.version==31):
    granularity='0-11:4,12-23:8'
    noise='0-11:0.14,12-23:0.2'
    threshold='0-23:5'
elif (opt.version==28 or opt.version==32):
    granularity='0-11:8'
    noise='0-11:0.2'
    threshold='0-11:5'
elif (opt.version==34 or opt.version==40 or opt.version==41 or opt.version==42 or opt.version==43 or opt.version==44 or opt.version==45):
    granularity='0-23:4'
    noise='0-23:0.14'
    threshold='0-23:%f'%(opt.threshold)
# elif (opt.version==40):
#     granularity='0-23:4'
#     noise='0-23:0.14'
#     threshold='0-23:%f'%(opt.threshold)
# elif (opt.version==41):
#     granularity='0-23:4'
#     noise='0-23:0.14'
#     threshold='0-23:%f'%(opt.threshold)
# elif (opt.version==42):
#     granularity='0-23:4'
#     noise='0-23:0.14'
#     threshold='0-23:%f'%(opt.threshold)
# elif (opt.version==43):
#     granularity='0-23:4'
#     noise='0-23:0.14'
#     threshold='0-23:%f'%(opt.threshold)
# elif (opt.version==44):
#     granularity='0-23:4'
#     noise='0-23:0.14'
#     threshold='0-23:%f'%(opt.threshold)
# elif (opt.version==45):
#     granularity='0-23:4'
#     noise='0-23:0.14'
#     threshold='0-23:%f'%(opt.threshold)
elif (opt.version==36):
    granularity='0-23:4,24-34:4,35-46:8'
    noise='0-34:0.14,35-46:0.2'
    threshold='0-46:5'
elif (opt.version==38):
    granularity='0-10:4,11-22:8'
    noise='0-10:0.14,11-22:0.2'
    threshold='0-22:5'
elif (opt.version==35):
    granularity='0-17:4'
    noise='0-17:0.14'
    threshold='0-17:5'
elif (opt.version==37):
    granularity='0-17:4,18-26:4,27-38:8'
    noise='0-26:0.14,27-38:0.2'
    threshold='0-38:5'
elif (opt.version==39):
    granularity='0-8:4,9-20:8'
    noise='0-8:0.14,9-20:0.2'
    threshold='0-20:5'
else:
    granularity='0-51:4'
    noise='0-51:0.15'
    threshold='0-51:5'

# interCalibList=[3];

suffix='IC3'

for et in enlist :

    nevents=opt.nevts
    # if en>150: nevents=nevents/2
    
    myqueue=opt.lqueue
    if et>0 and et<60 : myqueue=opt.squeue
    # if et>0 and et<20 : myqueue=opt.squeue
    
    bval="BOFF"
    if opt.Bfield>0 : bval="BON" 
    
    outDir='%s/git_%s/version_%d/model_%d/rwcuf%.1f/rwcum%.1f/%s/%s'%(opt.out,opt.gittag,opt.version,opt.model,opt.rwcuf,opt.rwcum,opt.datatype,bval)
    outDir='%s/%s'%(outDir,label)
    if et>0 : outDir='%s/et_%d'%(outDir,et)
    # eosDir='%s/git%s/%s'%(opt.eos,opt.gittag,opt.datatype)
    if opt.eta>0 : outDir='%s/eta_%3.3f/'%(outDir,opt.eta)
    # outDir='%s/run_no_%i'%(outDir,opt.runno)
    if (opt.run>=0) : outDir='%s/run_%d/'%(outDir,opt.run)

    eosDirIn=opt.out
    outlog='%s/digitizer.log'%(outDir)
    g4log='digijob%s.log'%(suffix)
    os.system('mkdir -p %s'%outDir)

    # wrapper
    scriptFile = open('%s/runJob.sh'%(outDir), 'w')
    scriptFile.write('#!/bin/bash\n')
    scriptFile.write('source %s/g4env.sh\n'%(os.getcwd()))
    # scriptFile.write('cd %s\n'%(outDir))
    scriptFile.write('cp %s/g4steer.mac .\n'%(outDir))
    scriptFile.write('PFCalEE g4steer.mac %d %d %f %s %s %s %s | tee g4.log\n'%(opt.version,opt.model,opt.eta,wthick,wcuthick,pbthick,droplayers))
    # outTag='%s_version%d_model%d_%s_runno_%i'%(label,opt.version,opt.model,bval,opt.runno)
    outTag='%s_version%d_model%d_%s_rwcuf_%.1f_rwcum_%.1f'%(label,opt.version,opt.model,bval,opt.rwcuf,opt.rwcum)
    if et>0 : outTag='%s_et%d'%(outTag,et)
    if opt.eta>0 : outTag='%s_eta%3.3f'%(outTag,opt.eta) 
    if (opt.run>=0) : outTag='%s_run%d'%(outTag,opt.run)
    scriptFile.write('scp -i /afs/cern.ch/user/t/tmudholk/.ssh/id_dsa PFcal.root tmudholk@cms02.phys.cmu.edu:/export/cmss2/tmudholk/HGCal/HGcal_%s.root\n'%(outTag))
    scriptFile.write('mv PFcal.root %s/HGcal_%s.root\n'%(eosDirIn,outTag))
    # scriptFile.write('cp %s/DigiPFcal.root %s/Digi_thr%f_%s.root\n'%(outDir,eosDirIn,opt.threshold,outTag))
    scriptFile.write('localdir=`pwd`\n')
    # scriptFile.write('%s/userlib/bin/digitizer 0 $localdir/HGcal_%s.root $localdir/ %s %s %s 3 0 | tee %s\n'%(os.getcwd(),outTag,granularity,noise,threshold,outlog))
    nPuVtx=0
    interCalib=3
    scriptFile.write('%s/userlib/bin/digitizer %d %s/HGcal_%s.root $localdir/ %s %s %s %d %d %s | tee %s\n'%(os.getcwd(),opt.nevts,eosDirIn,outTag,granularity,noise,threshold,interCalib,nPuVtx,INPATHPU,outlog))
    # scriptFile.write('rm %s/HGcal_%s.root\n'%(eosDirIn,outTag))
    scriptFile.write('echo "--Local directory is " $localdir >> g4.log\n')
    scriptFile.write('ls * >> g4.log\n')
    scriptFile.write('echo "--deleting core files: too heavy!!"\n')
    scriptFile.write('rm core.*\n')
    scriptFile.write('cp * %s/\n'%(outDir))
    scriptFile.write('scp -i /afs/cern.ch/user/t/tmudholk/.ssh/id_dsa %s/DigiPFcal.root tmudholk@cms02.phys.cmu.edu:/export/cmss2/tmudholk/HGCal/Digi_%s.root\n'%(outDir,outTag))
    scriptFile.write('cp %s/DigiPFcal.root %s/Digi_%s.root\n'%(outDir,eosDirIn,outTag))
    # scriptFile.write('rm %s/DigiPFcal.root\n'%(outDir))
    scriptFile.write('echo "All done"\n')
    scriptFile.close()
    
    # write geant 4 macro
    g4Macro = open('%s/g4steer.mac'%(outDir), 'w')
    g4Macro.write('/control/verbose 0\n')
    g4Macro.write('/control/saveHistory\n')
    g4Macro.write('/run/verbose 0\n')
    g4Macro.write('/event/verbose 0\n')
    g4Macro.write('/tracking/verbose 0\n')
    g4Macro.write('/N03/det/setField %1.1f T\n'%opt.Bfield)
    g4Macro.write('/N03/det/setModel %d\n'%opt.model)
    g4Macro.write('/random/setSeeds %d %d\n'%( random.uniform(0,100000), random.uniform(0,100000) ) )
    if opt.dogun :
        g4Macro.write('/generator/select particleGun\n')
        g4Macro.write('/gun/particle %s\n'%(opt.datatype))
        en=et*math.cosh(opt.eta)
        g4Macro.write('/gun/energy %f GeV\n'%(en))
        # g4Macro.write('/gun/direction %f %f %f\n'%(math.cos(math.pi*opt.phi)*math.sin(opt.alpha),math.sin(math.pi*opt.phi)*math.sin(opt.alpha),math.cos(opt.alpha)))
        # g4Macro.write('/gun/direction %f %f %f\n'%(random.uniform(0,1000)/100.-5.,math.sin(opt.alpha),math.cos(opt.alpha)))
    else :
        g4Macro.write('/generator/select hepmcAscii\n')
        g4Macro.write('/generator/hepmcAscii/open %s\n'%(opt.datafile))
        g4Macro.write('/generator/hepmcAscii/verbose 0\n')
    g4Macro.write('/run/beamOn %d\n'%(nevents))
    g4Macro.close()
    
    # submit
    os.system('chmod u+rwx %s/runJob.sh'%outDir)
    if opt.nosubmit : os.system('echo bsub -q %s %s/runJob.sh'%(myqueue,outDir)) 
    else: os.system("bsub -q %s \'%s/runJob.sh\'"%(myqueue,outDir))
