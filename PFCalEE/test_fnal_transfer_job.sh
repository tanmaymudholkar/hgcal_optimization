#!/bin/bash

xrdcp /afs/cern.ch/user/t/tmudholk/working_directory/test_fnal_transfer_file.dat root://cmseos.fnal.gov//store/user/tmudholk/test_fnal_transfer_file.dat

# kinit tmudholk@FNAL.GOV -k -t /afs/cern.ch/user/t/tmudholk/tanmay_kerberos.keytab
# rsync -ave ssh /afs/cern.ch/user/t/tmudholk/public/research/PFCal/PFCalEE/test_fnal_transfer_file.dat tmudholk@cmslpc-sl6.fnal.gov:/eos/uscms/store/user/tmudholk/test_fnal_transfer_file.dat
# kdestroy
