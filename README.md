# HGCstandalone

#1. Install ROOT:
   see: https://root.cern.ch/drupal/content/production-version-534
   Download ROOT source file and unpack it
      cd root
      ./configure
      make

#2. Download the Standalone code:
     git clone https://github.com/mengleisun/HGCstandalone.git

#3. Setup ROOT environment:
     source PATH_TO_YOUR_ROOT/root/bin/thisroot.sh
     the path to the root can be found by "which root"

#4. Compile the userlib
      cd PFCal/PFCalEE/
      cd userlib
      mkdir {obj, lib, bin}
      make dictionary
      make

#5. run the example:
      cd analysis
      cd macros
      root plotE.C++ 
