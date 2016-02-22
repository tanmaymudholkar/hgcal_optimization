#include<string>
#include<iostream>
#include<fstream>
#include<sstream>
#include<cmath>
#include<stdlib.h>

#include "TFile.h"
#include "TTree.h"
#include "TF1.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TString.h"
#include "TChain.h"
#include "TSystem.h"
#include "TGraph.h"
#include "TVector.h"

#include "../../../userlib/include/HGCSSEvent.hh"
#include "../../../userlib/include/HGCSSInfo.hh"
#include "../../../userlib/include/HGCSSRecoHit.hh"
#include "../../../userlib/include/HGCSSSimHit.hh"
#include "../../../userlib/include/HGCSSSamplingSection.hh"

void droplayer_resolution(){
  TF1 *f1=new TF1("f1","0.00060+0.315/sqrt(x)",1,100);
  TF1 *f2=new TF1("f1","0.00401+0.234/sqrt(x)",1,100);
  TF1 *f3=new TF1("f1","0.00375+0.240/sqrt(x)",1,100);
  TF1 *f4=new TF1("f1","0.00040+0.317/sqrt(x)",1,100);
  
  TCanvas *myc=new TCanvas("myc","myc",1000,1000);
}
