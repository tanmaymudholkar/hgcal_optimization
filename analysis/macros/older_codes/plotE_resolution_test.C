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
#include "TGraphErrors.h"
#include "TVector.h"

#include "../../userlib/include/HGCSSEvent.hh"
#include "../../userlib/include/HGCSSInfo.hh"
#include "../../userlib/include/HGCSSRecoHit.hh"
#include "../../userlib/include/HGCSSSimHit.hh"
#include "../../userlib/include/HGCSSSamplingSection.hh"

Double_t resolutions_fit(Double_t *x,Double_t *par) {
  return (par[0] + par[1]/sqrt(x[0]));
}

void plotE_resolution_test(){//main
  ofstream outfile;
  outfile.open(Form("test.txt"));
  outfile << "testing... 1... 2... 3..." << std::endl;
  outfile.close();
}
