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

#include "../../userlib/include/HGCSSEvent.hh"
#include "../../userlib/include/HGCSSInfo.hh"
#include "../../userlib/include/HGCSSRecoHit.hh"
#include "../../userlib/include/HGCSSSimHit.hh"
#include "../../userlib/include/HGCSSSamplingSection.hh"

void droplayer_resolution(){
  const int layers=28;
  //double layer_array[layers];
  std::vector<Double_t> layer_array;
  for(unsigned ilayer=0; ilayer<layers; ilayer++)
      //layer_array[ilayer]=(double)ilayer;
      layer_array.push_back((double)ilayer);
  gSystem->Load("/export/home/tmudholk/research/HGCstandalone/userlib/lib/libPFCalEEuserlib.so");

  std::vector<double> resol_b; // vector with double components
  ifstream f_resol_b; // input file stream
  f_resol_b.open("resol_b_eta2pt1_droplayer.dat"); //fopen
  std::string line; // a char string 
  double sig_b;
  if (f_resol_b.is_open()) {
    while (getline(f_resol_b,line)) { //=0 if the last line, line++
      sig_b = std::atof(line.c_str()); // line.c_str(): a string read from the line, atof: convert string to float
      resol_b.push_back(sig_b); // fill the last compoent of vector "weights", with double value "weight" and extend the array by one (like a stack)
    }
  }

  TVectorD Tresol_b(resol_b.size(),&resol_b[0]);
  TVectorD Tlayer_array(layer_array.size(),&layer_array[0]);
  //TGraph *resol_b_layer=new TGraph(28,layer_array,resol_b);
  TGraph *resol_b_layer=new TGraph(Tlayer_array,Tresol_b);
  TCanvas *myc=new TCanvas("droplayer_resol","droplayer_resol",800,600);
  myc->cd();
  resol_b_layer->Draw("AC*");
}
