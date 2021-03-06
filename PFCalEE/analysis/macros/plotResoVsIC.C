#include<string>
#include<iostream>
#include<fstream>
#include<sstream>
#include<iomanip>

#include "TFile.h"
#include "TTree.h"
#include "TNtuple.h"
#include "TH2F.h"
#include "TH1F.h"
#include "TLine.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TGraphErrors.h"
#include "TLegend.h"
#include "TF1.h"
#include "TString.h"
#include "TLatex.h"
#include "TGaxis.h"

#include "TDRStyle.h"

int plotResoVsIC(){

  SetTdrStyle();


  const unsigned nIC = 10;
  const unsigned ICval[nIC] = {0,1,2,3,4,5,10,15,20,50};


  std::ostringstream label;
  TFile *fcalib[nIC];
  
  TGraphErrors *constant = new TGraphErrors();
  constant->SetName("constant");
  constant->SetTitle(";intercalib. smearing");
  constant->SetMarkerStyle(20);
  constant->SetMarkerColor(1);
  constant->SetLineColor(1);
  TGraphErrors *constantSR7 =  new TGraphErrors();
  constantSR7->SetName("constantSR7");
  constantSR7->SetTitle(";intercalib. smearing");
  constantSR7->SetMarkerStyle(23);
  constantSR7->SetMarkerColor(2);
  constantSR7->SetLineColor(2);

  TGraphErrors *noise = (TGraphErrors *) constant->Clone("noise");
  TGraphErrors *sampling = (TGraphErrors *) constant->Clone("sampling");
  TGraphErrors *samplingSR7 = (TGraphErrors *) constantSR7->Clone("samplingSR7");

  TCanvas *mycReso = new TCanvas("mycReso","mycReso",1500,1000);
  mycReso->Divide(2,5);
  TCanvas *mycR = new TCanvas("mycR","Sampling",1500,1000);
  TCanvas *mycC = new TCanvas("mycC","Constant",1500,1000);
  TCanvas *mycN = new TCanvas("mycN","Noise",1500,1000);

  gStyle->SetOptFit(1111);
  gStyle->SetOptStat(0);

  gStyle->SetStatW(0.2);
  gStyle->SetStatH(0.5);

  TLatex lat;
  char buf[500];

  TGraphErrors *gr[nIC][2];
  double x0,y0;
  double x0_7,y0_7;

  for (unsigned ic(0);ic<nIC;++ic){//loop on intercalib
    label.str("");
    label << "PLOTS/CalibReso";
    label << "_vsE";
    label << "_IC" << ICval[ic];
    label << ".root";
    fcalib[ic] = TFile::Open(label.str().c_str());
    if (!fcalib[ic]) {
      std::cout << " -- failed to open file: " << label.str() << std::endl;
      continue;
    }
    else {
      std::cout << " -- file " << label.str() << " successfully opened." << std::endl;
    }
    fcalib[ic]->cd("SR2");
    gr[ic][0] = (TGraphErrors *)gDirectory->Get("resoRecoFit2eta21pu1");
    fcalib[ic]->cd("SR7");
    gr[ic][1] = (TGraphErrors *)gDirectory->Get("resoRecoFit7eta21pu1");


    TF1 *fit = gr[ic][0]->GetFunction("reso");
    TF1 *fit7 = gr[ic][1]->GetFunction("reso");
    mycReso->cd(ic+1);
    gr[ic][0]->Draw("APE");
    fit->SetLineColor(6);
    fit->Draw("same");
    lat.SetTextSize(0.1);
    sprintf(buf,"Single #gamma, #eta=2.1, 3#times3 cm^{2}");
    lat.DrawLatexNDC(0.2,0.8,buf);
    sprintf(buf,"ICsmear = %d %%",ICval[ic]);
    lat.DrawLatexNDC(0.2,0.7,buf);

    double cval = sqrt(pow(fit->GetParameter(1),2)-pow(y0,2));
    constant->SetPoint(ic,ICval[ic]/100.,cval);
    constant->SetPointError(ic,0,fit->GetParameter(1)*fit->GetParError(1)/cval);
    noise->SetPoint(ic,ICval[ic]/100.,fit->GetParameter(2));
    noise->SetPointError(ic,0,fit->GetParError(2));
    sampling->SetPoint(ic,ICval[ic]/100.,fit->GetParameter(0));
    sampling->SetPointError(ic,0,fit->GetParError(0));
    cval = sqrt(pow(fit7->GetParameter(1),2)-pow(y0_7,2));
    constantSR7->SetPoint(ic,ICval[ic]/100.,cval);
    constantSR7->SetPointError(ic,0,fit7->GetParameter(1)*fit7->GetParError(1)/cval);    
    //constantSR7->SetPoint(ic,ICval[ic]/100.,fit7->GetParameter(1));
    //constantSR7->SetPointError(ic,0,fit7->GetParError(1));
    samplingSR7->SetPoint(ic,ICval[ic]/100.,fit7->GetParameter(0));
    samplingSR7->SetPointError(ic,0,fit7->GetParError(0));

    if (ic==0) {
      constant->GetPoint(0,x0,y0);
      constantSR7->GetPoint(0,x0_7,y0_7);
      cval = sqrt(pow(fit->GetParameter(1),2)-pow(y0,2));
      constant->SetPoint(ic,ICval[ic]/100.,cval);
      constant->SetPointError(ic,0,fit->GetParameter(1)*fit->GetParError(1)/cval);
      cval = sqrt(pow(fit7->GetParameter(1),2)-pow(y0_7,2));
      constantSR7->SetPoint(ic,ICval[ic]/100.,cval);
      constantSR7->SetPointError(ic,0,fit7->GetParameter(1)*fit7->GetParError(1)/cval);
    }

  }

  mycReso->Update();
  mycReso->Print("PLOTS/ResolutionFitvsIC.pdf");
  
  TLegend *leg = new TLegend(0.6,0.6,0.8,0.8);
  leg->SetFillColor(10);
  leg->AddEntry(sampling,"3#times3 cm^{2}","P");
  leg->AddEntry(samplingSR7,"All detector","P");
  mycR->cd();
  gPad->SetGridy(1);
  sampling->GetYaxis()->SetTitle("Sampling term (GeV^{1/2})");
  sampling->SetMinimum(0.2);
  sampling->SetMaximum(0.3);
  sampling->Draw("APE");
  samplingSR7->Draw("PEsame");
  lat.SetTextSize(0.04);
  sprintf(buf,"Single #gamma, #eta=2.1");
  lat.DrawLatexNDC(0.2,0.87,buf);
  lat.DrawLatexNDC(0.01,0.01,"HGCAL G4 standalone");

  leg->Draw("same");
  mycR->Update();

  mycR->Print("PLOTS/SamplingvsIC.pdf");

  mycC->cd();
  gPad->SetLogx(1);
  gPad->SetGridy(1);
  gStyle->SetOptFit(0);
  //gStyle->SetStatH(0.1);
  //gStyle->SetStatW(0.2);

  constant->GetYaxis()->SetTitle("Constant from intercalib.");
  constant->SetMinimum(0);
  constant->SetMaximum(0.08);
  constant->Draw("APE");
  //constantSR7->Draw("PEsame");
  sprintf(buf,"Single #gamma, #eta=2.1");
  lat.DrawLatexNDC(0.2,0.87,buf);
  lat.DrawLatexNDC(0.01,0.01,"HGCAL G4 standalone");

  TF1 *BE = new TF1("BE","sqrt([0]*[0] + pow(x*1/sqrt([1]),2))",0,1);
  BE->SetParameters(0,30);
  //  BE->SetParLimits(0,1,1);
  BE->FixParameter(0,0);
  BE->SetLineColor(1);
  constant->Fit("BE","BI");

  lat.SetTextColor(1);
  //sprintf(buf,"c #propto c_{0} #oplus #frac{x}{#sqrt{n}}, n=%3.1f #pm %3.1f",BE->GetParameter(1),BE->GetParError(1));
  sprintf(buf,"c_{ic}=#frac{x}{#sqrt{n}}, n=%3.1f #pm %3.1f",BE->GetParameter(1),BE->GetParError(1));
  lat.DrawLatexNDC(0.2,0.77,"c=c_{0} #oplus c_{ic}");
  lat.DrawLatexNDC(0.2,0.67,buf);

  //BE->SetParameter(0,y0_7);
  //constantSR7->Fit("BE","BI","same");
  //BE->SetLineColor(2);
  //BE->Draw("same");
  //lat.SetTextColor(2);
  //sprintf(buf,"c #propto c_{0} #oplus #frac{x}{#sqrt{n}}, n=%3.1f #pm %3.1f",BE->GetParameter(1),BE->GetParError(1));
  //lat.DrawLatexNDC(0.2,0.6,buf);
  //lat.SetTextColor(1);

  //leg->Draw("same");
  mycC->Update();
  mycC->Print("PLOTS/ConstantvsIC.pdf");


  mycN->cd();
  noise->GetYaxis()->SetTitle("Noise term (GeV)");
  noise->Draw("APE");
  sprintf(buf,"Single #gamma, #eta=2.1, 3#times3 cm^{2}");
  lat.DrawLatexNDC(0.2,0.87,buf);
  lat.DrawLatexNDC(0.01,0.01,"HGCAL G4 standalone");
  mycN->Print("PLOTS/NoisevsIC.pdf");

  




  return 0;

}//main
