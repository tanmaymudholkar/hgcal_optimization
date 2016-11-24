#include<string>
#include<set>
#include<iostream>
#include<fstream>
#include<sstream>
#include <boost/algorithm/string.hpp>

#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TNtuple.h"
#include "TH2D.h"
#include "TH2Poly.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TRandom3.h"
#include "Math/Vector4D.h"

#include "fastjet/ClusterSequence.hh"

#include "HGCSSEvent.hh"
#include "HGCSSInfo.hh"
#include "HGCSSSimHit.hh"
#include "HGCSSRecoHit.hh"
#include "HGCSSRecoJet.hh"
#include "HGCSSCalibration.hh"
#include "HGCSSDigitisation.hh"
#include "HGCSSDetector.hh"
#include "HGCSSGeometryConversion.hh"

using namespace fastjet;

template <class T>
void extractParameterFromStr(std::string aStr,T & vec){ 
  if (aStr == "") return;
  std::vector<std::string> layVec;
  boost::split( layVec, aStr, boost::is_any_of(","));

  for (unsigned iE(0); iE<layVec.size(); ++iE){//loop on elements
    std::vector<std::string> lPair;
    boost::split( lPair, layVec[iE], boost::is_any_of(":"));
    if (lPair.size() != 2) {
      std::cout << " -- Wrong string for parameter given as input:" << layVec[iE] << " Try again, expecting exactly one symbol \":\" between two \",\" ..." << std::endl;
      exit(1);
    }
    std::vector<std::string> lLay;
    boost::split( lLay, lPair[0], boost::is_any_of("-"));
    if (lLay.size() > 2) {
      std::cout << " -- Wrong string for granularities given as input:" << lPair[0] << " Try again, expecting at most one symbol \"-\"." << std::endl;
      exit(1);
    }
    unsigned beginIdx =  atoi(lLay[0].c_str());
    unsigned endIdx = lLay.size() == 1 ? beginIdx :  atoi(lLay[1].c_str());
    for (unsigned iL(beginIdx); iL<endIdx+1; ++iL){
      if (iL < vec.size())
	std::istringstream(lPair[1])>>vec[iL];
      else {
	std::cout << " -- WARNING! Input parameter has more layers: " << endIdx << " than detector : " 
		  << vec.size()
		  << ". Ignoring additional layer #" << iL << "... PLEASE CHECK SETTINGS ARE CORRECT FOR EXISTING LAYERS!!"
		  << std::endl;
      }
    }
  }//loop on elements
}

/*
void processHist(const unsigned iL,
		 std::map<unsigned,MergeCells> & histE,
		 HGCSSDigitisation & myDigitiser,
		 TH1F* & p_noise,
		 const TH2Poly* histZ,
		 const double & meanZpos,
		 const bool isTBsetup,
		 const HGCSSSubDetector & subdet,
		 const std::vector<unsigned> & pThreshInADC,
		 const bool pSaveDigis,
		 HGCSSRecoHitVec & lDigiHits,
		 HGCSSRecoHitVec & lRecoHits,
		 const bool pMakeJets,
		 std::vector<PseudoJet> & lParticles
		 ){

  bool doSaturation=true;

  DetectorEnum adet = subdet.type;
  bool isScint = subdet.isScint;
  bool isSi = subdet.isSi;

  //double rLim = subdet.radiusLim;

  TIter next(histE->GetBins());
  TObject *obj=0; 
  TH2PolyBin *polyBin = 0;
  
  while ((obj=next())){
    polyBin=(TH2PolyBin*)obj;
    double x = (polyBin->GetXMax()+polyBin->GetXMin())/2.;
    double y = (polyBin->GetYMax()+polyBin->GetYMin())/2.;
    //double radius = sqrt(pow(x,2)+pow(y,2));
    
    double digiE = 0;
    double simE = polyBin->GetContent();
    //double time = 0;
    //if (simE>0) time = histTime->GetBinContent(iX,iY)/simE;
    
    //fill vector with neighbours and calculate cross-talk
    double xtalkE = simE;
    if (isScint){
      std::vector<double> simEvec;
      simEvec.push_back(simE);
      double side = polyBin->GetXMax()-polyBin->GetXMin();
      simEvec.push_back(histE->GetBinContent(histE->FindBin(x-side,y)));
      simEvec.push_back(histE->GetBinContent(histE->FindBin(x+side,y)));
      simEvec.push_back(histE->GetBinContent(histE->FindBin(x,y-side)));
      simEvec.push_back(histE->GetBinContent(histE->FindBin(x,y+side)));
      xtalkE = myDigitiser.ipXtalk(simEvec);
    }
      
    //bool passTime = myDigitiser.passTimeCut(adet,time);
    //if (!passTime) continue;
      
    double posz = meanZpos;
    //for noise only hits
    //if (simE>0) posz = histZ->GetBinContent(iX,iY)/simE;
    //else posz = meanZpos;
    
    
    //if (fabs(x) > 500 || fabs(y)>500) std::cout << " x=" << x << ", y=" << y << std::endl;
    
    
    //correct for particle angle in conversion to MIP
    //not necessary, if not done for aborber thickness either
    double simEcor = xtalkE;//isTBsetup ? xtalkE : myDigitiser.mipCor(xtalkE,x,y,posz);
    digiE = simEcor;
    
    if (isScint && simEcor>0 && doSaturation) {
      digiE = myDigitiser.digiE(simEcor);
    }
    myDigitiser.addNoise(digiE,iL,p_noise);
    
    double noiseFrac = 1.0;
    if (simEcor>0) noiseFrac = (digiE-simEcor)/simEcor;
    
    //for silicon-based Calo
    unsigned adc = 0;
    if (isSi){
      adc = myDigitiser.adcConverter(digiE,adet);
      digiE = myDigitiser.adcToMIP(adc,adet);
    }
    bool aboveThresh = //digiE > 0.5;
      (isSi && adc > pThreshInADC[iL]) ||
      (isScint && digiE > pThreshInADC[iL]*myDigitiser.adcToMIP(1,adet,false));
    //histE->SetBinContent(iX,iY,digiE);
    if ((!pSaveDigis && aboveThresh) ||
	pSaveDigis)
      {//save hits
	//double calibE = myDigitiser.MIPtoGeV(subdet,digiE);
	HGCSSRecoHit lRecHit;
	lRecHit.layer(iL);
	lRecHit.energy(digiE);
	lRecHit.adcCounts(adc);
	lRecHit.x(x);
	lRecHit.y(y);
	lRecHit.z(posz);
	lRecHit.noiseFraction(noiseFrac);
	//unsigned x_cell = static_cast<unsigned>(fabs(x)/(cellSize*granularity[iL]));
	//unsigned y_cell = static_cast<unsigned>(fabs(y)/(cellSize*granularity[iL]));
	//lRecHit.encodeCellId(x>0,y>0,x_cell,y_cell,granularity[iL]);
	
	if (pSaveDigis) lDigiHits.push_back(lRecHit);
	
	lRecoHits.push_back(lRecHit);
	
	if (pMakeJets){
	  if (posz>0) lParticles.push_back( PseudoJet(lRecHit.px(),lRecHit.py(),lRecHit.pz(),lRecHit.E()));
	}
	
      }//save hits
  }//loop on bins
     
}//processHist
*/

void processHist(const unsigned iL,
		 std::map<unsigned,MergeCells> & histE,
		 std::map<int,std::pair<double,double> > & geom,
		 HGCSSDigitisation & myDigitiser,
		 TH1F* & p_noise,
		 //const TH2Poly* histZ,
		 const double & meanZpos,
		 const bool isTBsetup,
		 const HGCSSSubDetector & subdet,
		 const std::vector<unsigned> & pThreshInADC,
		 const bool pSaveDigis,
		 HGCSSRecoHitVec & lDigiHits,
		 HGCSSRecoHitVec & lRecoHits,
		 const bool pMakeJets,
		 std::vector<PseudoJet> & lParticles
		 ){

  bool doSaturation=false;

  DetectorEnum adet = subdet.type;
  bool isScint = subdet.isScint;
  bool isSi = subdet.isSi;

  //double rLim = subdet.radiusLim;
  for (unsigned iB(0); iB<histE.size();++iB){//loop on bins
    std::pair<double,double> xy = geom[iB+1];
    
    double digiE = 0;
    double simE = histE[iB].energy;
    //double time = 0;
    //if (simE>0) time = histE[iB].time/simE;
    
    //fill vector with neighbours and calculate cross-talk
    double xtalkE = simE;
    //CAMM @TODO
    /*if (isScint){
      std::vector<double> simEvec;
      simEvec.push_back(simE);
      double side = polyBin->GetXMax()-polyBin->GetXMin();
      simEvec.push_back(histE->GetBinContent(histE->FindBin(x-side,y)));
      simEvec.push_back(histE->GetBinContent(histE->FindBin(x+side,y)));
      simEvec.push_back(histE->GetBinContent(histE->FindBin(x,y-side)));
      simEvec.push_back(histE->GetBinContent(histE->FindBin(x,y+side)));
      xtalkE = myDigitiser.ipXtalk(simEvec);
      }*/
      
    //bool passTime = myDigitiser.passTimeCut(adet,time);
    //if (!passTime) continue;
      
    double posz = meanZpos;
    //for noise only hits
    //if (simE>0) posz = histZ->GetBinContent(iX,iY)/simE;
    //else posz = meanZpos;
    
    
    //if (fabs(x) > 500 || fabs(y)>500) std::cout << " x=" << x << ", y=" << y << std::endl;
    
    
    //correct for particle angle in conversion to MIP
    //not necessary, if not done for aborber thickness either
    double simEcor = xtalkE;//isTBsetup ? xtalkE : myDigitiser.mipCor(xtalkE,x,y,posz);
    digiE = simEcor;
    
    if (isScint && simEcor>0 && doSaturation) {
      digiE = myDigitiser.digiE(simEcor);
    }
    myDigitiser.addNoise(digiE,iL,p_noise);
    
    double noiseFrac = 1.0;
    if (simEcor>0) noiseFrac = (digiE-simEcor)/simEcor;
    
    //for silicon-based Calo
    unsigned adc = 0;
    if (isSi){
      adc = myDigitiser.adcConverter(digiE,adet);
      digiE = myDigitiser.adcToMIP(adc,adet);
    }
    bool aboveThresh = //digiE > 0.5;
      (isSi && adc > pThreshInADC[iL]) ||
      (isScint && digiE > pThreshInADC[iL]*myDigitiser.adcToMIP(1,adet,false));
    //histE->SetBinContent(iX,iY,digiE);
    if ((!pSaveDigis && aboveThresh) ||
	pSaveDigis)
      {//save hits
	//double calibE = myDigitiser.MIPtoGeV(subdet,digiE);
	HGCSSRecoHit lRecHit;
	lRecHit.layer(iL);
	lRecHit.energy(digiE);
	lRecHit.adcCounts(adc);
	lRecHit.x(xy.first);
	lRecHit.y(xy.second);
	lRecHit.z(posz);
	lRecHit.noiseFraction(noiseFrac);
	//unsigned x_cell = static_cast<unsigned>(fabs(x)/(cellSize*granularity[iL]));
	//unsigned y_cell = static_cast<unsigned>(fabs(y)/(cellSize*granularity[iL]));
	//lRecHit.encodeCellId(x>0,y>0,x_cell,y_cell,granularity[iL]);
	
	if (pSaveDigis) lDigiHits.push_back(lRecHit);
	
	lRecoHits.push_back(lRecHit);
	
	if (pMakeJets){
	  if (posz>0) lParticles.push_back( PseudoJet(lRecHit.px(),lRecHit.py(),lRecHit.pz(),lRecHit.E()));
	}
	
      }//save hits
  }//loop on bins
     
}//processHist


int main(int argc, char** argv){//main  

  /////////////////////////////////////////////////////////////
  //parameters
  /////////////////////////////////////////////////////////////
  const unsigned nReqA = 11;
  const unsigned nPar = static_cast<unsigned>(argc);
  if (nPar < nReqA-1) {
    std::cout << " Usage: "
              << argv[0] << " <nEvts to process (0=all)>"<< std::endl
              << "<full path to input file>"<< std::endl
              << "<full path to output file>"<< std::endl
              << "<granularities \"layer_i-layer_j:factor,layer:factor,...\">"<< std::endl
              << "<noise (in Mips) \"layer_i-layer_j:factor,layer:factor,...\">"<< std::endl
              << "<threshold (in ADC counts) \"layer_i-layer_j:factor,layer:factor,...\">"<< std::endl
              << "<intercalib factor in %>" << std::endl
	      << "<Number of si layers for TB setups>" << std::endl
              << "<number of PU>" << std::endl
              << "<full path to MinBias file if nPU!=0>" << std::endl
              << std::endl
              << "<optional: randomSeed (default=0)> "  << std::endl
              << "<optional: debug (default=0)>" << std::endl
              << "<optional: save sim hits (default=0)> " << std::endl
              << "<optional: save digi hits (default=0)> " << std::endl
              << "<optional: make jets (default=0)> " << std::endl
              << std::endl;
    return 1;
  }

  const unsigned pNevts = atoi(argv[1]);
  std::string inFilePath = argv[2];
  std::string outFilePath = argv[3];
  std::string granulStr = argv[4];
  std::string noiseStr = argv[5];
  std::string threshStr = argv[6];
  const unsigned interCalib = atoi(argv[7]);
  const unsigned nSiLayers = atoi(argv[8]);
  const unsigned nPU = atoi(argv[9]);
  std::string puPath;
  if (nPar > nReqA-1) puPath = argv[10];
  if (nPU>0 && puPath.size()==0) {
    std::cout << " -- Error! Missing full path to minbias file. Exiting." << std::endl;
    return 1;
  }

  std::cout << " ----------------------------------------" << std::endl
            << " -- Input parameters: " << std::endl
            << " -- Input file path: " << inFilePath << std::endl
            << " -- Output file path: " << outFilePath << std::endl
            << " -- Processing " ;
  if (pNevts>0) std::cout << pNevts;
  else std::cout << "all";
  std::cout << " events." << std::endl
            << " -- Granularities: " << granulStr << std::endl
            << " -- noise: " << noiseStr << std::endl
            << " -- thresholds: " << threshStr << std::endl
            << " -- intercalibration factor (in %): " << interCalib << std::endl
	    << " -- number of Si layers: " << nSiLayers << std::endl
            << " -- number of PU: " << nPU << std::endl
	    << " -- pu file path: " << puPath << std::endl
    ;

	    
  //std::string pModel = "model2";
  unsigned debug = 0;
  unsigned pSeed = 0;
  bool pSaveDigis = 0;
  bool pSaveSims = 0;
  bool pMakeJets = false;
  //if (nPar > nReqA-1) pModel = argv[nReqA];
  if (nPar > nReqA) std::istringstream(argv[nReqA])>>pSeed;
  if (nPar > nReqA+1) {
    debug = atoi(argv[nReqA+1]);
    std::cout << " -- DEBUG output is set to " << debug << std::endl;
  }
  if (nPar > nReqA+2) std::istringstream(argv[nReqA+2])>>pSaveDigis;
  if (nPar > nReqA+3) std::istringstream(argv[nReqA+3])>>pSaveSims;
  if (nPar > nReqA+4) std::istringstream(argv[nReqA+4])>>pMakeJets;
  
  //try to get model automatically
  //if (inFilePath.find("model0")!=inFilePath.npos) pModel = "model0";
  //else if (inFilePath.find("model1")!=inFilePath.npos) pModel = "model1";
  //else if (inFilePath.find("model2")!=inFilePath.npos) pModel = "model2";
  //else if (inFilePath.find("model3")!=inFilePath.npos) pModel = "model3";


  //std::cout << " -- Model is set to : " << pModel << std::endl;
  std::cout<< " -- Random seed will be set to : " << pSeed << std::endl;
  if (pSaveDigis) std::cout << " -- DigiHits are saved." << std::endl;
  if (pSaveSims) std::cout << " -- SimHits are saved." << std::endl;
  if (pMakeJets) std::cout << " -- Making jets." << std::endl;
  std::cout << " ----------------------------------------" << std::endl;
  
  //////////////////////////////////////////////////////////
  //// Hardcoded config ////////////////////////////////////
  //////////////////////////////////////////////////////////
  //for HGCAL, true means only 12 FHCAL layers considered (24 are simulated)
  bool concept = true;

  // choose a jet definition
  double R = 0.5;
  JetDefinition jet_def(antikt_algorithm, R);

  //////////////////////////////////////////////////////////
  //// End Hardcoded config ////////////////////////////////////
  //////////////////////////////////////////////////////////


  /////////////////////////////////////////////////////////////
  //input
  /////////////////////////////////////////////////////////////

  std::string inputStr = inFilePath ;
  TFile *inputFile = TFile::Open(inputStr.c_str());

  std::cout << "Opening input file " << inputStr << std::flush;

  if (!inputFile) {
    std::cout << " -- Error, input file " << inputStr << " cannot be opened. Exiting..." << std::endl;
    return 1;
  }

  std::cout << "...Done." << std::endl <<  "Getting HGCSSTree" << std::flush;

  TTree *inputTree = (TTree*)inputFile->Get("HGCSSTree");
  if (!inputTree){
    std::cout << " -- Error, tree HGCSSTree  cannot be opened. Exiting..." << std::endl;
    return 1;
  }

  std::cout << "...Done." <<  std::endl;

  /////////////////////////////////////////////////////////////
  //  //input PU hits
  /////////////////////////////////////////////////////////////

  bool signalIsPu = false;

  TChain *puTree = new TChain("HGCSSTree");
  unsigned nPuVtx = 0;
  unsigned nPuEvts = 0;
  std::vector<HGCSSSimHit> * puhitvec = 0;
  if(nPU!=0){
     ofstream myscript;
     myscript.open("eosls.sh");
     myscript<<"#!/bin/bash" << std::endl;
     myscript<<"/afs/cern.ch/project/eos/installation/0.3.15/bin/eos.select ls " << puPath << std::endl; 
     myscript.close();
     FILE *script = popen("bash eosls.sh", "r");
     char eoslsName[100];
     while(fgets(eoslsName, 100, script)) {
        std::ostringstream puInput;
        std::string temp = std::string(eoslsName).substr(0,strlen(eoslsName)-1);
        if(temp.find("Digi")!=temp.npos)continue;
        puInput << puPath << temp;
	if (puInput.str()==inFilePath) {
	  std::cout << " -- duplicate: file already used as signal. Removing." << std::endl; 
	  signalIsPu = true;
	  continue;
	}
        puTree->AddFile(puInput.str().c_str());
        std::cout << "Adding MinBias file:" << puInput.str().c_str() << std::endl;
     }
     pclose(script);  
     // system("rm ./eosls.sh");

    puTree->SetBranchAddress("HGCSSSimHitVec",&puhitvec);

    nPuEvts = puTree->GetEntries();
    std::cout << "- Number of PU events available: " << nPuEvts  << std::endl;
  }
  /////////////////////////////////////////////////////////////
  //input tree
  /////////////////////////////////////////////////////////////

  HGCSSInfo * info=(HGCSSInfo*)inputFile->Get("Info");

  assert(info);

  const double calorSizeXY = info->calorSizeXY();
  const double cellSize = info->cellSize();
  const unsigned versionNumber = info->version();
  const unsigned model = info->model();
  
  bool isTBsetup = (model != 2);
  if (isTBsetup) std::cout << " -- Number of Si layers: " << nSiLayers << std::endl;
  else std::cout << " -- Number of Si layers ignored: hardcoded as a function of radius in HGCSSGeometryConversion class." << std::endl;


  HGCSSEvent * event=0;
  std::vector<HGCSSSimHit> * hitvec = 0;

  inputTree->SetBranchAddress("HGCSSEvent",&event);
  inputTree->SetBranchAddress("HGCSSSimHitVec",&hitvec);
    
  //initialise detector
  HGCSSDetector & myDetector = theDetector();
  bool isCaliceHcal = versionNumber==23;//inFilePath.find("version23")!=inFilePath.npos || inFilePath.find("version_23")!=inFilePath.npos;

  //unsigned versionNumber = 0;
  //unsigned nchar = 0;
  //if (inFilePath.substr(inFilePath.find("version")+8,1)=="_") nchar = 1;
  //else nchar = 2;
  //std::string lvers = inFilePath.substr(inFilePath.find("version")+7,nchar);
  //std::istringstream(lvers)>>versionNumber;
  
  std::cout << " -- Calor size XY = " << calorSizeXY
	    << ", version number = " << versionNumber 
	    << ", model = " << model
	    << ", cellSize = " << cellSize
	    << std::endl;

  bool bypassR = false;
  if (isTBsetup) bypassR = true;
  myDetector.buildDetector(versionNumber,concept,isCaliceHcal,bypassR);

  //initialise calibration class
  HGCSSCalibration mycalib(inFilePath,bypassR,nSiLayers);

  const unsigned nLayers = myDetector.nLayers();

  HGCSSGeometryConversion geomConv(model,cellSize,bypassR,nSiLayers);
  geomConv.setXYwidth(calorSizeXY);
  //const double xWidth = geomConv.getXYwidth();
  geomConv.initialiseHoneyComb(calorSizeXY,cellSize);
  //square map for BHCAL
  geomConv.initialiseSquareMap(calorSizeXY/2.,10.);

  HGCSSDigitisation myDigitiser;

  std::vector<unsigned> granularity;
  granularity.resize(nLayers,1);
  std::vector<double> pNoiseInMips;
  pNoiseInMips.resize(nLayers,0.12);
  std::vector<unsigned> pThreshInADC;
  pThreshInADC.resize(nLayers,5);

  extractParameterFromStr<std::vector<unsigned> >(granulStr,granularity);
  extractParameterFromStr<std::vector<double> >(noiseStr,pNoiseInMips);
  extractParameterFromStr<std::vector<unsigned> >(threshStr,pThreshInADC);

  //unsigned nbCells = 0;

  std::cout << " -- Granularities and noise are setup like this:" << std::endl;
  for (unsigned iL(0); iL<nLayers; ++iL){
    std::cout << "Layer " ;
    if (iL<10) std::cout << " ";
    std::cout << iL << " : " << granularity[iL] << ", " << pNoiseInMips[iL] << " mips, " << pThreshInADC[iL] << " adc - ";
    if (iL%5==4) std::cout << std::endl;
    myDigitiser.setNoise(iL,pNoiseInMips[iL]);
    //nbCells += N_CELLS_XY_MAX/(granularity[iL]*granularity[iL]);
  }
  std::cout << std::endl;     
  //std::cout << " -- Total number of cells = " << nbCells << std::endl;

  geomConv.setGranularity(granularity);
  geomConv.initialiseHistos();

  TRandom3 *lRndm = new TRandom3();
  lRndm->SetSeed(pSeed);
  myDigitiser.setRandomSeed(pSeed);

  std::cout << " -- Random3 seed = " << lRndm->GetSeed() << std::endl
	    << " ----------------------------------------" << std::endl;

  myDigitiser.setIntercalibrationFactor(interCalib);

  /////////////////////////////////////////////////////////////
  //output
  /////////////////////////////////////////////////////////////

  std::ostringstream outputStr;
  outputStr << outFilePath << "/DigiPFcal" ;
  if (pSaveDigis)  outputStr << "_withDigiHits";
  if (pSaveSims)  outputStr << "_withSimHits";
  outputStr << ".root";
  
  TFile *outputFile = TFile::Open(outputStr.str().c_str(),"RECREATE");

  if (!outputFile) {
    std::cout << " -- Error, output file " << outputStr.str() << " cannot be opened. Exiting..." << std::endl;
    return 1;
  }
  else {
    std::cout << " -- File will be saved as " << outputStr.str() << std::endl;
  }

  HGCSSInfo *lInfo = new HGCSSInfo();
  lInfo->calorSizeXY(calorSizeXY);
  lInfo->cellSize(cellSize);
  lInfo->version(versionNumber);
  lInfo->model(model);

  TTree *outputTree = new TTree("RecoTree","HGC Standalone simulation reco tree");
  HGCSSSimHitVec lSimHits;
  HGCSSRecoHitVec lDigiHits;
  HGCSSRecoHitVec lRecoHits;
  HGCSSRecoJetVec lCaloJets;
  unsigned maxSimHits = 0;
  unsigned maxRecHits = 0;
  unsigned maxRecJets = 0;
  HGCSSEvent lEvent;
  outputTree->Branch("HGCSSEvent",&lEvent);
  if (nPU!=0) outputTree->Branch("nPuVtx",&nPuVtx);
  if (pSaveSims) outputTree->Branch("HGCSSSimHitVec","std::vector<HGCSSSimHit>",&lSimHits);
  if (pSaveDigis) outputTree->Branch("HGCSSDigiHitVec","std::vector<HGCSSRecoHit>",&lDigiHits);
  outputTree->Branch("HGCSSRecoHitVec","std::vector<HGCSSRecoHit>",&lRecoHits);
  if (pMakeJets) outputTree->Branch("HGCSSRecoJetVec","std::vector<HGCSSRecoJet>",&lCaloJets);
  TH1F * p_noise = new TH1F("noiseCheck",";noise (MIPs)",100,-5,5);


  /////////////////////////////////////////////////////////////
  //Loop on events
  /////////////////////////////////////////////////////////////

  const unsigned nEvts = (pNevts > inputTree->GetEntries() || pNevts==0) ? static_cast<unsigned>(inputTree->GetEntries()) : pNevts;

  std::cout << "- Processing = " << nEvts  << " events out of " << inputTree->GetEntries() << std::endl;

  std::vector<PseudoJet> lParticles;

  for (unsigned ievt(0); ievt<nEvts; ++ievt){//loop on entries

    inputTree->GetEntry(ievt);
    lEvent.eventNumber(event->eventNumber());
    lEvent.vtx_x(event->vtx_x());
    lEvent.vtx_y(event->vtx_y());
    lEvent.vtx_z(event->vtx_z());
    //unsigned layer = volNb;
    
    mycalib.setVertex(lEvent.vtx_x(),lEvent.vtx_y(),lEvent.vtx_z());

    if (debug>0) {
      std::cout << " **DEBUG** Processing evt " << ievt << std::endl;
    }
    else if (ievt%50 == 0) std::cout << "... Processing event: " << ievt << std::endl;
    

    unsigned prevLayer = 10000;
    DetectorEnum type = DetectorEnum::FECAL;
    unsigned subdetLayer=0;
    bool isScint = false;
    for (unsigned iH(0); iH<(*hitvec).size(); ++iH){//loop on hits
      HGCSSSimHit lHit = (*hitvec)[iH];

      //do not save hits with 0 energy...
      if (lHit.energy()>0 && pSaveSims) lSimHits.push_back(lHit);
      unsigned layer = lHit.layer();
      if (layer != prevLayer){
	const HGCSSSubDetector & subdet = myDetector.subDetectorByLayer(layer);
	isScint = subdet.isScint;
	type = subdet.type;
	subdetLayer = layer-subdet.layerIdMin;
	prevLayer = layer;
	if (debug > 1) std::cout << " - layer " << layer << " " << subdet.name << " " << subdetLayer << std::endl;
      }
      std::pair<double,double> xy = lHit.get_xy(isScint,geomConv);
      double posx = xy.first;//lHit.get_x(cellSize);
      double posy = xy.second;//lHit.get_y(cellSize);
      double radius = sqrt(pow(posx,2)+pow(posy,2));
      double posz = lHit.get_z();
      double energy = lHit.energy()*mycalib.MeVToMip(layer,radius); // if (energy > 0) std::cout << "sim energy = "<<lHit.energy()<<", reco energy = "<<energy<<std::endl;
      double realtime = mycalib.correctTime(lHit.time(),posx,posy,posz);
      bool passTime = myDigitiser.passTimeCut(type,realtime);
      if (!passTime) continue;
      if (energy>0 && 
	  lHit.silayer() < geomConv.getNumberOfSiLayers(type,radius) 
	  ){
	if (debug > 1) std::cout << " hit " << iH 
				 << " lay " << layer  
				 << " x " << posx 
				 << " y " << posy
				 << " z " << posz
				 << " t " << lHit.time() << " " << realtime
				 << std::endl;
	//geomConv.fill(type,subdetLayer,energy,realtime,posx,posy,posz);
	geomConv.fill(layer,energy,realtime,lHit.cellid(),posz);
      }

    }//loop on input simhits

    if(nPU!=0){
      //get PU events
      //std::vector<unsigned> ipuevt;

      //get poisson <140>
      nPuVtx = lRndm->Poisson(nPU);
      //ipuevt.resize(nPuVtx,1);
      if (signalIsPu) nPuVtx -= 1;
      std::cout << " -- Adding " << nPuVtx << " events to signal event: " << ievt << std::endl;
      std::set<unsigned> lidxSet;
      for (unsigned iV(0); iV<nPuVtx; ++iV){//loop on interactions
        unsigned ipuevt = 0;
	while (1){
	  ipuevt = lRndm->Integer(nPuEvts);
	  if (lidxSet.find(ipuevt)==lidxSet.end()){
	    lidxSet.insert(ipuevt);
	    break;
	  }
          else {
            std::cout << " -- Found duplicate ! Taking another shot." << std::endl;
	  }
	}
	//std::cout << " ---- adding evt " << ipuevt << std::endl;
        puTree->GetEntry(ipuevt);
        //lRecoHits.reserve(lRecoHits.size()+(*puhitvec).size());
        prevLayer = 10000;
        type = DetectorEnum::FECAL;
        subdetLayer=0;
	isScint = false;
        for (unsigned iH(0); iH<(*puhitvec).size(); ++iH){//loop on hits
          HGCSSSimHit lHit = (*puhitvec)[iH];
          if (lHit.energy()>0){
	    unsigned layer = lHit.layer();
	    if (layer != prevLayer){
	      const HGCSSSubDetector & subdet = myDetector.subDetectorByLayer(layer);
	      isScint = subdet.isScint;
	      type = subdet.type;
	      subdetLayer = layer-subdet.layerIdMin;
	      prevLayer = layer;
	      //std::cout << " - layer " << layer << " " << subdet.name << " " << subdetLayer << std::endl;
            }
	    std::pair<double,double> xy = lHit.get_xy(isScint,geomConv);
	    double posx = xy.first;//lHit.get_x(cellSize);
	    double posy = xy.second;//lHit.get_y(cellSize);
            double posz = lHit.get_z();
            double radius = sqrt(pow(posx,2)+pow(posy,2));
	    if (lHit.silayer() < geomConv.getNumberOfSiLayers(type,radius)){
	      double energy = lHit.energy()*mycalib.MeVToMip(layer,radius);
	      double realtime = mycalib.correctTime(lHit.time(),posx,posy,posz);
	      bool passTime = myDigitiser.passTimeCut(type,realtime);
              if (!passTime) continue;
	      
	      if (debug > 1) std::cout << " hit " << iH
				       << " lay " << layer
				       << " x " << posx
				       << " y " << posy
				       << " z " << posz
				       << " t " << lHit.time() << " " << realtime
				       << std::endl;
	      //geomConv.fill(type,subdetLayer,energy,realtime,posx,posy,posz);
	      geomConv.fill(layer,energy,realtime,lHit.cellid(),posz);
	    }
          }
        }//loop on hits
      }//loop on interactions
    }//add PU
    
    if (debug>0) {
      std::cout << " **DEBUG** simhits = " << (*hitvec).size() << " " << lSimHits.size() << std::endl;
    }

    //create hits, everywhere to have also pure noise
    //digitise
    //apply threshold
    //save
    unsigned nTotBins = 0;
    for (unsigned iL(0); iL<nLayers; ++iL){//loop on layers
      //TH2Poly *histE = geomConv.get2DHist(iL,"E");
      std::map<unsigned,MergeCells> & histE = geomConv.get2DHist(iL);
      //TH2Poly *histTime = geomConv.get2DHist(iL,"Time");
      //TH2Poly *histZ = geomConv.get2DHist(iL,"Z");
      const HGCSSSubDetector & subdet = myDetector.subDetectorByLayer(iL);
      bool isScint = subdet.isScint;
      //nTotBins += histE->GetNumberOfBins();
      unsigned nBins = isScint?geomConv.squareMap()->GetNumberOfBins() : geomConv.hexagonMap()->GetNumberOfBins();
      nTotBins += nBins;
      if (pSaveDigis) lDigiHits.reserve(nTotBins);
      
      double meanZpos = geomConv.getAverageZ(iL);
      
      //extend map to include all cells
      for (unsigned iB(1); iB<nBins+1;++iB){
	MergeCells tmpCell;
	tmpCell.energy = 0;
	tmpCell.time = 0;
	histE.insert(std::pair<unsigned,MergeCells>(iB,tmpCell));
      }

      //std::cout << iL << " " << meanZpos << std::endl;

      if (debug>0){
	std::cout << " -- Layer " << iL << " " << subdet.name << " z=" << meanZpos
		  << " bins = " << nBins << " histE entries = " << histE.size() << std::endl;
      }

      //cell-to-cell cross-talk for scintillator
      if (isScint){
	//2.5% per 30-mm edge
	myDigitiser.setIPCrossTalk(0.025*10.*granularity[iL]/30.);
      }
      else {
	myDigitiser.setIPCrossTalk(0);
      }

      //processHist(iL,histE,myDigitiser,p_noise,histZ,meanZpos,isTBsetup,subdet,pThreshInADC,pSaveDigis,lDigiHits,lRecoHits,pMakeJets,lParticles);
      processHist(iL,histE,isScint?geomConv.squareGeom:geomConv.hexaGeom,myDigitiser,p_noise,meanZpos,isTBsetup,subdet,pThreshInADC,pSaveDigis,lDigiHits,lRecoHits,pMakeJets,lParticles);
 
    }//loop on layers

    if (debug) {
      std::cout << " **DEBUG** sim-digi-reco hits = " << (*hitvec).size() << "-" << lDigiHits.size() << "-" << lRecoHits.size() << std::endl;
    }
    
    
    if (pMakeJets){//pMakeJets
      
      // run the clustering, extract the jets
      ClusterSequence cs(lParticles, jet_def);
      std::vector<PseudoJet> jets = sorted_by_pt(cs.inclusive_jets());
      
      // print the jets
      std::cout <<   "-- evt " << ievt << ": found " << jets.size() << " Jets." << std::endl;
      for (unsigned i = 0; i < jets.size(); i++) {
	const PseudoJet & lFastJet = jets[i];
	//TOFIX // inverted y and z...
	HGCSSRecoJet ljet(lFastJet.px(),
			  lFastJet.py(),
			  lFastJet.pz(),
			  lFastJet.E());
	if (lFastJet.has_constituents()) ljet.nConstituents(lFastJet.constituents().size());
	if (lFastJet.has_area()){
	  ljet.area(lFastJet.area());
	  ljet.area_error(lFastJet.area_error());
	}
	
	lCaloJets.push_back(ljet);
	std::cout << " -------- jet " << i << ": "
		  << lFastJet.E() << " " 
		  << lFastJet.perp() << " " 
		  << lFastJet.rap() << " " << lFastJet.phi() << " "
		  << lFastJet.constituents().size() << std::endl;
	// std::vector<PseudoJet> constituents = lFastJet.constituents();
	// for (unsigned j = 0; j < constituents.size(); j++) {
	//   std::cout << "    constituent " << j << "'s pt: " << constituents[j].perp()
	// 	      << std::endl;
	// }
      }
      
    }//pMakeJets
    
    outputTree->Fill();
    //reserve necessary space and clear vectors.
    if (lSimHits.size() > maxSimHits) {
      maxSimHits = 2*lSimHits.size();
      std::cout << " -- INFO: event " << ievt << " maxSimHits updated to " << maxSimHits << std::endl;
    }
    if (lRecoHits.size() > maxRecHits) {
      maxRecHits = 2*lRecoHits.size();
      std::cout << " -- INFO: event " << ievt << " maxRecHits updated to " << maxRecHits << std::endl;
    }
    if (lCaloJets.size() > maxRecJets) {
      maxRecJets = 2*lCaloJets.size();
      std::cout << " -- INFO: event " << ievt << " maxRecJets updated to " << maxRecJets << std::endl;
    }
    lSimHits.clear();
    lDigiHits.clear();
    lRecoHits.clear();
    lCaloJets.clear();
    geomConv.initialiseHistos();
    lParticles.clear();
    if (pSaveSims) lSimHits.reserve(maxSimHits);
    lRecoHits.reserve(maxRecHits);
    lParticles.reserve(maxRecHits);
    lCaloJets.reserve(maxRecJets);
    
  }//loop on entries

  outputFile->cd();
  outputFile->WriteObjectAny(lInfo,"HGCSSInfo","Info");
  outputTree->Write();
  p_noise->Write();
  outputFile->Close();

  return 0;

}//main
