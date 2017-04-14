#include "HGCSSSimHit.hh"

#include <iomanip>
#include <cmath>
#include <stdlib.h>

HGCSSSimHit::HGCSSSimHit(const G4SiHit & aSiHit, 
			 const unsigned & asilayer, 
			 TH2Poly* map, 
			 const float ){
  energy_ = aSiHit.energy;
  //energy weighted time
  //PS: need to call calculateTime() after all hits 
  //have been added to have divided by totalE!!
  time_ = aSiHit.time*aSiHit.energy;
  zpos_ = aSiHit.hit_z;
  setLayer(aSiHit.layer,asilayer);
  bool radialMapToBeUsed = (aSiHit.layer >= ANNULARGEOMETRYFIRSTLAYER);

  //coordinates in mm
  //double z = aSiHit.hit_x;
  double x = aSiHit.hit_x;
  double y = aSiHit.hit_y;
  //cellid encoding:
  //map->Reset("");
  //map->Fill(x,y);
  //GetMaximumBin doesn't work :(

  assert(map);
  if (radialMapToBeUsed) {
    double r, phi;
    r = sqrt(x*x + y*y);
    phi = (ROOT::Math::XYZPoint(x, y, 0.).phi());
    if (phi < 0) phi += 2*TMath::Pi();
    cellid_ = map->FindBin(phi, r);
  }
  else cellid_ = map->FindBin(x, y);

  //for (int ix(1);ix<map->GetNumberOfBins()+1; ++ix){
  //if (map->GetBinContent(ix)!=0)
    //cellid_ = ix;
    //std::cout << ix << " " << map->GetBinContent(ix) << std::endl;
  //}

  /*
  TIter next(map->GetBins());
  TObject *obj=0; 
  TH2PolyBin *polyBin = 0;

  while ((obj=next())){
    polyBin=(TH2PolyBin*)obj;
    int id = polyBin->GetBinNumber();
    if (id==cellid_) break; 
  }

  std::cout << " - Sanity check: x,y = " << x << " " << y 
	    << " cellid=" << cellid_ 
	    << " polybin# " << polyBin->GetBinNumber() << " area " << polyBin->GetArea() << std::endl
	    << " x bin = " << polyBin->GetXMin() << "-" << polyBin->GetXMax() << std::endl
	    << " y bin = " << polyBin->GetYMin() << "-" << polyBin->GetYMax() << std::endl
	    << " middle = " << (polyBin->GetXMax()+polyBin->GetXMin())/2 
	    << " " << (polyBin->GetYMax()+polyBin->GetYMin())/2 
	    << std::endl;
*/


  //bool x_side = x>0 ? true : false;
  //bool y_side = y>0 ? true : false;
  //unsigned x_cell = static_cast<unsigned>(fabs(x)/(cellSize*getGranularity()));
  //unsigned y_cell = static_cast<unsigned>(fabs(y)/(cellSize*getGranularity()));

  //encodeCellId(x_side,y_side,x_cell,y_cell);

  nGammas_= 0;
  nElectrons_ = 0;
  nMuons_ = 0;
  nNeutrons_ = 0;
  nProtons_ = 0;
  nHadrons_ = 0;
  if(abs(aSiHit.pdgId)==22) nGammas_++;
  else if(abs(aSiHit.pdgId)==11) nElectrons_++;
  else if(abs(aSiHit.pdgId)==13) nMuons_++;
  else if(abs(aSiHit.pdgId)==2112) nNeutrons_++;
  else if(abs(aSiHit.pdgId)==2212) nProtons_++;
  else nHadrons_++;

  trackIDMainParent_ = aSiHit.parentId;
  energyMainParent_ = aSiHit.energy;

}

/*void HGCSSSimHit::encodeCellId(const bool x_side,const bool y_side,const unsigned x_cell,const unsigned y_cell){
  cellid_ = 
    x_side | (x_cell<<1) |
    (y_side<<16) | (y_cell<<17);

  // std::cout << " Cross-check of encoding: cellid=" << cellid_ << std::endl
  // 	    << " x_side " << x_side << " " << get_x_side() << std::endl
  // 	    << " y_side " << y_side << " " << get_y_side() << std::endl
  // 	    << " x_cell " << x_cell << " " << get_x_cell() << std::endl
  // 	    << " y_cell " << y_cell << " " << get_y_cell() << std::endl
  //   ;
  }*/

void HGCSSSimHit::Add(const G4SiHit & aSiHit){

  time_ = time_ + aSiHit.time*aSiHit.energy;
  //PS: need to call calculateTime() after all hits 
  //have been added to have divided by totalE!!

  if(abs(aSiHit.pdgId)==22) nGammas_++;
  else if(abs(aSiHit.pdgId)==11) nElectrons_++;
  else if(abs(aSiHit.pdgId)==13) nMuons_++;
  else if(abs(aSiHit.pdgId)==2112) nNeutrons_++;
  else if(abs(aSiHit.pdgId)==2212) nProtons_++;
  else nHadrons_++;

  energy_ += aSiHit.energy;
  if (aSiHit.energy > energyMainParent_){
    trackIDMainParent_ = aSiHit.parentId;
    energyMainParent_ = aSiHit.energy;
  }

}

/*double HGCSSSimHit::eta() const {
  double x = get_x();
  double y = get_y();
  double theta = acos(fabs(zpos_)/sqrt(zpos_*zpos_+x*x+y*y));
  double leta = -log(tan(theta/2.));
  if (zpos_>0) return leta;
  else return -leta;
  }*/



// std::pair<double,double> HGCSSSimHit::get_xy(const bool isScintillator,
// 					     const HGCSSGeometryConversion & aGeom) const {
std::pair<double,double> HGCSSSimHit::get_xy(const bool radialMapToBeUsed, const HGCSSGeometryConversion & aGeom, unsigned layerCounter) const {
  // if (isScintillator) return aGeom.squareGeom.find(cellid_)->second;
  if (radialMapToBeUsed) {
    // const std::map<int,std::map<int,std::pair<double,double> > > &fhbhGeomsTemp = aGeom.fhbhGeoms;
    // const std::map<int,std::pair<double,double> > &fhbhGeom = fhbhGeomsTemp[layerCounter];
    // const std::map<int,std::pair<double,double> > &fhbhGeom = aGeom.fhbhGeom(layerCounter);
    const std::pair<double,double> phir = aGeom.fhbhGeoms[layerCounter-ANNULARGEOMETRYFIRSTLAYER].find(cellid_)->second;
    double phi = phir.first;
    double r = phir.second;
    std::pair<double,double> xypair(r*cos(phi),r*sin(phi));
    return xypair;
  }
  else return aGeom.hexaGeom.find(cellid_)->second;
}

// ROOT::Math::XYZPoint HGCSSSimHit::position(const bool isScintillator,
// 					   const HGCSSGeometryConversion & aGeom) const{
ROOT::Math::XYZPoint HGCSSSimHit::position(const bool radialMapToBeUsed, const HGCSSGeometryConversion & aGeom, unsigned layerCounter) const {
  std::pair<double,double> xy = get_xy(radialMapToBeUsed, aGeom, layerCounter);
  return ROOT::Math::XYZPoint(xy.first/10.,xy.second/10.,zpos_/10.);
}

double HGCSSSimHit::theta(const bool isScintillator,
			  const HGCSSGeometryConversion & aGeom) const {
  return 2*atan(exp(-1.*eta(isScintillator,aGeom)));
}

double HGCSSSimHit::eta(const bool isScintillator,
			const HGCSSGeometryConversion & aGeom) const {
  bool radialMapToBeUsed = (layer_ >= 3*ANNULARGEOMETRYFIRSTLAYER);
  return position(radialMapToBeUsed,aGeom,layer_/3).eta();
}

double HGCSSSimHit::phi(const bool isScintillator,
			const HGCSSGeometryConversion & aGeom) const {
  bool radialMapToBeUsed = (layer_ >= 3*ANNULARGEOMETRYFIRSTLAYER);
  return position(radialMapToBeUsed,aGeom,layer_/3).phi();
}

void HGCSSSimHit::Print(std::ostream & aOs) const{
  aOs << "====================================" << std::endl
      << " = Layer " << layer() << " siLayer " << silayer() << " cellid " << cellid_ << std::endl
      << " = Energy " << energy_ << " time " << time_ << std::endl
      << " = g " << nGammas_ 
      << " e " << nElectrons_ 
      << " mu " << nMuons_ 
      << " neutron " << nNeutrons_ 
      << " proton " << nProtons_ 
      << " had " << nHadrons_ 
      << std::endl
      << " = main parent: trackID " << trackIDMainParent_ << " efrac " << mainParentEfrac()
      << std::endl
      << "====================================" << std::endl;

}

