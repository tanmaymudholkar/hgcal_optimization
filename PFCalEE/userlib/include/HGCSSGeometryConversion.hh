#ifndef HGCSSGeometryConversion_h
#define HGCSSGeometryConversion_h


#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <cmath>
#include <cstdlib>
#include "TH2D.h"
#include "TH2Poly.h"
#include "TMath.h"
#include "HGCSSDetector.hh"
#include "HGCSSHardcodedConstants.hh"



struct MergeCells {
  double energy;
  double time;
  //double z;
};

class HGCSSGeometryConversion{
  
public:
  HGCSSGeometryConversion(){};
  HGCSSGeometryConversion(const unsigned model, const double cellsize, const bool bypassR=false, const unsigned nSiLayers=3);

  ~HGCSSGeometryConversion();

  TH2Poly *hexagonMap(){
    static TH2Poly hc;
    return &hc;
  };

  // TH2Poly *squareMap(){
  //   static TH2Poly hsq;
  //   return &hsq;
  // };

  // TH2Poly *fhbhMaps(unsigned layerCounter) {
  //   static std::map<int, TH2Poly> hfhbh;
  //   return &(hfhbh[layerCounter]);
  // }

  TH2Poly *fhbhMaps(unsigned fhbhCounter) {
    static TH2Poly hfhbh[NFHBHLAYERS];
    return &(hfhbh[fhbhCounter]);
  }

  std::map<int,std::pair<double,double> > hexaGeom;
  // std::map<int,std::pair<double,double> > squareGeom;
  // std::map<int,std::map<int,std::pair<double,double> > > fhbhGeoms; // first index for layer counter
  std::map<int,std::pair<double,double> > fhbhGeoms[NFHBHLAYERS]; // first index for layer counter
  
  // void initialiseSquareMap(const double xymin, const double side);

  // void initialiseSquareMap(TH2Poly *map, const double xymin, const double side, bool print);

  std::vector<Double_t> readInGeometryFromFile(std::string inputFileName);

  void initialiseFHBHMaps(std::string geometryInputFolder);

  void initialiseFHBHMap(TH2Poly *map, std::vector<Double_t> geometryVector, unsigned nPhiDivisions);

  void initialiseHoneyComb(const double xymin, const double side);

  void initialiseHoneyComb(TH2Poly *map, const double xymin, const double side, bool print);

  void fillXY(TH2Poly* hist, std::map<int,std::pair<double,double> > & geom);

  void setGranularity(const std::vector<unsigned> & granul);

  unsigned getGranularity(const unsigned aLayer, const HGCSSSubDetector & adet);

  // const std::map<int,std::pair<double,double> > & fhbhGeom(unsigned layerCounter) {
  //   return fhbhGeoms[layerCounter];
  // }

  inline double getXYwidth() const {
    return width_;
  };
  
  inline void setXYwidth(double width) {
    width_ = width;
  };
  
  inline double cellSize() const{
    return cellSize_;
  };

  //hardcode fine granularity at high eta ?
  //  inline double cellSize(const unsigned aLayer, const double aEta) const{
  //  if (fabs(aEta)<10) 
  //    return cellSize_*granularity_[aLayer];
  //  return cellSize_*3;
  //};
  //  inline double cellSizeInCm(const unsigned aLayer, const double aEta) const{
  // return cellSize(aLayer, aEta)/10.;
  //};
  double cellSize(const unsigned aLayer, const double aR) const;

  double cellSizeInCm(const unsigned aLayer, const double aR) const;

  //unsigned getNumberOfSiLayers(const DetectorEnum type,
  //const double & eta=0) const;
  unsigned getNumberOfSiLayers(const DetectorEnum type,
			       const double & radius=10000) const;

  void initialiseHistos(const bool recreate=false,
			std::string uniqStr="",
			const bool print=true);

  //void fill(const DetectorEnum type,
  //const unsigned newlayer,
  //const double & weightedE,
  //const double & aTime,
  //const double & posx,
  //const double & posy,
  //const double & posz);

  void fill(const unsigned layer,
	    const double & weightedE,
	    const double & aTime,
	    const unsigned & cellid,
	    const double & posz);


  double getAverageZ(const unsigned layer);

  //with TH2Poly
  /*double sumBins(const std::vector<TH2Poly *> & aHistVec,
		 const double & aMipThresh=0.);

  void resetVector(std::vector<TH2Poly *> & aVec,
		   std::string aVar,
		   std::string aString,
		   const HGCSSSubDetector & aDet,
		   const unsigned nLayers,
		   bool recreate=false,
		   bool print=true);


  void deleteHistos(std::vector<TH2Poly *> & aVec);

  TH2Poly * get2DHist(const unsigned layer,std::string name);

  inline std::vector<TH2Poly *> & get2DEnergyVec(const DetectorEnum aDet){
    return HistMapE_[aDet];
  };

  inline std::vector<TH2Poly *> & get2DTimeVec(const DetectorEnum aDet){
    return HistMapTime_[aDet];
  };

  inline std::vector<TH2Poly *> & get2DZposVec(const DetectorEnum aDet){
    return HistMapZ_[aDet];
  };
  */

  //with maps of MergeCells struct
  double sumBins(const std::map<unsigned,MergeCells> & aHistVec,
		 const double & aMipThresh=0.);

  void resetVector(std::map<unsigned,MergeCells> & aVec);

  void deleteHistos(std::map<unsigned,MergeCells> & aVec);

  std::map<unsigned,MergeCells> & get2DHist(const unsigned layer);

private:

void myHoneycomb(TH2Poly* map,
		 Double_t xstart,
		 Double_t ystart,
		 Double_t a,  // side length
		 Int_t k,     // # hexagons in a column
		 Int_t s);    // # columns

  bool dopatch_;
  double width_;
  double cellSize_;
  std::vector<unsigned> granularity_;
  unsigned model_;
  bool bypassRadius_;
  unsigned nSiLayers_;
  //std::map<DetectorEnum,std::vector<TH2Poly *> > HistMapE_;
  //std::map<DetectorEnum,std::vector<TH2Poly *> > HistMapTime_;
  //std::map<DetectorEnum,std::vector<TH2Poly *> > HistMapZ_;
  //std::map<DetectorEnum,std::vector<double> > avgMapZ_;
  //std::map<DetectorEnum,std::vector<double> > avgMapE_;
  std::map<unsigned,std::map<unsigned,MergeCells> > HistMap_;
  std::map<unsigned,double> avgMapZ_;
  std::map<unsigned,double> avgMapE_;



};



#endif
