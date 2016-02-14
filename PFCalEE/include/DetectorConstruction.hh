#ifndef DetectorConstruction_h
#define DetectorConstruction_h 1

#include "SamplingSection.hh"

#include "G4VUserDetectorConstruction.hh"
#include "globals.hh"

#include <map>
#include <string>

class G4CSGSolid;
class G4Box;
class G4LogicalVolume;
class G4VPhysicalVolume;
class G4Material;
class G4UniformMagField;
class DetectorMessenger;
class G4Colour;

/**
   @class DetectorConstruction
   @short builds a simple detector
 */
class DetectorConstruction : public G4VUserDetectorConstruction
{
public:

  enum DetectorVersion { 
    v_CALICE=0,
    v_HGCALEE_Si80=1,
    v_HGCALEE_Si120=2,
    v_HGCALEE_Si200=3,
    v_HGCALEE_Si500=4,
    v_HGCALEE_gap1=5,
    v_HGCALEE_CALICE=6,
    v_HGCALEE_inverted=7,
    v_HGCALEE_concept=8,
    v_HGCALEE_W=9,
    v_HGCALEE_gap4=10,
    v_HGCALEE_prePCB=11,
    v_HGCALEE_v5=12,
    v_HGCALEE_v5_gap4=13,
    v_HGCAL=20,
    v_HGCALHE=21,
    v_HGCALHEScint=22,
    v_HGCALHE_CALICE=23,
    v_HGCALHE_CMSSWv4=24,
    v_HGCAL_v5=25,
    v_HGCAL_v5_gap4=26,
    v_HGCALHE_v5=27,
    v_HGCALBE_v5=28,
    v_HGCALEE_v6=30,
    v_HGCALHE_v6=31,
    v_HGCALBE_v6=32,
    v_HGCAL_v6=33,
    v_HGCALEE_v624_cnncal=34,
    v_HGCALEE_v618=35,
    v_HGCAL_v624=36,
    v_HGCAL_v618=37,
    v_HGCALHE_v624=38,
    v_HGCALHE_v618=39,
    v_HGCALEE_v624_frac1=40,
    v_HGCALEE_v624_frac75=41,
    v_HGCALEE_v624_frac50=42,
    v_HGCALEE_v624_frac25=43,
    v_HGCALEE_v624_flat=44,
    v_HGCALEE_v624_syst_opt=45
  };

  enum DetectorModel {
    m_SIMPLE_20=0,
    m_SIMPLE_50=1,
    m_FULLSECTION=2,
    m_SIMPLE_100=3
  };

  /**
     @short CTOR
   */
  DetectorConstruction(G4int ver=DetectorConstruction::v_CALICE, 
		       G4int mod=DetectorConstruction::m_SIMPLE_20,
		       std::string absThickW="2.6,2.6,2.6,2.6,3.6,3.6,3.6,3.6,4.2,4.2,4.2,4.2",
		       std::string absThickWCu="1.0,1.0,1.0,1.0,1.75,1.75,1.75,1.75,2.2,2.2,2.2,2.2",
		       std::string absThickPb="1,1,1,1,1,2.1,2.1,2.1,2.1,2.1,4.4,4.4,4.4,4.4",
		       std::string dropLayer="");

  void buildHGCALFHE(const unsigned aVersion);
  void buildHGCALBHE(const unsigned aVersion);
  /**
     @short calorimeter structure (sampling sections)
   */
  std::vector<SamplingSection> m_caloStruct;
  std::vector<SamplingSection> *getStructure() { return &m_caloStruct; }

  int getModel() const { return model_; }
  int getVersion() const { return version_; }

  const std::vector<G4LogicalVolume*>  & getSiLogVol() {return m_logicSi; }
  const std::vector<G4LogicalVolume*>  & getAbsLogVol() {return m_logicAbs; }


  /**
     @short define the calorimeter materials
   */
  void DefineMaterials(); 
  std::map<std::string, G4Material *> m_materials;
  std::map<std::string, G4double > m_dEdx;
  std::map<std::string, G4Colour > m_colours;
  
  /**
     @short set magnetic field
   */
  void SetMagField(G4double fieldValue);
  G4UniformMagField* m_magField;      //pointer to the magnetic field

  /**
     @short set detector model
   */

  void SetDetModel(G4int model);

  void SetWThick(std::string thick);
  void SetWCuThick(std::string thick);
  void SetPbThick(std::string thick);
  void SetDropLayers(std::string layers);

  /**
     @short DTOR
   */
  ~DetectorConstruction();
  

  /**
     @short getters
   */
  G4double GetCalorSizeXY() { return m_CalorSizeXY; }
  G4double GetCalorSizeZ()  { return m_CalorSizeZ; }
  G4double GetWorldSizeXY() { return m_WorldSizeXY; }
  G4double GetWorldSizeZ()  { return m_WorldSizeZ; }

  /**
     @short build the detector
   */

  G4VPhysicalVolume* Construct();

private:

  //detector version
  int version_;
  //integer to define detector model
  int model_;

  //add a pre PCB plate
  bool addPrePCB_;

  std::vector<G4double> absThickW_;
  std::vector<G4double> absThickWCu_;
  std::vector<G4double> absThickPb_;
  std::vector<G4bool> dropLayer_;

  /**
     @short compute the calor dimensions
   */
  void UpdateCalorSize(); 

  /**
     @short build the calorimeter
   */
  G4VPhysicalVolume* ConstructCalorimeter();     

  G4CSGSolid *constructSolid (std::string baseName, G4double thick, G4double zpos);
 
  std::vector<G4Material* > m_SensitiveMaterial;
  
  G4double           m_CalorSizeXY, m_CalorSizeZ;
  G4double           m_minRadius,m_maxRadius;
  G4double           m_minEta,m_maxEta;
  G4double           m_z0pos;
  G4double           m_WorldSizeXY, m_WorldSizeZ;

            
  G4CSGSolid*        m_solidWorld;    //pointer to the solid World 
  G4LogicalVolume*   m_logicWorld;    //pointer to the logical World
  G4VPhysicalVolume* m_physWorld;     //pointer to the physical World  
  
  std::vector<G4LogicalVolume*>   m_logicSi;    //pointer to the logical Si volumes
  std::vector<G4LogicalVolume*>   m_logicAbs;    //pointer to the logical absorber volumes situated just before the si

  DetectorMessenger* m_detectorMessenger;  //pointer to the Messenger
};


#endif

