// ==================== DetectorConstruction.hh Head File ====================== //
// Developed by Bryan V. Egner, Darren E. Holland, and Julie V. Logan
// Modified by Darren Holland 2020-11-02
// ============================================================================= //
// This file is the mandatory initialization class for detector setup.
// ============================================================================= //
#ifndef B4DetectorConstruction_h  // Required
#define B4DetectorConstruction_h 1  // Required
#include "G4VUserDetectorConstruction.hh"
#include "globals.hh"
// ============================================================================= //
class G4VPhysicalVolume;
class B4DetectorConstruction : public G4VUserDetectorConstruction
{
  public:
    B4DetectorConstruction();
    virtual ~B4DetectorConstruction();
  public:  
    virtual G4VPhysicalVolume* Construct();
    // Detector class (NOTE won't work properly if function name is changed)
    virtual void ConstructSDandField();
  private:
    G4VPhysicalVolume* DefineVolumes();
    // Option to activate checking of volumes overlaps
    G4bool  fCheckOverlaps;
};
#endif
