// ==================== ActionInitialization.hh Head File ====================== //
// Developed by Bryan V. Egner, Darren E. Holland, and Julie V. Logan
// Modified by Darren Holland 2020-11-02
// ============================================================================= //
// This file calls/initializes the following classes:
//      - B4PrimaryGeneratorAction
//      - B4RunAction
//      - B4aEventAction
//      - B4aSteppingAction
// ============================================================================= //
#ifndef B4aActionInitialization_h
#define B4aActionInitialization_h 1
#include "G4VUserActionInitialization.hh"
// ========== Detector Construction Class:
class B4DetectorConstruction;
class PrimaryGeneratorAction;
// ========== Action Initialization Class:
class B4aActionInitialization : public G4VUserActionInitialization
{
  public:
    B4aActionInitialization();
    virtual ~B4aActionInitialization();
    virtual void Build() const;
};
#endif
