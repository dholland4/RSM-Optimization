// ==================== SteppingAction.hh Header File =========================== //
// Developed by Bryan V. Egner, Darren E. Holland, and Julie V. Logan
// Modified by Darren Holland 2020-11-02
// ============================================================================= //
// This file initiates the stepping class
// ============================================================================= //
#ifndef B4aSteppingAction_h
#define B4aSteppingAction_h 1
#include "G4UserSteppingAction.hh"
// ============================================================================= //
class B4DetectorConstruction;
class B4aEventAction;
// ============================================================================= //
class B4aSteppingAction : public G4UserSteppingAction
{
public:
  B4aSteppingAction(B4aEventAction* eventAction);
  virtual ~B4aSteppingAction();
  virtual void UserSteppingAction(const G4Step* step);
private:
  B4aEventAction*  fEventAction;  
};
#endif
// ============================================================================= //
