// ==================== EventAction.hh Head Filer ============================== //
// Developed by Bryan V. Egner, Darren E. Holland, and Julie V. Logan
// Modified by Darren Holland 2020-11-02
// ============================================================================= //
// This file initiates tracking the energy deposition per thread
// ============================================================================= //
#ifndef B4aEventAction_h
#define B4aEventAction_h 1
#include "G4THitsMap.hh"
#include "G4UserEventAction.hh"
#include "globals.hh"
// ==================== Event Action Class ===================================== //
class B4RunAction;
class B4aEventAction : public G4UserEventAction
{
  public:
    B4aEventAction(B4RunAction* runAction);
    virtual ~B4aEventAction();

    virtual void  BeginOfEventAction(const G4Event* event);
    virtual void    EndOfEventAction(const G4Event* event);
    
    G4double     eventTotEdepDetector;
  private:
    // data members
    G4int  fAbsoEdepHCID;
};
#endif
