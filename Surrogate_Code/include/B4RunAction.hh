// ==================== RunAction.hh Header File =============================== //
// Developed by Bryan V. Egner, Darren E. Holland, and Julie V. Logan
// Modified by Darren Holland 2020-11-02
// ============================================================================= //
// This file initiates the run class
// ============================================================================= //
#ifndef B4RunAction_h         // Required
#define B4RunAction_h 1       // Required
#include "G4UserRunAction.hh" // Required
#include "globals.hh"         // Required
// ============================================================================= //
class G4Run;
class B4RunAction : public G4UserRunAction
{
  public:
    B4RunAction();
    virtual ~B4RunAction();

    virtual void BeginOfRunAction(const G4Run* run);
    virtual void   EndOfRunAction(const G4Run* run);
};
#endif
