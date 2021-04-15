// ==================== ActionInitialization.cc Class ========================== //
// Developed by Bryan V. Egner, Darren E. Holland, and Julie V. Logan
// Modified by Darren Holland 2020-11-02
// ============================================================================= //
// This file calls/initializes the following classes:
//	- B4PrimaryGeneratorAction
//	- B4RunAction
//	- B4aEventAction
//	- B4aSteppingAction
// ============================================================================= //
#include "B4aActionInitialization.hh"
#include "B4PrimaryGeneratorAction.hh"
#include "B4RunAction.hh"
#include "B4aEventAction.hh"
#include "B4aSteppingAction.hh"

// Function to create class
B4aActionInitialization::B4aActionInitialization()
 : G4VUserActionInitialization()
{}

// Function to destroy class when finished
B4aActionInitialization::~B4aActionInitialization()
{}

void B4aActionInitialization::BuildForMaster() const
{
  SetUserAction(new B4RunAction);
}

void B4aActionInitialization::Build() const
{
  // Initialize particle generator
  SetUserAction(new B4PrimaryGeneratorAction);
  // Start run (each "global" run is a given source position)
  B4RunAction* runAction = new B4RunAction();
  SetUserAction(runAction);
  // Find events for a local run (single particle) and record total energy deposited in detector
  B4aEventAction* eventAction = new B4aEventAction(runAction);
  SetUserAction(eventAction);
  // Step through local run
  SetUserAction(new B4aSteppingAction(eventAction));
  // Nothing occurs after stepaction, and so it loops back to the event action
}  
