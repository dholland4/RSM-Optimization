// ==================== SteppingAction.cc Class ================================= //
// Developed by Bryan V. Egner, Darren E. Holland, and Julie V. Logan
// Modified by Darren Holland 2020-11-02
// ============================================================================= //
// Doesn't need to do anything for the surrogate model
// ============================================================================= //
#include "B4aSteppingAction.hh"
#include "B4aEventAction.hh"
#include "B4DetectorConstruction.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"
//
#include "G4Step.hh"
#include "G4RunManager.hh"
//
#include "Randomize.hh"
#include <iomanip>
#include <iostream>
#include <fstream>
#include <stdio.h> 
#include <math.h>
#include <Settings.hh>
//
#include "G4UImanager.hh"
#include "G4UIcommand.hh"

using namespace std;

// Start step instance
B4aSteppingAction::B4aSteppingAction(B4aEventAction* eventAction) : G4UserSteppingAction(),fEventAction(eventAction)
{}

// Destroy step instance when completed
B4aSteppingAction::~B4aSteppingAction()
{}

// Options to occur during event
void B4aSteppingAction::UserSteppingAction(const G4Step* theStep)
{}