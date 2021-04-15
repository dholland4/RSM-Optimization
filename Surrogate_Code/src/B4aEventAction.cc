// ==================== EventAction.cc Class =================================== //
// Developed by Bryan V. Egner, Darren E. Holland, and Julie V. Logan
// Modified by Darren Holland 2020-11-02
// ============================================================================= //
// This file doesn't need to track anything for the surrogate model
// ============================================================================= //
#include "B4aEventAction.hh"
#include "B4RunAction.hh"
#include "G4Event.hh"
//
#include "G4RunManager.hh"
#include "G4Event.hh"
#include "G4Run.hh"
#include "G4SDManager.hh"
#include "G4HCofThisEvent.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"
//
#include "Randomize.hh"
#include <iomanip>
#include <iostream>
#include <fstream>
#include <stdio.h> 
#include <math.h>
//
#include <Settings.hh>

using namespace std;

// Create instance
B4aEventAction::B4aEventAction(B4RunAction* runAction)
 : G4UserEventAction(),
   fAbsoEdepHCID(-1),  // Absolute Energy Deposition
 eventTotEdepDetector(0.)
{}

// When finished, destroy the instance
B4aEventAction::~B4aEventAction()
{}

// Options for beginning of event
void B4aEventAction::BeginOfEventAction(const G4Event* event)
{}

//Options for end of event
void B4aEventAction::EndOfEventAction(const G4Event* event)
{}