// ==================== RunAction.cc Class ===================================== //
// Developed by Bryan V. Egner, Darren E. Holland, and Julie V. Logan
// Modified by Darren Holland 2020-11-02
// ============================================================================= //
// This file starts a global run.
// ============================================================================= //
#include "B4RunAction.hh"
#include "G4Run.hh"
#include "G4RunManager.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"
#include "B4PrimaryGeneratorAction.hh"
#include "G4SDManager.hh"
#include <iomanip>
#include <iostream>
#include <fstream>
#include <stdio.h>
#include "Settings.hh"

using namespace std;

// Create instance
B4RunAction::B4RunAction()
 : G4UserRunAction() 
{ 
  G4RunManager::GetRunManager()->SetPrintProgress(0);
}

// When completed, destroy instance
B4RunAction::~B4RunAction()
{}

// Code for beginning of run
void B4RunAction::BeginOfRunAction(const G4Run* /*run*/)
{}

// Code for end of run
void B4RunAction::EndOfRunAction(const G4Run* run)
{}