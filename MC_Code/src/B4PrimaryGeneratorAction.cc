// ==================== PrimaryGeneratorAction.cc Class ======================== //
// Developed by Bryan V. Egner, Darren E. Holland, and Julie V. Logan
// Modified by Darren Holland 2020-11-02
// ============================================================================= //
// This file generates the particle with the desired energy/energies and 
// releases it into a cone toward the detector
// ============================================================================= //
#include "B4PrimaryGeneratorAction.hh"
#include "G4RunManager.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4LogicalVolume.hh"
#include "G4Box.hh"
#include "G4Orb.hh"
#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh"
#include "G4GeneralParticleSource.hh"
#include "Settings.hh"

// Create particle gun
B4PrimaryGeneratorAction::B4PrimaryGeneratorAction()
 : G4VUserPrimaryGeneratorAction(),
   fParticleGun(0)
{
  // Source Characterstics:
  fParticleGun = new G4GeneralParticleSource;
  // Load particle type (gamma, neutron)
  G4ParticleDefinition* particleDefinition = G4ParticleTable::GetParticleTable()->FindParticle(Settings::PartType); // Particle Type
  fParticleGun->SetParticleDefinition(particleDefinition);  // Particle Definition
  // Only one particle at a time
  fParticleGun->SetNumberOfParticles (1);  // Number of Particles
}

// Destroy instance
B4PrimaryGeneratorAction::~B4PrimaryGeneratorAction()
{
  delete fParticleGun;
}

void B4PrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
  // This function is called at the begining of local run (aka begin by shooting particle)
  // Set particle direction
  fParticleGun->GeneratePrimaryVertex(anEvent);
}
