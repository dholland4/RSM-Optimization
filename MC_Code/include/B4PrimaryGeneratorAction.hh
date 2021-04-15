// ==================== PrimaryGeneratorAction.hh Header File ================== //
// Developed by Bryan V. Egner, Darren E. Holland, and Julie V. Logan
// Modified by Darren Holland 2020-11-02
// ============================================================================= //
// This file initiates the particle source generation
// ============================================================================= //
#ifndef B4PrimaryGeneratorAction_h
#define B4PrimaryGeneratorAction_h 1
#include "G4VUserPrimaryGeneratorAction.hh"
#include "G4ParticleGun.hh"
#include "globals.hh"
#include "G4GeneralParticleSource.hh"  // Create General Particle
// ============================================================================= //
class G4ParticleGun;
class G4Event;
class B4PrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
{
public:
  B4PrimaryGeneratorAction();    
  virtual ~B4PrimaryGeneratorAction();
  virtual void GeneratePrimaries(G4Event* event);
  // Set methods
  void SetRandomFlag(G4bool value);
  // Method to access particle gun:
  const G4GeneralParticleSource* GetParticleGun() const { return fParticleGun; }
private:
  G4GeneralParticleSource*  fParticleGun; // G4 particle gun
};
#endif
