// ================= Settings Header File (important variables) ================ //
// Developed by Bryan V. Egner, Darren E. Holland, and Julie V. Logan
// Modified by Darren Holland 2020-11-02
// ============================================================================= //
// This file initiates the settings needed for the Geant geometry creation and
// subseqent analysis
// ============================================================================= //
#ifndef Settings_h
#define Settings_h
#include <vector>
// ============================================================================= //
namespace Settings
{
	// ========== Important Parameters:(Length(cm), Energy(MeV)
	extern const int RSMTet;		    // RSM Tet or Tess
	extern std::string fname_nodes;		// Mask nodes filename
	extern std::string fname_ele;		// Mask element filename
	extern std::string fname_out;		// Output filename
    extern const int nParts2Run;		// Number of Particles
    extern std::string SourceEnergyType;// Source Energy Type
	extern std::string PartType;		// Particle Type
    extern const int SourceDiv;		    // Source position sub-divisions
    extern const double SourceDist;		// Source Distance
    extern const double coneangle;		// Cone Angle
    extern const int numEnergies;		// Number of Energies (remove)
    extern const double energiesMeV;	// Energy in MeV
    extern const double deltatheta;     // Theta increment
    extern const double deltaphi;	    // Phi increment
    extern const double DetRad;       	// Detector radius
    extern const double DetHeight;    	// Detector half height
    extern const double SleeveOuterRad; // Sleeve radius
    extern const double SleeveHeight;   // Top of sleeve (extends past detector)
    extern const double SleeveBottom;   // Total length of sleeve
    extern const double StartPhi;	    // Initial source angle
    extern const double EndPhi;		    // Final source angle    
}
#endif
