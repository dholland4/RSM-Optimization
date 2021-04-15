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
#include <vector>
#include <string>
#include "G4Threading.hh"
#include <boost/filesystem.hpp>
#include <thread>
// ==================== Event Action Class ===================================== //
using namespace std;
class B4RunAction;
class B4aEventAction : public G4UserEventAction
{
  public:
    B4aEventAction(B4RunAction* runAction);
    virtual ~B4aEventAction();

    virtual void  BeginOfEventAction(const G4Event* event);
    virtual void    EndOfEventAction(const G4Event* event);
    
    // Set number of threads
    static const G4int NumThreads = 64;
    // Track run number, total energy deposited, and event number
    G4int runNum[NumThreads];
    G4double eventTotEdepDetector[NumThreads]; 
    G4int Evt[NumThreads];
    // Sum energy deposited
    void AddEdepDetector(G4double edep, G4int Threadindex) {
         eventTotEdepDetector[Threadindex] += edep; };
    // Track thread number to save to unique file (avoids race condition)
    std::vector<boost::filesystem::path> ThreadNum;
  private:
    // Get file names and path
    void GetFilesOfTypeInDirectory(const boost::filesystem::path &directoryPath, const string &fileExtension, std::vector<boost::filesystem::path> &list);
    bool PathSort(const boost::filesystem::path &first, const boost::filesystem::path &second);
    void PrintEventStatistics(G4double absoEdep) const;
    // Write results to file
    virtual void writetofile(G4String fname, const G4double IPE, const G4double Edep, G4int runNum, G4int eid);
};
#endif
