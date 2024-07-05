#ifndef EvaluateVertex_H
#define EvaluateVertex_H

#include <string>
#include <iostream>

#include "Tool.h"
#include "Hit.h"
#include "Particle.h"
#include "Position.h"

#include "TTree.h"
#include "TFile.h"


/**
 * \class EvaluateVertex
 *
 * This is a blank template for a Tool used by the script to generate a new custom tool. Please fill out the description and author information.
*
* $Author: B.Richards $
* $Date: 2019/05/28 10:44:00 $
* Contact: b.richards@qmul.ac.uk
*/
class EvaluateVertex: public Tool {


 public:

  EvaluateVertex(); ///< Simple constructor
  bool Initialise(std::string configfile,DataModel &data); ///< Initialise Function for setting up Tool resources. @param configfile The path and name of the dynamic configuration file to read in. @param data A reference to the transient data class used to pass information between Tools.
  bool Execute(); ///< Execute function used to perform Tool purpose.
  bool Finalise(); ///< Finalise function used to clean up resources.

  bool LoadFromStores(); ///< Does all the loading so I can move it away from the Execute function
  void SetupTTree();

 private:

  // Configuration variables
  std::string fClusterMapName; // The name of the cluster map in the ANNIEEvent
  std::string fVertexMapName;  // The name of the vertex map in the ANNIEEvent

  // Pointers to load from the ANNIE Event
  std::map<double, std::vector<MCHit>> *fClusterMap         = nullptr; // The clusters
  std::map<double, Position>           *fVertexMap          = nullptr; // The vertices
  std::vector<MCParticle>              *fMCParticles        = nullptr; // The true particles from the event
  std::map<int, int>                   *fMCParticleIndexMap = nullptr; // Map between the particle Id and it's position in MCParticles vector

  // Backtracker results
  std::map<double, int>    *fClusterToBestParticleID  = nullptr;
  std::map<double, int>    *fClusterToBestParticlePDG = nullptr; 
  std::map<double, double> *fClusterEfficiency        = nullptr;
  std::map<double, double> *fClusterPurity            = nullptr;
  std::map<double, double> *fClusterTotalCharge       = nullptr;
  std::map<double, double> *fClusterNeutronCharge     = nullptr;

  // Output ROOT file things
  TFile *fOutFile;
  TTree *fOutTree;
  double fTrueVtxX, fTrueVtxY, fTrueVtxZ;
  double fRecoVtxX, fRecoVtxY, fRecoVtxZ;
  double fDistX, fDistY, fDistZ, fDist; 
  int fBestPDG, fMoreNeutronQ;
  double fEff, fPur, fTotalQ, fNeutronQ;


  /// \brief verbosity levels: if 'verbosity' < this level, the message type will be logged.
  int verbosity;
  int v_error=0;
  int v_warning=1;
  int v_message=2;
  int v_debug=3;
  std::string logmessage;

};


#endif
