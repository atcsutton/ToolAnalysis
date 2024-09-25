#ifndef CutsOptimizer_H
#define CutsOptimizer_H

#include <string>
#include <iostream>

#include "Tool.h"
#include "Hit.h"
#include "Particle.h"
#include "Position.h"

#include "TTree.h"
#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TROOT.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TSystem.h"
#include "TMath.h"

/**
 * \class CutsOptimizer
 *
 * This is a blank template for a Tool used by the script to generate a new custom tool. Please fill out the description and author information.
*
* $Author: D.Ajana $
* $Date: 2024/08/29 10:44:00 $
* Contact: dja23@fsu.edu
*/
class CutsOptimizer: public Tool {


 public:

  CutsOptimizer(); ///< Simple constructor
  bool Initialise(std::string configfile,DataModel &data); ///< Initialise Function for setting up Tool resources. @param configfile The path and name of the dynamic configuration file to read in. @param data A reference to the transient data class used to pass information between Tools.
  bool Execute(); ///< Execute function used to perform Tool purpose.
  bool Finalise(); ///< Finalise function used to clean up resources.

  bool LoadFromStores();
  void SetupTTree();
  void SetupHist();
  void InitHist(double max);
  void WriteHist();
  bool LoadTankClusterClassifiers(double cluster_time);
 private:

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

  TFile *fOutFile;
  TTree *fOutTree, *fOutMCParticles;
  double fTrueVtxX, fTrueVtxY, fTrueVtxZ;
  int fBestPDG, fMoreNeutronQ;
  double fEff, fPur, fTotalQ, fNeutronQ;

  int verbosity;
  int v_error=0;
  int v_warning=1;
  int v_message=2;
  int v_debug=3;
  std::string logmessage;
  Geometry *fGeo = nullptr;
  bool IsInTankMC;
  bool IsInTank;
  int nTotalTrueNeutronsPrompt = 0;
  int nTotalTrueNeutronsDelayed = 0;
  int nTotalNonTrueNeutronsPrompt = 0;
  int nTotalNonTrueNeutronsDelayed = 0;
  int nSelectedTrueNeutronsPrompt = 0;
  int nAllSelectedClusterPrompt = 0;
  int nSelectedTrueNeutronsDelayed = 0;
  int nAllSelectedClusterDelayed = 0;
  bool isPrompt;
  double ChargeTotal = 0;
  double chargeCut = 0;
  double cbCut = 0;

  double fClusterChargeBalance;
  std::map<double, double> cluster_CB;
};


#endif
