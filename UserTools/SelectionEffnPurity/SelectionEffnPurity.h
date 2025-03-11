#ifndef SelectionEffnPurity_H
#define SelectionEffnPurity_H

#include <string>
#include <iostream>
#include <unordered_map>

#include "Tool.h"
#include "Hit.h"
#include "Particle.h"
#include "Position.h"
#include "ADCPulse.h"
#include "CalibratedADCWaveform.h"

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
 * \class SelectionEffnPurity
 *
*
* $Author: D.Ajana $
* $Date: 2024/07/09 10:44:00 $
* Contact: dja23@fsu.edu
*/
class SelectionEffnPurity: public Tool {


 public:

  SelectionEffnPurity(); ///< Simple constructor
  bool Initialise(std::string configfile,DataModel &data); ///< Initialise Function for setting up Tool resources. @param configfile The path and name of the dynamic configuration file to read in. @param data A reference to the transient data class used to pass information between Tools.
  bool Execute(); ///< Execute function used to perform Tool purpose.
  bool Finalise(); ///< Finalise function used to clean up resources.
  bool LoadFromStores();
  void SetupTTree();
  void SetupHist();
  void InitHist(double max);
  void WriteHist();
  //  bool LoadTankClusterClassifiers(double cluster_time);
  
private:
  
    // Configuration variables
  std::string fClusterMapName; // The name of the cluster map in the ANNIEEvent
  std::string fVertexMapName;  // The name of the vertex map in the ANNIEEvent

  // Pointers to load from the ANNIE Event
  std::map<double, std::vector<Hit>> *fClusterMap         = nullptr; // The clusters
  std::map<double, Position>           *fVertexMap          = nullptr; // The vertices
  std::vector<MCParticle>              *fMCParticles        = nullptr; // The true particles from the event
  std::map<int, int>                   *fMCParticleIndexMap = nullptr; // Map between the particle Id and it's position in MCParticles vector
  std::map<std::string,std::vector<double>> fMCNeutCap;

  // Backtracker results
  std::map<double, int>    *fClusterToBestParticlePDG = nullptr; 
  std::map<double, double> *fClusterEfficiency        = nullptr;
  std::map<double, double> *fClusterPurity            = nullptr;
  std::map<double, double> *fClusterTotalCharge       = nullptr;
  std::map<double, double> *fClusterNeutronCharge     = nullptr;
  std::map<double, int> *fClusterToBestParticleIdx = nullptr;
  std::map<double, double> *fClusterEarliestMCTime    = nullptr;
  std::map<double, double> *fClusterMeanMCTime        = nullptr;
  std::map<double, double> *fClusterMedianMCTime      = nullptr;
  std::map<unsigned long, std::vector<MCHit>> *fMCHitsMap = nullptr;
  
  //Experiment
  std::map<int, int> ParticleCountsDelayed;
  std::map<int, int> ParticleCountsPrompt;
  std::map<int, double> map_chankey2spe;
  std::map<int, double> parentIdxToEarliestTime;
  std::map<int, double> parentIdxToWeightedTimeNumerator;
  std::map<int, double> parentIdxToTotalCharge;
  
  // Output ROOT file things
  TFile *fOutFile;
  TTree *fOutTree;
  TTree *fOutTreeContDelayed;
  TTree *fOutTreeContPrompt;
  double fTrueVtxX, fTrueVtxY, fTrueVtxZ;
  double fMCX, fMCY, fMCZ;
  double fRecoVtxX, fRecoVtxY, fRecoVtxZ;
  double fDistX, fDistY, fDistZ, fDist; 
  int fBestPDG, fMoreNeutronQ;
  double fEff, fPur, fTotalQ, fNeutronQ;

  //combination delayed and prompt 
  TH2F *h_nSelectedTrueNeutronsTVtxXZ = nullptr;
  TH2F *h_nTotalTrueNeutronsTVtxXZ = nullptr;
  //Histograms
  //Neutrino Energy
  TH1F *h_nTotalTrueNeutronsPromptNE       = nullptr;
  TH1F *h_nTotalTrueNeutronsDelayedNE     = nullptr;
  TH1F *h_nAllSelectedClustersPromptNE = nullptr;
  TH1F *h_nSelectedTrueNeutronsPromptNE = nullptr;
  TH1F *h_nAllSelectedClustersDelayedNE = nullptr;
  TH1F *h_nSelectedTrueNeutronsDelayedNE = nullptr;

  //True vertex X-Z
  TH2F *h_nTotalTrueNeutronsPromptTVtxXZ = nullptr;
  TH2F *h_nTotalTrueNeutronsDelayedTVtxXZ = nullptr;
  TH2F *h_nAllSelectedClustersPromptTVtxXZ = nullptr;
  TH2F *h_nSelectedTrueNeutronsPromptTVtxXZ = nullptr;
  TH2F *h_nAllSelectedClustersDelayedTVtxXZ = nullptr;
  TH2F *h_nSelectedTrueNeutronsDelayedTVtxXZ = nullptr;

  //True vertex X-Y
  TH2F *h_nTotalTrueNeutronsPromptTVtxXY = nullptr;
  TH2F *h_nTotalTrueNeutronsDelayedTVtxXY = nullptr;
  TH2F *h_nAllSelectedClustersPromptTVtxXY = nullptr;
  TH2F *h_nSelectedTrueNeutronsPromptTVtxXY = nullptr;
  TH2F *h_nAllSelectedClustersDelayedTVtxXY = nullptr;
  TH2F *h_nSelectedTrueNeutronsDelayedTVtxXY = nullptr;

  //True vertex Y-Z
  TH2F *h_nTotalTrueNeutronsPromptTVtxYZ = nullptr;
  TH2F *h_nTotalTrueNeutronsDelayedTVtxYZ = nullptr;
  TH2F *h_nAllSelectedClustersPromptTVtxYZ = nullptr;
  TH2F *h_nSelectedTrueNeutronsPromptTVtxYZ = nullptr;
  TH2F *h_nAllSelectedClustersDelayedTVtxYZ = nullptr;
  TH2F *h_nSelectedTrueNeutronsDelayedTVtxYZ = nullptr;
  //Cluster Time
  TH1F *h_nTotalTrueNeutronsPromptCT = nullptr;
  TH1F *h_nTotalTrueNeutronsDelayedCT = nullptr;
  TH1F *h_nAllSelectedClustersPromptCT = nullptr;
  TH1F *h_nSelectedTrueNeutronsPromptCT = nullptr;
  TH1F *h_nAllSelectedClustersDelayedCT = nullptr;
  TH1F *h_nSelectedTrueNeutronsDelayedCT = nullptr;

  //Particle PDGs
  TH1F *h_nTotalTrueNeutronsPromptPDG = nullptr;
  TH1F *h_nTotalTrueNeutronsDelayedPDG = nullptr;
  TH1F *h_nAllSelectedClustersPromptPDG = nullptr;
  TH1F *h_nSelectedTrueNeutronsPromptPDG = nullptr;
  TH1F *h_nAllSelectedClustersDelayedPDG = nullptr;
  TH1F *h_nSelectedTrueNeutronsDelayedPDG = nullptr;

  //N hits
  TH1F *h_nTotalTrueNeutronsPromptNhits = nullptr;
  TH1F *h_nTotalTrueNeutronsDelayedNhits = nullptr;
  TH1F *h_nAllSelectedClustersPromptNhits = nullptr;
  TH1F *h_nSelectedTrueNeutronsPromptNhits = nullptr;
  TH1F *h_nAllSelectedClustersDelayedNhits = nullptr;
  TH1F *h_nSelectedTrueNeutronsDelayedNhits = nullptr;

  TH1F *h_allContaminationPrompt = nullptr;
  TH1F *h_allContaminationDelayed = nullptr;

  TH1F *h_nTotalTrueNeutronsPromptCC = nullptr;

  //ClusterCharge vs ChargeBalance
  TH1F *h_nAllSelectedClusterPromptClusterCharge = nullptr;
  TH1F *h_nAllSelectedClusterPromptChargeBalance = nullptr;
  TH1F *h_nSelectedTrueNeutronsPromptClusterCharge = nullptr;
  TH1F *h_nSelectedTrueNeutronsPromptBalance = nullptr;

  TH1F *h_nAllSelectedClusterDelayedClusterCharge = nullptr;
  TH1F *h_nAllSelectedClusterDelayedChargeBalance = nullptr;
  TH1F *h_nSelectedTrueNeutronsDelayedClusterCharge = nullptr;
  TH1F *h_nSelectedTrueNeutronsDelayedBalance = nullptr;

  TH2F *h_nAllSelectedClusterPromptCBCC = nullptr;
  TH2F *h_nSelectedTrueNeutronsPromptCBCC = nullptr;

  TH2F *h_nAllSelectedClusterDelayedCBCC = nullptr;
  TH2F *h_nSelectedTrueNeutronsDelayedCBCC = nullptr;
  /// \brief verbosity levels: if 'verbosity' < this level, the message type will be logged.
  int verbosity;
  int v_error=0;
  int v_warning=1;
  int v_message=2;
  int v_debug=3;
  std::string logmessage;
  Geometry *fGeo = nullptr;
  
  int nSelectedTrueNeutronsPrompt  = 0;                                      
  int nTotalTrueNeutrons           = 0;                        
  int nTotalTrueNeutronsPrompt     = 0;
  int nTotalTrueNeutronsDelayed    = 0;
  int nAllSelectedClustersPrompt   = 0;                                   
  int nTotalTrueNeutronsWorld      = 0; //Total True Neutron in MCParticle before geometry cut.
  int nAllSelectedClustersWorld    = 0; //All selected cluster in the backtracker before charge balance etc cut
  int nAllSelectedClustersDelayed  = 0;
  int nSelectedTrueNeutronsDelayed = 0;
  int nTotalTrueNeutronsDelayedReq = 0;
  int nTotalTrueNeutronsDelayedMichel = 0;
  int nSelectedTrueNeutronsDelayedMichel = 0;
  int nSelectedTrueNeutronsDelayedReq = 0;
  int promptTrueVisibleNeutrons = 0;
  int delayedTrueVisibleNeutrons = 0;
  int totalTrueVisibleNeutrons = 0;
  //experimenting
  int nSumingAllTrueNeutron = 0;
  int nSumingTotalTrueNeutron = 0;
  
  
  double fClusterChargeBalance;
  std::map<double, double> ClusterChargeBalances;
  //  std::map<double, double> cluster_CB;


};


#endif
