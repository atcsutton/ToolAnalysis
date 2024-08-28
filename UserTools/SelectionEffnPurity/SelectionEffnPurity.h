#ifndef SelectionEffnPurity_H
#define SelectionEffnPurity_H

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
  bool LoadTankClusterClassifiers(double cluster_time);
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
  TTree *fOutTreeContDelayed;
  TTree *fOutTreeContPrompt;
  double fTrueVtxX, fTrueVtxY, fTrueVtxZ;
  double fRecoVtxX, fRecoVtxY, fRecoVtxZ;
  double fDistX, fDistY, fDistZ, fDist; 
  int fBestPDG, fMoreNeutronQ;
  double fEff, fPur, fTotalQ, fNeutronQ;

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

  int nPromptProton = 0;
  int nPromptAntiProton = 0;
  int nPromptElectron = 0;
  int nPromptPositron = 0;
  int nPromptElectronNeutrino = 0;
  int nPromptAntiElectronNeutrino = 0;
  int nPromptGamma = 0;
  int nPromptNeutron = 0;
  int nPromptAntiNeutron = 0;
  int nPromptMuonPlus = 0;
  int nPromptMuonMinus = 0;
  int nPromptKaonlong = 0;
  int nPromptPionPlus = 0;
  int nPromptPionMinus = 0;
  int nPromptKaonPlus = 0;
  int nPromptKaonMinus = 0;
  int nPromptKaonshort = 0;
  int nPromptPion0 = 0;
  int nPromptKaon0 = 0;
  int nPromptMuonNeutrino = 0;
  int nPromptAntiMuonNeutrino = 0;
  int nPromptTauPlus = 0;
  int nPromptTauMinus = 0;


  int nPromptLambda = 0;
  int nPromptAntiLambda = 0;
  int nPromptSigmaMinus = 0;
  int nPromptSigmaPlus = 0;
  int nPromptSigma0 = 0;
  int nPromptAntiKaon0 = 0;
  int nPromptAntiSigmaMinus = 0;
  int nPromptAntiSigma0 = 0;
  int nPromptAntiSigmaPlus = 0;
  int nPromptXsi0 = 0;
  int nPromptAntiXsi0 = 0;
  int nPromptXsiMinus = 0;
  int nPromptXsiPlus = 0;
  int nPromptOmegaMinus = 0;
  int nPromptOmegaPlus = 0;
  int nPromptOpticalPhoton = 0;
  int nPromptAlpha = 0;
  int nPromptDeuteron = 0;
  int nPromptTriton = 0;
  int nPromptLi7 = 0;
  int nPromptC10 = 0;
  int nPromptB11 = 0;
  int nPromptC12 = 0;
  int nPromptC13 = 0;
  int nPromptN13 = 0;
  int nPromptN14 = 0;
  int nPromptN15 = 0;
  int nPromptN16 = 0;
  int nPromptO16 = 0;
  int nPromptAl27 = 0;
  int nPromptFe54 = 0;
  int nPromptMn54 = 0;
  int nPromptMn55 = 0;
  int nPromptMn56 = 0;
  int nPromptFe56 = 0;
  int nPromptFe57 = 0;
  int nPromptFe58 = 0;
  int nPromptEu154 = 0;
  int nPromptGd158 = 0;
  int nPromptGd156 = 0;
  int nPromptGd157 = 0;
  int nPromptGd155 = 0;

  
  int nDelayedProton = 0;
  int nDelayedAntiProton = 0;
  int nDelayedElectron = 0;
  int nDelayedPositron = 0;
  int nDelayedElectronNeutrino = 0;
  int nDelayedAntiElectronNeutrino = 0;
  int nDelayedGamma = 0;
  int nDelayedNeutron = 0;
  int nDelayedAntiNeutron = 0;
  int nDelayedMuonPlus = 0;
  int nDelayedMuonMinus = 0;
  int nDelayedKaonlong = 0;
  int nDelayedPionPlus = 0;
  int nDelayedPionMinus = 0;
  int nDelayedKaonPlus = 0;
  int nDelayedKaonMinus = 0;
  int nDelayedKaonshort = 0;
  int nDelayedPion0 = 0;
  int nDelayedKaon0 = 0;
  int nDelayedMuonNeutrino = 0;
  int nDelayedAntiMuonNeutrino = 0;
  int nDelayedTauPlus = 0;
  int nDelayedTauMinus = 0;


  int nDelayedLambda = 0;
  int nDelayedAntiLambda = 0;
  int nDelayedSigmaMinus = 0;
  int nDelayedSigmaPlus = 0;
  int nDelayedSigma0 = 0;
  int nDelayedAntiKaon0 = 0;
  int nDelayedAntiSigmaMinus = 0;
  int nDelayedAntiSigma0 = 0;
  int nDelayedAntiSigmaPlus = 0;
  int nDelayedXsi0 = 0;
  int nDelayedAntiXsi0 = 0;
  int nDelayedXsiMinus = 0;
  int nDelayedXsiPlus = 0;
  int nDelayedOmegaMinus = 0;
  int nDelayedOmegaPlus = 0;
  int nDelayedOpticalPhoton = 0;
  int nDelayedAlpha = 0;
  int nDelayedDeuteron = 0;
  int nDelayedTriton = 0;
  int nDelayedLi7 = 0;
  int nDelayedC10 = 0;
  int nDelayedB11 = 0;
  int nDelayedC12 = 0;
  int nDelayedC13 = 0;
  int nDelayedN13 = 0;
  int nDelayedN14 = 0;
  int nDelayedN15 = 0;
  int nDelayedN16 = 0;
  int nDelayedO16 = 0;
  int nDelayedAl27 = 0;
  int nDelayedFe54 = 0;
  int nDelayedMn54 = 0;
  int nDelayedMn55 = 0;
  int nDelayedMn56 = 0;
  int nDelayedFe56 = 0;
  int nDelayedFe57 = 0;
  int nDelayedFe58 = 0;
  int nDelayedEu154 = 0;
  int nDelayedGd158 = 0;
  int nDelayedGd156 = 0;
  int nDelayedGd157 = 0;
  int nDelayedGd155 = 0;
  
  double fClusterChargeBalance;
  std::map<double, double> cluster_CB;




};


#endif
