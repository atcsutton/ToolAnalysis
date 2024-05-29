#ifndef PhaseIINeutronBG_H
#define PhaseIINeutronBG_H

#include <string>
#include <iostream>

#include "Tool.h"
#include "ADCPulse.h"
#include "Hit.h"
#include "Particle.h"
#include "BeamStatus.h"


#include "TApplication.h"
#include "TCanvas.h"
#include "TF1.h"
#include "TFile.h"
#include "TGraph.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TMath.h"
#include "TROOT.h"
#include "TSystem.h"
#include "TTree.h"


/**
 * \class PhaseIINeutronBG
 *
 *
 * $ Author: Julie He
 * $ Date: March 30, 2021
 * $ Contact: juhe@ucdavis.edu
 *
 * Adoptive parent: Emily Pottebaum
 * Contact: egpottebaum@gmail.edu
**/

class PhaseIINeutronBG: public Tool {


 public:

  PhaseIINeutronBG(); ///< Simple constructor
  bool Initialise(std::string configfile, DataModel &data); ///< Initialise Function for setting up Tool resources. @param configfile The path and name of the dynamic configuration file to read in. @param data A reference to the transient data class used to pass information between Tools.
  bool Execute(); ///< Execute function used to perform Tool purpose.
  bool Finalise(); ///< Finalise function used to clean up resources.

  void InitTree();
  void InitHist();
  void WriteHist();
  bool LoadTankClusterClassifiers(double cluster_time);

 private:

  // config variables
  int verbosity;
  std::string outpfile_prefix;
  int min_clusterPE;
  int max_clusterPE;
  int BeamOn;
  int BeamOff;
  Geometry *fGeo = nullptr;

  std::map<int, double> map_chankey2spe;
  std::map<double, std::vector<Hit>> *m_all_clusters = nullptr;
  std::map<double, Position> cluster_CP;
  std::map<double, double> cluster_CB;
  std::map<double, double> cluster_maxPEs;
  std::map<unsigned long, vector<Hit>> *tdcdata = nullptr;

  uint32_t trigword;
  int trigext;
  int goodnoncc;
  int goodnoncctrig;
  int pethresh;
  int PMTthresh;
  bool pot_ok;
  int cc_lt = 0;
  int nc_lt = 0;
  int beamoff_lt = 0;
  int woext_lt = 0; 
  uint64_t totalpot = 0;                                                                                                                     
  uint64_t beamoff_pot = 0;                                                                                                                  
  uint64_t woext_pot = 0;                                                                                                                   
  uint64_t cc_pot = 0;                                                                                                       
  uint64_t nc_pot = 0;

  // ROOT 
  TFile *p2nbg_root_outp = nullptr;

  // ROOT histograms
  TH1F *h_clusterCharge = nullptr;
  TH1F *h_beamFlag = nullptr;

  TH1F *h_beamoff_lt = nullptr;
  TH1F *h_woext_lt = nullptr;
  TH1F *h_cc_lt = nullptr;
  TH1F *h_nc_lt = nullptr;

  TH1F *h_totalpot = nullptr;

  TH1F *h_beamoff_pot = nullptr;
  TH1F *h_woext_pot = nullptr;
  TH1F *h_cc_pot = nullptr;
  TH1F *h_nc_pot = nullptr;

  TH1F *h_clusterTime = nullptr;
  TH1F *h_clusterTime_beamoff = nullptr;

  TH1F *h_clusterTime_beam = nullptr;
  TH1F *h_clusterTime_prompt = nullptr;
  TH1F *h_clusterTime_delayed = nullptr;
  TH1F *h_clusterTime_nonCCbeam = nullptr;
  TH1F *h_clusterTime_nonCCbeamoff = nullptr;

  TH1F *h_clusterTime_nonCCbeam_prompt = nullptr;
  TH1F *h_clusterTime_nonCCbeamoff_prompt = nullptr;

  TH1F *h_clusterTime_nonCCbeam_prompt_noVeto = nullptr;
  TH1F *h_clusterTime_nonCCbeam_prompt_noVetoMRD = nullptr;
  TH1F *h_clusterTime_nonCCbeam_delayed = nullptr;
  TH1F *h_clusterTime_nonCCbeamoff_delayed = nullptr;

  TH1F *h_clusterTime_nonCCbeam_delayed_noVeto = nullptr;
  TH1F *h_clusterTime_nonCCbeam_delayed_noVetoMRD = nullptr;
  TH1F *h_clusterTime_background_neutrons = nullptr;
  TH1F *h_clusterTime_hitPMT = nullptr;
  TH1F *h_clusterTime_delayed_hitPMT = nullptr;
  TH1F *h_clusterTime_prompt_hitPMT = nullptr;
  TH1F *h_clusterPE = nullptr;
  TH1F *h_clusterPE_beamoff = nullptr;

  TH1F *h_clusterPE_beam = nullptr;
  TH1F *h_clusterPE_prompt = nullptr;
  TH1F *h_clusterPE_delayed = nullptr;
  TH1F *h_clusterPE_nonCCbeam = nullptr;
  TH1F *h_clusterPE_nonCCbeamoff = nullptr;

  TH1F *h_clusterPE_nonCCbeam_prompt = nullptr;
  TH1F *h_clusterPE_nonCCbeamoff_prompt = nullptr;

  TH1F *h_clusterPE_nonCCbeam_prompt_noVeto = nullptr;
  TH1F *h_clusterPE_nonCCbeam_prompt_noVetoMRD = nullptr;
  TH1F *h_clusterPE_nonCCbeam_delayed = nullptr;
  TH1F *h_clusterPE_nonCCbeamoff_delayed = nullptr;

  TH1F *h_clusterPE_nonCCbeam_delayed_noVeto = nullptr;
  TH1F *h_clusterPE_nonCCbeam_delayed_noVetoMRD = nullptr;
  TH1F *h_clusterPE_background_neutrons = nullptr;
  TH1F *h_clusterCB_michel = nullptr;
  TH1F *h_clusterCB_michel_goodNonCC = nullptr;
  TH1F *h_clusterCB_afterpulse = nullptr;
  TH1F *h_clusterCB_afterpulse_goodNonCC = nullptr;
  TH1F *h_clusterCB_neutrons = nullptr;
  TH1F *h_clusterCB_neutrons_goodNonCC = nullptr;
  TH1F *h_clusterCB_muons = nullptr;
  TH1F *h_clusterCB_muons_goodNonCC = nullptr;

  TH1F *h_clusterCB_beamoff_michel = nullptr;
  TH1F *h_clusterCB_michel_goodNonCC_beamoff = nullptr;
  TH1F *h_clusterCB_afterpulse_beamoff = nullptr;
  TH1F *h_clusterCB_afterpulse_goodNonCC_beamoff = nullptr;
  TH1F *h_clusterCB_neutrons_beamoff = nullptr;
  TH1F *h_clusterCB_neutrons_goodNonCC_beamoff = nullptr;
  TH1F *h_clusterCB_muons_beamoff = nullptr;
  TH1F *h_clusterCB_muons_goodNonCC_beamoff = nullptr;


  // ROOT trees
  TTree *t_TankCluster = nullptr;
  TTree *t_Trigger = nullptr;
  uint32_t fRunNumber;
  uint32_t fSubrunNumber;
  uint32_t fRunType;
  uint64_t fStartTime;
  double fpot;
  int beamflag;
  uint32_t fEventNumber;
  uint64_t fEventTimeTank;
  int fClusterNumber;
  double fClusterCharge;
  double fClusterTime;
  double fClusterPE;
  double fClusterMaxPE;
  int fClusterHits;
  double fClusterChargeBalance;
  int fVetoHit;
  int fnMRDTracks;
  BeamStatus beamstat;
  std::vector<double> fHitPE;
  std::vector<unsigned long> fHitPMT;
  std::vector<double> fClusterHitPE;
  std::vector<int> fClusterPMT;

  int nevents = 0;

  // verbosity variables
  int v_error = 0;
  int v_warning = 1;
  int v_message = 2;
  int v_debug = 3;

};


#endif
