#include "CutsOptimizer.h"
#include <iostream>
#include <fstream>
#include <map>
#include <tuple>

CutsOptimizer::CutsOptimizer():Tool(){}


bool CutsOptimizer::Initialise(std::string configfile, DataModel &data){

  /////////////////// Useful header ///////////////////////
  if(configfile!="") m_variables.Initialise(configfile); // loading config file
  //m_variables.Print();

  m_data= &data; //assigning transient data pointer
  /////////////////////////////////////////////////////////////////

  bool gotVerbosity = m_variables.Get("verbosity", verbosity);
  if (!gotVerbosity){
    verbosity = 0;
    Log("CutsOptimizer: \"verbosity\" not set in the config, defaulting to 0", v_error, verbosity);
  }
  
  bool gotClusterMapName = m_variables.Get("ClusterMapName", fClusterMapName);
  if (!gotClusterMapName) {
    Log("CutsOptimizer: \"ClusterMapName\" not set in the config! Aborting!", v_error, verbosity);
    return false;
  }
  
  bool gotVertexMapName = m_variables.Get("VertexMapName", fVertexMapName);
  if (!gotVertexMapName) {
    Log("CutsOptimizer: \"VertexMapName\" not set in the config! Aborting!", v_error, verbosity);
    return false;
  }
  
  bool gotGeometry = m_data->Stores.at("ANNIEEvent")->Header->Get("AnnieGeometry", fGeo);
  if(!gotGeometry){
    Log("CutsOptimizer: Error retrieving Geometry from ANNIEEvent! Aborting!", v_error, verbosity);
    return false;
  }
  
  
  SetupTTree();
  return true;
}


bool CutsOptimizer::Execute(){
  if (!LoadFromStores())
    return false;

  nTotalTrueNeutronsPrompt = 0;
  nTotalTrueNeutronsDelayed = 0;
  nTotalNonTrueNeutronsPrompt = 0;
  nTotalNonTrueNeutronsDelayed = 0;
  
  for (auto& MCkey : *fMCParticles){
    int ParticlePDG = MCkey.GetPdgCode();
    int ParentPdg = MCkey.GetParentPdg();
    Position Positionvtx = MCkey.GetStopVertex();
    //Geo cut apply!!!!!!
    fTrueVtxX = Positionvtx.X();
    fTrueVtxY = Positionvtx.Y();
    fTrueVtxZ = Positionvtx.Z();
    IsInTank = fGeo->GetTankContained(Positionvtx);
    double clttime = MCkey.GetStopTime();
    // fTotalQ = fClusterTotalCharge->at(clttime);
    //    std::cout << fTotalQ << "<<Charge Total" << std::endl;
    
    if (IsInTank && ParticlePDG==2112 && ParentPdg ==0){
      if (clttime <= 2000.0){
	//	h_nTotalTrueNeutronsPrompt_CC->Fill(fTotalQ);
	//	h_nTotalTrueNeutronsPrompt_CT->Fill()
	nTotalTrueNeutronsPrompt++;}
      
      else if (clttime > 2000.0){nTotalTrueNeutronsDelayed++;}
    }
  }
  
  double minQ = 60.0;
  double maxQ = 180.0;
  double minCB = 0.2;
  double maxCB = 0.8;
  
  // Define the number of steps for each loop
  double chargeStepSize = 10.0;
  double cbStepSize = 0.05;
  
  for (chargeCut = minQ; chargeCut < maxQ; chargeCut += chargeStepSize) {
    for (cbCut = minCB; cbCut < maxCB; cbCut +=cbStepSize) {
      
      nAllSelectedClusterPrompt = 0;
      nSelectedTrueNeutronsPrompt = 0;
      nAllSelectedClusterDelayed = 0;
      nSelectedTrueNeutronsDelayed = 0;
      
      for (auto& clusterKey : *fClusterMap){
	double clusterTime = clusterKey.first;
	int bestPrtID = fClusterToBestParticleID->at(clusterTime);
	int bestPrtIdx = fMCParticleIndexMap->at(bestPrtID);
	fTotalQ = fClusterTotalCharge->at(clusterTime);
	fBestPDG = fClusterToBestParticlePDG->at(clusterTime);
	MCParticle bestPrt = fMCParticles->at(bestPrtIdx);
	Position pos = bestPrt.GetStopVertex();
	IsInTankMC = fGeo->GetTankContained(pos);
	
	bool good_class = this->LoadTankClusterClassifiers(clusterTime);
	if (!good_class) { Log("CutsOptimizer tool: NO cluster classifiers..", v_debug, verbosity); }
	
	/*	isPrompt = false;
		if (clusterTime <= 2000.0)
		isPrompt = true; */
	
	bool passCut = true;
	
	if ( fTotalQ > chargeCut ) passCut = false;
	
	if ( fClusterChargeBalance > cbCut ) passCut = false;
	
	if (fTotalQ < chargeCut && fClusterChargeBalance < cbCut) {
	  if (clusterTime <=2000.0){
	    nAllSelectedClusterPrompt++;
	    if (fBestPDG == 2112){
	      nSelectedTrueNeutronsPrompt++;
	    }
	  }
	  else if (clusterTime > 2000.0){
	    nAllSelectedClusterDelayed++;
            if (fBestPDG == 2112){
              nSelectedTrueNeutronsDelayed++;
            }
	  }
	}
	
	fOutTree->Fill();
      }
    }
  }
  return true;
}



bool CutsOptimizer::Finalise(){

  // fOutTree->Fill();
  fOutFile->cd();
  fOutTree->Write();
  //this->WriteHist();
  fOutFile->Close();
  return true;
}

void CutsOptimizer::SetupTTree()
{
  
  fOutFile = new TFile("CutsOptimizer.root", "RECREATE");
  fOutTree = new TTree("tree", "tree");
  fOutTree->Branch("nTotalTrueNeutronsPrompt",                 &nTotalTrueNeutronsPrompt);
  fOutTree->Branch("nTotalTrueNeutronsDelayed",                 &nTotalTrueNeutronsDelayed);
  fOutTree->Branch("nTotalNonTrueNeutronsPrompt",                 &nTotalNonTrueNeutronsPrompt);
  fOutTree->Branch("nTotalNonTrueNeutronsDelayed",                 &nTotalNonTrueNeutronsDelayed);
  fOutTree->Branch("nAllSelectedClusterPrompt", &nAllSelectedClusterPrompt);
  fOutTree->Branch("nSelectedTrueNeutronsPrompt", &nSelectedTrueNeutronsPrompt);
  fOutTree->Branch("nAllSelectedClusterDelayed", &nAllSelectedClusterDelayed);
  fOutTree->Branch("nSelectedTrueNeutronsDelayed", &nSelectedTrueNeutronsDelayed);
  fOutTree->Branch("isPrompt", &isPrompt);
  fOutTree->Branch("chargeCut", &chargeCut);
  fOutTree->Branch("cbCut", &cbCut);
  gROOT->cd();
  
}
bool CutsOptimizer::LoadFromStores()
{
  
  bool goodAnnieEvent = m_data->Stores.count("ANNIEEvent");
  if (!goodAnnieEvent) {
    logmessage = "CutsOptimizer: no ANNIEEvent store!";
    Log(logmessage, v_error, verbosity);
    return false;
  }
  
  bool goodClusterMap = m_data->Stores.at("ANNIEEvent")->Get(fClusterMapName, fClusterMap);
  if (!goodClusterMap) {
    logmessage = "CutsOptimizer: no " + fClusterMapName + " in the ANNIEEvent!";
    Log(logmessage, v_error, verbosity);
    return false;
  }
  
  bool goodClusterToBestParticleID = m_data->Stores.at("ANNIEEvent")->Get("ClusterToBestParticleID", fClusterToBestParticleID);
  if (!goodClusterToBestParticleID) {
    logmessage = "CutsOptimizer: no ClusterToBestParticleID in the ANNIEEvent!";
    Log(logmessage, v_error, verbosity);
    return false;
  }
  
  bool goodMCParticleIndexMap = m_data->Stores.at("ANNIEEvent")->Get("TrackId_to_MCParticleIndex", fMCParticleIndexMap);
  if (!goodMCParticleIndexMap) {
    logmessage = "CutsOptimizer: no TrackId_to_MCParticleIndex in the ANNIEEvent!";
    Log(logmessage, v_error, verbosity);
    return false;
  }
  
  bool goodClusterTotalCharge = m_data->Stores.at("ANNIEEvent")->Get("ClusterTotalCharge", fClusterTotalCharge);
  if (!goodClusterTotalCharge) {
    logmessage = "CutsOptimizer: no ClusterTotalCharge in the ANNIEEvent!";
    Log(logmessage, v_error, verbosity);
    return false;
  }
  
  bool goodClusterToBestParticlePDG = m_data->Stores.at("ANNIEEvent")->Get("ClusterToBestParticlePDG", fClusterToBestParticlePDG);
  if (!goodClusterToBestParticlePDG) {
    logmessage = "CutsOptimizer: no ClusterToBestParticlePDG in the ANNIEEvent!";
    Log(logmessage, v_error, verbosity);
    return false;
  }
  
  bool goodMCParticles = m_data->Stores.at("ANNIEEvent")->Get("MCParticles", fMCParticles);
  if (!goodMCParticles) {
    logmessage = "CutsOptimizer: no MCParticles in the ANNIEEvent!";
    Log(logmessage, v_error, verbosity);
    return false;
  }
  
  return true;
  
}  

bool CutsOptimizer::LoadTankClusterClassifiers(double clusterTime)
{
  
  bool got_ccb = m_data->Stores["ANNIEEvent"]->Get("ClusterChargeBalances", cluster_CB);
  
  bool good_class = got_ccb;
  if (!good_class) { Log("CutsOptimizer tool: One of the charge cluster classifiers is not available", v_debug, verbosity); }
  else
    {
      Log("CutsOptimizer tool:11 Setting fCluster variables to classifier parameters", v_debug, verbosity);
      fClusterChargeBalance = cluster_CB.at(clusterTime);
    }
  return good_class;
}
