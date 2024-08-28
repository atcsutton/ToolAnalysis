#include "SelectionEffnPurity.h"
#include <iostream>
#include <fstream>


SelectionEffnPurity::SelectionEffnPurity():Tool(){}


bool SelectionEffnPurity::Initialise(std::string configfile, DataModel &data){

  /////////////////// Useful header ///////////////////////
  if(configfile!="") m_variables.Initialise(configfile); // loading config file
  //m_variables.Print();

  m_data= &data; //assigning transient data pointer
  /////////////////////////////////////////////////////////////////
  // this->InitHist();
  //hahaha
  bool gotVerbosity = m_variables.Get("verbosity", verbosity);
  if (!gotVerbosity){
    verbosity = 0;
    Log("SelectionEffnPurity: \"verbosity\" not set in the config, defaulting to 0", v_error, verbosity);
  }

  bool gotClusterMapName = m_variables.Get("ClusterMapName", fClusterMapName);
  if (!gotClusterMapName) {
    Log("EvaluateVertex: \"ClusterMapName\" not set in the config! Aborting!", v_error, verbosity);
    return false;
  }
  
  bool gotVertexMapName = m_variables.Get("VertexMapName", fVertexMapName);
  if (!gotVertexMapName) {
    Log("EvaluateVertex: \"VertexMapName\" not set in the config! Aborting!", v_error, verbosity);
    return false;
  }

  bool gotGeometry = m_data->Stores.at("ANNIEEvent")->Header->Get("AnnieGeometry", fGeo);
  if(!gotGeometry){
    Log("SelectionEffnPurity:Error retrieving Geometry from ANNIEEvent! Aborting!", v_error, verbosity);
    return false;
  }

  SetupTTree();
  InitHist(6000.);
  return true;
}


bool SelectionEffnPurity::Execute(){
  if (!LoadFromStores())
    return false;

  Float_t NeutrinoEnergy=-2;
  bool isok;
  MCParticle neutrino;
  isok = m_data->Stores["ANNIEEvent"]->Get("NeutrinoParticle", neutrino);
  if (isok) NeutrinoEnergy = neutrino.GetStartEnergy();
  bool IsInTank;
  //  double mchits = 0;
  for (auto& MCkey : *fMCParticles){

    int ParticlePDG = MCkey.GetPdgCode();
    //    std::cout <<"MCParticle PDG code:" << ParticlePDG << std::endl;
    int ParentPdg = MCkey.GetParentPdg();
    //    std::cout << "MCParticle Parent PDG" << ParentPdg << std::endl;
    Position Positionvtx = MCkey.GetStopVertex();
    //Geo cut apply!!!!!!
    fTrueVtxX = Positionvtx.X();
    fTrueVtxY = Positionvtx.Y();
    fTrueVtxZ = Positionvtx.Z();
    IsInTank = fGeo->GetTankContained(Positionvtx);
    double clttime = MCkey.GetStopTime();

    std::cout << "Vertex X:-" << fTrueVtxX << "; Vertex Y:-" << fTrueVtxY << "; Vertex Z:-" << fTrueVtxZ << std::endl;
    
    if (ParticlePDG==2112 && ParentPdg ==0){
      nTotalTrueNeutronsWorld++;
      
      if (IsInTank){ //Selecting only Inside the tank events
	nTotalTrueNeutrons++;
	
	if (clttime < 2000.0){
	  //	  h_nTotalTrueNeutronsPromptNhits->Fill(MCNhits);
	  h_nTotalTrueNeutronsPromptPDG->Fill(ParticlePDG);
	  h_nTotalTrueNeutronsPromptCT->Fill(clttime);
	  h_nTotalTrueNeutronsPromptTVtxXY->Fill(fTrueVtxX, fTrueVtxY);
	  h_nTotalTrueNeutronsPromptTVtxYZ->Fill(fTrueVtxY, fTrueVtxZ);
	  h_nTotalTrueNeutronsPromptTVtxXZ->Fill(fTrueVtxX, fTrueVtxZ);
	  h_nTotalTrueNeutronsPromptNE ->Fill(NeutrinoEnergy);
	  nTotalTrueNeutronsPrompt++;
	}
	else if (clttime > 2000.0){
	  // h_nTotalTrueNeutronsDelayedNhits->Fill(MCNhits);
	  h_nTotalTrueNeutronsDelayedPDG->Fill(ParticlePDG);
	  h_nTotalTrueNeutronsDelayedCT->Fill(clttime);
	  h_nTotalTrueNeutronsDelayedTVtxXY->Fill(fTrueVtxX, fTrueVtxY);
	  h_nTotalTrueNeutronsDelayedTVtxYZ->Fill(fTrueVtxY, fTrueVtxZ);
	  h_nTotalTrueNeutronsDelayedTVtxXZ->Fill(fTrueVtxX, fTrueVtxZ);
	  h_nTotalTrueNeutronsDelayedNE->Fill(NeutrinoEnergy);
	  nTotalTrueNeutronsDelayed++;
	}
      }
    }
  }
  double hits = 0;
  for (auto& clusterKey : *fClusterMap){
    double clusterTime = clusterKey.first;
    int bestPrtID = fClusterToBestParticleID->at(clusterTime);
    int bestPrtIdx = fMCParticleIndexMap->at(bestPrtID);
    fTotalQ = fClusterTotalCharge->at(clusterTime);
    fBestPDG = fClusterToBestParticlePDG->at(clusterTime);
    
    const std::vector<MCHit>& hits = clusterKey.second;
    size_t Nhits = hits.size();
    
    //    std::cout << "Geo contained" << IsInTank << std::endl;
    //    std::cout << "Vertex X:-" << fTrueVtxX << std::endl; 
    
    bool good_class = this->LoadTankClusterClassifiers(clusterTime);
    if (!good_class) { Log("PhaseIINeutronBG tool: NO cluster classifiers..", v_debug, verbosity); }
    
    
    //looping over all the particles from ClusterMap
    if (IsInTank){
      if (fBestPDG != 0){
	nAllSelectedClustersWorld++;
	if (clusterTime < 2000.0){ //Prompt window && neutron selection cuts
	  if (fClusterChargeBalance < 0.4 && fTotalQ < 120 && fClusterChargeBalance < 0.5 - fTotalQ / 300){
	    h_nAllSelectedClustersPromptNhits->Fill(Nhits);
	    h_nAllSelectedClustersPromptPDG->Fill(fBestPDG);
	    h_nAllSelectedClustersPromptCT->Fill(clusterTime);
	    h_nAllSelectedClustersPromptTVtxXY->Fill(fTrueVtxX, fTrueVtxY);
	    h_nAllSelectedClustersPromptTVtxYZ->Fill(fTrueVtxY, fTrueVtxZ);
	    h_nAllSelectedClustersPromptTVtxXZ->Fill(fTrueVtxX, fTrueVtxZ);
	    h_nAllSelectedClustersPromptNE->Fill(NeutrinoEnergy);
	    nAllSelectedClustersPrompt++; 
	    
	    
	    if (fBestPDG == 2112){
	      h_nSelectedTrueNeutronsPromptNhits->Fill(Nhits);
	      h_nSelectedTrueNeutronsPromptPDG->Fill(fBestPDG);
	      h_nSelectedTrueNeutronsPromptCT->Fill(clusterTime);
	      //	      h_nSelectedTrueNeutronsPromptTVtxXY->Fill(fTrueVtxX, fTrueVtxY);
	      h_nSelectedTrueNeutronsPromptTVtxYZ->Fill(fTrueVtxY, fTrueVtxZ);
	      h_nSelectedTrueNeutronsPromptTVtxXZ->Fill(fTrueVtxX, fTrueVtxZ);
	      h_nSelectedTrueNeutronsPromptNE->Fill(NeutrinoEnergy);
	      nSelectedTrueNeutronsPrompt++;
	    }   	
	  }	
	}
	else if (clusterTime > 2000.0){ //Delayed Window 
	  if (fClusterChargeBalance < 0.4 && fTotalQ < 120 && fClusterChargeBalance < 0.5 - fTotalQ / 300){
	    h_nAllSelectedClustersDelayedNhits->Fill(Nhits);
	    h_nAllSelectedClustersDelayedPDG->Fill(fBestPDG);
	    h_nAllSelectedClustersDelayedCT->Fill(clusterTime);
	    h_nAllSelectedClustersDelayedTVtxXY->Fill(fTrueVtxX, fTrueVtxY);
	    h_nAllSelectedClustersDelayedTVtxYZ->Fill(fTrueVtxY, fTrueVtxZ);
	    h_nAllSelectedClustersDelayedTVtxXZ->Fill(fTrueVtxX, fTrueVtxZ);
	    h_nAllSelectedClustersDelayedNE->Fill(NeutrinoEnergy);
	    nAllSelectedClustersDelayed++;
	    
	    if (fBestPDG == 2112){
	      h_nSelectedTrueNeutronsDelayedNhits->Fill(Nhits);
	      h_nSelectedTrueNeutronsDelayedPDG->Fill(fBestPDG);
	      h_nSelectedTrueNeutronsDelayedCT->Fill(clusterTime);
	      h_nSelectedTrueNeutronsDelayedTVtxXY->Fill(fTrueVtxX, fTrueVtxY);
	      h_nSelectedTrueNeutronsDelayedTVtxYZ->Fill(fTrueVtxY, fTrueVtxZ);
	      h_nSelectedTrueNeutronsDelayedTVtxXZ->Fill(fTrueVtxX, fTrueVtxZ);
	      h_nSelectedTrueNeutronsDelayedNE->Fill(NeutrinoEnergy);
	      nSelectedTrueNeutronsDelayed++;
	    }
	  }
	}  
      } 
    }
  }
  
  return true;
}


bool SelectionEffnPurity::Finalise(){
  std::cout << "Total True Neutrons before geo cut:-" << nTotalTrueNeutronsWorld << std::endl;
  std::cout << "Total True Neutrons:-" << nTotalTrueNeutrons << std::endl;
  std::cout << "Total True Neutrons Prompt:-" << nTotalTrueNeutronsPrompt << std::endl;
  std::cout << "Total True Neutrons Delayed:-" << nTotalTrueNeutronsDelayed << std::endl;
  std::cout << "All Selected cluster before charge/balance cut:-" << nAllSelectedClustersWorld << std::endl;
  std::cout << "All Selected Cluster Prompt:-" << nAllSelectedClustersPrompt << std::endl; 
  std::cout << "Selected True Neutrons Prompt:-" << nSelectedTrueNeutronsPrompt << std::endl;
  std::cout << "All Selected Cluster Delayed:-" << nAllSelectedClustersDelayed << std::endl;
  std::cout << "Selected True Neutrons Delayed:-" << nSelectedTrueNeutronsDelayed << std::endl;
  
  fOutTree->Fill();
  fOutFile->cd();
  fOutTree->Write();
  //  this->WriteHist();
  fOutFile->Close();
  return true;
}

void SelectionEffnPurity::SetupTTree()
{
  
  fOutFile = new TFile("SelectionEffnPurity.root", "RECREATE");
  fOutTree = new TTree("tree", "tree");
  
  fOutTree->Branch("nTotalTrueNeutronsWorld",            &nTotalTrueNeutronsWorld);
  fOutTree->Branch("nTotalTrueNeutrons",                 &nTotalTrueNeutrons);
  fOutTree->Branch("nTotalTrueNeutronsPrompt",           &nTotalTrueNeutronsPrompt);
  fOutTree->Branch("nTotalTrueNeutronsDelayed",          &nTotalTrueNeutronsDelayed);
  fOutTree->Branch("nAllSelectedClustersWorld",          &nAllSelectedClustersWorld);
  fOutTree->Branch("nAllSelectedClustersPrompt",         &nAllSelectedClustersPrompt);
  fOutTree->Branch("nSelectedTrueNeutronsPrompt",        &nSelectedTrueNeutronsPrompt);
  fOutTree->Branch("nAllSelectedClustersDelayed",        &nAllSelectedClustersDelayed);
  fOutTree->Branch("nSelectedTrueNeutronsDelayed",       &nSelectedTrueNeutronsDelayed);
  
  gROOT->cd();
  
}

void SelectionEffnPurity::InitHist(double max)
{
  
  //  fOutFile->cd();
  h_nTotalTrueNeutronsPromptNE = new TH1F("h_nTotalTrueNeutronsPromptNE", "h_nTotalTrueNeutronsPromptNE", 12, 0, max);
  h_nTotalTrueNeutronsDelayedNE = new TH1F("h_nTotalTrueNeutronsDelayedNE", "h_nTotalTrueNeutronsDelayedNE", 12, 0 ,max);
  h_nAllSelectedClustersPromptNE =new TH1F("h_nAllSelectedClustersPromptNE", "h_nAllSelectedClustersPromptNE", 12, 0, max);
  h_nSelectedTrueNeutronsPromptNE = new TH1F("h_nSelectedTrueNeutronsPromptNE", "h_nSelectedTrueNeutronsPromptNE", 12, 0, max);
  h_nAllSelectedClustersDelayedNE = new TH1F("h_nAllSelectedClustersDelayedNE", "h_nAllSelectedClustersDelayedNE", 12, 0, max);
  h_nSelectedTrueNeutronsDelayedNE = new TH1F("h_nSelectedTrueNeutronsDelayedNE", "h_nSelectedTrueNeutronsDelayedNE", 12, 0, max);
  
  h_nTotalTrueNeutronsPromptTVtxXZ = new TH2F("h_nTotalTrueNeutronsPromptTVtxXZ", "h_nTotalTrueNeutronsPromptTVtxXZ", 40, -10, 10, 40, -10, 10);
  h_nTotalTrueNeutronsDelayedTVtxXZ = new TH2F("h_nTotalTrueNeutronsDelayedTVtxXZ", "h_nTotalTrueNeutronsDelayedTVtxXZ", 40, -10, 10, 40, -10, 10);
  h_nAllSelectedClustersPromptTVtxXZ = new TH2F("h_nAllSelectedClustersPromptTVtxXZ", "h_nAllSelectedClustersPromptTVtxXZ", 40, -10, 10, 40, -10, 10);
  h_nSelectedTrueNeutronsPromptTVtxXZ = new TH2F("h_nSelectedTrueNeutronsPromptTVtxXZ", "h_nSelectedTrueNeutronsPromptTVtxXZ", 40, -10, 10, 40, -10, 10);
  h_nAllSelectedClustersDelayedTVtxXZ = new TH2F("h_nAllSelectedClustersDelayedTVtxXZ", "h_nAllSelectedClustersDelayedTVtxXZ", 40, -10, 10, 40, -10, 10);
  h_nSelectedTrueNeutronsDelayedTVtxXZ = new TH2F("h_nSelectedTrueNeutronsDelayedTVtxXZ", "h_nSelectedTrueNeutronsDelayedTVtxXZ", 40, -10, 10, 40, -10, 10);

  h_nTotalTrueNeutronsPromptTVtxXY = new TH2F("h_nTotalTrueNeutronsPromptTVtxXY", "h_nTotalTrueNeutronsPromptTVtxXY", 40, -10, 10, 40, -10, 10);
  h_nTotalTrueNeutronsDelayedTVtxXY = new TH2F("h_nTotalTrueNeutronsDelayedTVtxXY", "h_nTotalTrueNeutronsDelayedTVtxXY", 40, -10, 10, 40, -10, 10);
  h_nAllSelectedClustersPromptTVtxXY = new TH2F("h_nAllSelectedClustersPromptTVtxXY", "h_nAllSelectedClustersPromptTVtxXY", 40, -10, 10, 40, -10, 10);
  h_nAllSelectedClustersDelayedTVtxXY = new TH2F("h_nAllSelectedClustersDelayedTVtxXY", "h_nAllSelectedClustersDelayedTVtxXY", 40, -10, 10, 40, -10, 10);
  h_nSelectedTrueNeutronsDelayedTVtxXY = new TH2F("h_nSelectedTrueNeutronsDelayedTVtxXY", "h_nSelectedTrueNeutronsDelayedTVtxXY", 40, -10, 10, 40, -10, 10);

  h_nTotalTrueNeutronsPromptTVtxYZ = new TH2F("h_nTotalTrueNeutronsPromptTVtxYZ", "h_nTotalTrueNeutronsPromptTVtxYZ", 40, -10, 10, 40, -10, 10);
  h_nTotalTrueNeutronsDelayedTVtxYZ = new TH2F("h_nTotalTrueNeutronsDelayedTVtxYZ", "h_nTotalTrueNeutronsDelayedTVtxYZ", 40, -10, 10, 40, -10, 10);
  h_nAllSelectedClustersPromptTVtxYZ = new TH2F("h_nAllSelectedClustersPromptTVtxYZ", "h_nAllSelectedClustersPromptTVtxYZ", 40, -10, 5, 40, -10, 10);
  h_nSelectedTrueNeutronsPromptTVtxYZ = new TH2F("h_nSelectedTrueNeutronsPromptTVtxYZ", "h_nSelectedTrueNeutronsPromptTVtxYZ", 40, -10, 5, 40, -10, 10);
  h_nAllSelectedClustersDelayedTVtxYZ = new TH2F("h_nAllSelectedClustersDelayedTVtxYZ", "h_nAllSelectedClustersDelayedTVtxYZ",40, -10, 10, 40, -10, 10);
  h_nSelectedTrueNeutronsDelayedTVtxYZ = new TH2F("h_nSelectedTrueNeutronsDelayedTVtxYZ", "h_nSelectedTrueNeutronsDelayedTVtxYZ", 40, -10, 10, 40, -10, 10);
  
  h_nTotalTrueNeutronsPromptCT = new TH1F("h_nTotalTrueNeutronsPromptCT", "h_nTotalTrueNeutronsPromptCT", 100, 0, 71000);
  h_nTotalTrueNeutronsDelayedCT = new TH1F("h_nTotalTrueNeutronsDelayedCT", "h_nTotalTrueNeutronsDelayedCT", 100, 0, 71000);
  h_nAllSelectedClustersPromptCT = new TH1F("h_nAllSelectedClustersPromptCT", "h_nAllSelectedClustersPromptCT", 100, 0, 71000);
  h_nSelectedTrueNeutronsPromptCT = new TH1F("h_nSelectedTrueNeutronsPromptCT", "h_nSelectedTrueNeutronsPromptCT", 100, 0, 71000);
  h_nAllSelectedClustersDelayedCT = new TH1F("h_nAllSelectedClustersDelayedCT", "h_nAllSelectedClustersDelayedCT", 100, 0, 71000);
  h_nSelectedTrueNeutronsDelayedCT = new TH1F("h_nSelectedTrueNeutronsDelayedCT", "h_nSelectedTrueNeutronsDelayedCT", 100, 0, 71000);

  h_nTotalTrueNeutronsPromptPDG = new TH1F("h_nTotalTrueNeutronsPromptPDG", "h_nTotalTrueNeutronsPromptPDG", 50, -3500, 3500);
  h_nTotalTrueNeutronsDelayedPDG = new TH1F("h_nTotalTrueNeutronsDelayedPDG", "h_nTotalTrueNeutronsDelayedPDG", 50, -3500, 3500);
  h_nAllSelectedClustersPromptPDG = new TH1F("h_nAllSelectedClustersPromptPDG", "h_nAllSelectedClustersPromptPDG", 50, -3500, 3500);
  h_nSelectedTrueNeutronsPromptPDG = new TH1F("h_nSelectedTrueNeutronsPromptPDG", "h_nSelectedTrueNeutronsPromptPDG", 50, -3500, 3500);
  h_nAllSelectedClustersDelayedPDG = new TH1F("h_nAllSelectedClustersDelayedPDG", "h_nAllSelectedClustersDelayedPDG", 50, -3500, 3500);
  h_nSelectedTrueNeutronsDelayedPDG = new TH1F("h_nSelectedTrueNeutronsDelayedPDG", "h_nSelectedTrueNeutronsDelayedPDG", 50, -3500, 3500);


  h_nTotalTrueNeutronsPromptNhits = new TH1F("h_nTotalTrueNeutronsPromptNhits", "h_nTotalTrueNeutronsPromptNhits", 50, 0, 200);
  h_nTotalTrueNeutronsDelayedNhits = new TH1F("h_nTotalTrueNeutronsDelayedNhits", "h_nTotalTrueNeutronsDelayedNhits", 50, 0, 200);
  h_nAllSelectedClustersPromptNhits = new TH1F("h_nAllSelectedClustersPromptNhits", "h_nAllSelectedClustersPromptNhits", 50, 0, 200);
  h_nSelectedTrueNeutronsPromptNhits = new TH1F("h_nSelectedTrueNeutronsPromptNhits", "h_nSelectedTrueNeutronsPromptNhits", 50, 0, 200);
  h_nAllSelectedClustersDelayedNhits = new TH1F("h_nAllSelectedClustersDelayedNhits", "h_nAllSelectedClustersDelayedNhits", 50, 0, 200);
  h_nSelectedTrueNeutronsDelayedNhits = new TH1F("h_nSelectedTrueNeutronsDelayedNhits", "h_nSelectedTrueNeutronsDelayedNhits", 50, 0, 200);

  gROOT->cd();

}

void SelectionEffnPurity::WriteHist()
{

  fOutFile->cd();
  TDirectory *dir_allhist = fOutFile->mkdir("Histograms");
  dir_allhist->cd();

  h_nTotalTrueNeutronsPromptNE->Write();
  h_nTotalTrueNeutronsDelayedNE->Write();
  h_nAllSelectedClustersPromptNE->Write();
  h_nSelectedTrueNeutronsPromptNE->Write();
  h_nAllSelectedClustersDelayedNE->Write();
  h_nSelectedTrueNeutronsDelayedNE->Write();

  h_nTotalTrueNeutronsPromptTVtxXZ->Write();
  h_nTotalTrueNeutronsDelayedTVtxXZ->Write();
  h_nAllSelectedClustersPromptTVtxXZ->Write();
  h_nSelectedTrueNeutronsPromptTVtxXZ->Write();
  h_nAllSelectedClustersDelayedTVtxXZ->Write();
  h_nSelectedTrueNeutronsDelayedTVtxXZ->Write();

  h_nTotalTrueNeutronsPromptTVtxXY->Write();
  h_nTotalTrueNeutronsDelayedTVtxXY->Write();
  h_nAllSelectedClustersPromptTVtxXY->Write();
  h_nSelectedTrueNeutronsPromptTVtxXY->Write();
  h_nAllSelectedClustersDelayedTVtxXY->Write();
  h_nSelectedTrueNeutronsDelayedTVtxXY->Write();

  h_nTotalTrueNeutronsPromptTVtxYZ->Write();
  h_nTotalTrueNeutronsDelayedTVtxYZ->Write();
  h_nAllSelectedClustersPromptTVtxYZ->Write();
  h_nSelectedTrueNeutronsPromptTVtxYZ->Write();
  h_nAllSelectedClustersDelayedTVtxYZ->Write();
  h_nSelectedTrueNeutronsDelayedTVtxYZ->Write();
  
  h_nTotalTrueNeutronsPromptCT->Write();
  h_nTotalTrueNeutronsDelayedCT->Write();
  h_nAllSelectedClustersPromptCT->Write();
  h_nSelectedTrueNeutronsPromptCT->Write();
  h_nAllSelectedClustersDelayedCT->Write();
  h_nSelectedTrueNeutronsDelayedCT->Write();

  h_nTotalTrueNeutronsPromptPDG->Write();
  h_nTotalTrueNeutronsDelayedPDG->Write();
  h_nAllSelectedClustersPromptPDG->Write();
  h_nSelectedTrueNeutronsPromptPDG->Write();
  h_nAllSelectedClustersDelayedPDG->Write();
  h_nSelectedTrueNeutronsDelayedPDG->Write();

  h_nTotalTrueNeutronsPromptNhits->Write();
  h_nTotalTrueNeutronsDelayedNhits->Write();
  h_nAllSelectedClustersPromptNhits->Write();
  h_nSelectedTrueNeutronsPromptNhits->Write();
  h_nAllSelectedClustersDelayedNhits->Write();
  h_nSelectedTrueNeutronsDelayedNhits->Write();
  
  gROOT->cd();

}

bool SelectionEffnPurity::LoadFromStores()
{
  bool goodAnnieEvent = m_data->Stores.count("ANNIEEvent");
  if (!goodAnnieEvent) {
    logmessage = "SelectionEffnPurity:0. no ANNIEEvent store!";
    Log(logmessage, v_error, verbosity);
    return false;
  }
  
  bool goodClusterMap = m_data->Stores.at("ANNIEEvent")->Get(fClusterMapName, fClusterMap);
  if (!goodClusterMap) {
    logmessage = "SelectionEffnPurity: 1.no " + fClusterMapName + " in the ANNIEEvent!";
    Log(logmessage, v_error, verbosity);
    return false;
  }
    
  bool goodMCParticles = m_data->Stores.at("ANNIEEvent")->Get("MCParticles", fMCParticles);
  if (!goodMCParticles) {
    logmessage = "SelectionEffnPurity:3. no MCParticles in the ANNIEEvent!";
    Log(logmessage, v_error, verbosity);
    return false;
  }
  
  bool goodMCParticleIndexMap = m_data->Stores.at("ANNIEEvent")->Get("TrackId_to_MCParticleIndex", fMCParticleIndexMap);
  if (!goodMCParticleIndexMap) {
    logmessage = "SelectionEffnPurity: 4.no TrackId_to_MCParticleIndex in the ANNIEEvent!";
    Log(logmessage, v_error, verbosity);
    return false;
  }
  
  bool goodClusterToBestParticleID = m_data->Stores.at("ANNIEEvent")->Get("ClusterToBestParticleID", fClusterToBestParticleID);
  if (!goodClusterToBestParticleID) {
    logmessage = "SelectionEffnPurity: 5.no ClusterToBestParticleID in the ANNIEEvent!";
    Log(logmessage, v_error, verbosity);
    return false;
  }
  
  bool goodClusterToBestParticlePDG = m_data->Stores.at("ANNIEEvent")->Get("ClusterToBestParticlePDG", fClusterToBestParticlePDG);
  if (!goodClusterToBestParticlePDG) {
    logmessage = "SelectionEffnPurity: 6.no ClusterToBestParticlePDG in the ANNIEEvent!";
    Log(logmessage, v_error, verbosity);
    return false;
  }
  
  bool goodClusterEfficiency = m_data->Stores.at("ANNIEEvent")->Get("ClusterEfficiency", fClusterEfficiency);
  if (!goodClusterEfficiency) {
    logmessage = "SelectionEffnPurity: 7.no ClusterEfficiency in the ANNIEEvent!";
    Log(logmessage, v_error, verbosity);
    return false;
  }
  
  bool goodClusterPurity = m_data->Stores.at("ANNIEEvent")->Get("ClusterPurity", fClusterPurity);
  if (!goodClusterPurity) {
    logmessage = "SelectionEffnPurity: 8.no ClusterPurity in the ANNIEEvent!";
    Log(logmessage, v_error, verbosity);
    return false;
  }
  
  bool goodClusterTotalCharge = m_data->Stores.at("ANNIEEvent")->Get("ClusterTotalCharge", fClusterTotalCharge);
  if (!goodClusterTotalCharge) {
    logmessage = "SelectionEffnPurity: 9.no ClusterTotalCharge in the ANNIEEvent!";
    Log(logmessage, v_error, verbosity);
    return false;
  }
  
  bool goodClusterNeutronCharge = m_data->Stores.at("ANNIEEvent")->Get("ClusterNeutronCharge", fClusterNeutronCharge);
  if (!goodClusterNeutronCharge) {
    logmessage = "SelectionEffnPurity: 10.no ClusterNeutronCharge in the ANNIEEvent!";
    Log(logmessage, v_error, verbosity);
    return false;
    }
  
  return true;
}

bool SelectionEffnPurity::LoadTankClusterClassifiers(double clusterTime)
{

  bool got_ccb = m_data->Stores["ANNIEEvent"]->Get("ClusterChargeBalances", cluster_CB);

  bool good_class = got_ccb;
  if (!good_class) { Log("SelectionEffnPurity tool: One of the charge cluster classifiers is not available", v_debug, verbosity); }
  else
  {
    Log("SelectionEffnPurity tool:11 Setting fCluster variables to classifier parameters", v_debug, verbosity);
    fClusterChargeBalance = cluster_CB.at(clusterTime);
  }
  return good_class;
  }
