#include "SelectionEffnPurity.h"
#include <iostream>
#include <fstream>
#include <map>

SelectionEffnPurity::SelectionEffnPurity():Tool(){}


bool SelectionEffnPurity::Initialise(std::string configfile, DataModel &data){

  /////////////////// Useful header ///////////////////////
  if(configfile!="") m_variables.Initialise(configfile); // loading config file
  //m_variables.Print();

  m_data= &data; //assigning transient data pointer
  /////////////////////////////////////////////////////////////////

  SetupTTree();

  bool gotVerbosity = m_variables.Get("verbosity", verbosity);
  if (!gotVerbosity){
    verbosity = 0;
    Log("1. SelectionEffnPurity: \"verbosity\" not set in the config, defaulting to 0", v_error, verbosity);
  }

  bool gotClusterMapName = m_variables.Get("ClusterMapName", fClusterMapName);
  if (!gotClusterMapName) {
    Log("2. SelectionEffnPurity: \"ClusterMapName\" not set in the config! Aborting!", v_error, verbosity);
    return false;
  }

  bool gotGeometry = m_data->Stores.at("ANNIEEvent")->Header->Get("AnnieGeometry", fGeo);
  if(!gotGeometry){
    Log("4. SelectionEffnPurity:Error retrieving Geometry from ANNIEEvent! Aborting!", v_error, verbosity);
    return false;
  }


  InitHist(6000.);
  return true;
}


bool SelectionEffnPurity::Execute(){
  bool skip = false;
  bool goodSkipStatus = m_data->Stores.at("ANNIEEvent")->Get("SkipExecute", skip);
  if (goodSkipStatus && skip) {
    logmessage = "BackTracker: An upstream tool told me to skip this event.";
    Log(logmessage, v_warning, verbosity);
    return true;
  }
  
  if (!LoadFromStores())
    return false;
  
  Float_t NeutrinoEnergy=-2;
  bool isok;
  MCParticle neutrino;
  isok = m_data->Stores["ANNIEEvent"]->Get("NeutrinoParticle", neutrino);
  if (isok) NeutrinoEnergy = neutrino.GetStartEnergy();
  bool IsInTank,IsInTankMC ;
  //  double mchits = 0;
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
    
    //    std::cout << "Vertex X:-" << fTrueVtxX << "; Vertex Y:-" << fTrueVtxY << "; Vertex Z:-" << fTrueVtxZ << std::endl;
    
    if (ParticlePDG==2112 && ParentPdg ==0){
      nTotalTrueNeutronsWorld++;
      
      if (IsInTank){ //Selecting only Inside the tank events
	nTotalTrueNeutrons++;
	
	if (clttime <= 2000.0){
	  //	  h_nTotalTrueNeutronsPromptNhits->Fill(MCNhits);
	  h_nTotalTrueNeutronsPromptPDG->Fill(ParticlePDG);
	  h_nTotalTrueNeutronsPromptCT->Fill(clttime);
	  h_nTotalTrueNeutronsPromptTVtxXY->Fill(fTrueVtxX, fTrueVtxY);
	  h_nTotalTrueNeutronsPromptTVtxYZ->Fill(fTrueVtxY, fTrueVtxZ);
	  h_nTotalTrueNeutronsPromptTVtxXZ->Fill(fTrueVtxX, fTrueVtxZ);
	  h_nTotalTrueNeutronsPromptNE ->Fill(NeutrinoEnergy);
	  //	  h_nTotalTrueNeutronsPromptCC->Fill(fTotalQ);
	  //	  std::cout << "Total True Neutron PROMPT" << std::endl;

	  //combination of both Delayed and Prompt                                                                                                                                                           
          h_nTotalTrueNeutronsTVtxXZ->Fill(fTrueVtxX, fTrueVtxZ);
	  nTotalTrueNeutronsPrompt++;
	  nSumingTotalTrueNeutron++;
	}
	else if (clttime > 2000.0){
	  // h_nTotalTrueNeutronsDelayedNhits->Fill(MCNhits);
	  h_nTotalTrueNeutronsDelayedPDG->Fill(ParticlePDG);
	  h_nTotalTrueNeutronsDelayedCT->Fill(clttime);
	  h_nTotalTrueNeutronsDelayedTVtxXY->Fill(fTrueVtxX, fTrueVtxY);
	  h_nTotalTrueNeutronsDelayedTVtxYZ->Fill(fTrueVtxY, fTrueVtxZ);
	  h_nTotalTrueNeutronsDelayedTVtxXZ->Fill(fTrueVtxX, fTrueVtxZ);
	  h_nTotalTrueNeutronsDelayedNE->Fill(NeutrinoEnergy);

	  //combination of both Delayed and Prompt
	  h_nTotalTrueNeutronsTVtxXZ->Fill(fTrueVtxX, fTrueVtxZ);
	  // std::cout << "Total True Neutron DELAYED" << std::endl;
	  nTotalTrueNeutronsDelayed++;
	  nSumingTotalTrueNeutron++;
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
    MCParticle bestPrt = fMCParticles->at(bestPrtIdx);
    Position pos = bestPrt.GetStopVertex();
    IsInTankMC = fGeo->GetTankContained(pos);
    fMCX = pos.X();
    fMCY = pos.Y();
    fMCZ = pos.Z();
    fClusterChargeBalance = ClusterChargeBalances.at(clusterTime);
    
    const std::vector<MCHit>& hits = clusterKey.second;
    size_t Nhits = hits.size();

    //    std::map<int, int> ParticleCountsDelayed;
    //    std::map<int, int> ParticleCountsPrompt;
    
    //looping over all the particles from ClusterMap
    if (IsInTankMC){
      if (fBestPDG != 0){
	nAllSelectedClustersWorld++;
	
	if (clusterTime <= 2000.0){ //Prompt window && neutron selection cuts
	  if (fClusterChargeBalance < 0.4 && fTotalQ < 120 && fClusterChargeBalance < 0.5 - fTotalQ / 300){
	    h_nAllSelectedClustersPromptNhits->Fill(Nhits);
	    h_nAllSelectedClustersPromptPDG->Fill(fBestPDG);
	    h_nAllSelectedClustersPromptCT->Fill(clusterTime);
	    h_nAllSelectedClustersPromptTVtxXY->Fill(fMCX, fMCY);
	    h_nAllSelectedClustersPromptTVtxYZ->Fill(fMCY, fMCZ);
	    h_nAllSelectedClustersPromptTVtxXZ->Fill(fMCX, fMCZ);
	    h_nAllSelectedClustersPromptNE->Fill(NeutrinoEnergy);
	    nAllSelectedClustersPrompt++; 
	    
	    if (fBestPDG == 2112){
	      h_nSelectedTrueNeutronsPromptNhits->Fill(Nhits);
	      h_nSelectedTrueNeutronsPromptPDG->Fill(fBestPDG);
	      h_nSelectedTrueNeutronsPromptCT->Fill(clusterTime);
	      h_nSelectedTrueNeutronsPromptTVtxXY->Fill(fMCX, fMCY);
	      h_nSelectedTrueNeutronsPromptTVtxYZ->Fill(fMCY, fMCZ);
	      h_nSelectedTrueNeutronsPromptTVtxXZ->Fill(fMCX, fMCZ);
	      h_nSelectedTrueNeutronsPromptNE->Fill(NeutrinoEnergy);
	      //  std::cout<< "Selected true neutron PROMPT" << std::endl;

	      //Combination of both delayed and prompt
	      h_nSelectedTrueNeutronsTVtxXZ->Fill(fMCX, fMCZ);
	      nSelectedTrueNeutronsPrompt++;
	      nSumingAllTrueNeutron++;
	    }
	
	    ParticleCountsPrompt[fBestPDG]++;
	  }	
	}
	else if (clusterTime > 2000.0){
	  //Delayed Window
	  if (fClusterChargeBalance < 0.4 && fTotalQ < 120 && fClusterChargeBalance < 0.5 - fTotalQ / 300){
	    h_nAllSelectedClustersDelayedNhits->Fill(Nhits);
	    h_nAllSelectedClustersDelayedPDG->Fill(fBestPDG);
	    h_nAllSelectedClustersDelayedCT->Fill(clusterTime);
	    h_nAllSelectedClustersDelayedTVtxXY->Fill(fMCX, fMCY);
	    h_nAllSelectedClustersDelayedTVtxYZ->Fill(fMCY, fMCZ);
	    h_nAllSelectedClustersDelayedTVtxXZ->Fill(fMCX, fMCZ);
	    h_nAllSelectedClustersDelayedNE->Fill(NeutrinoEnergy);
	    nAllSelectedClustersDelayed++;
	    
	    if (fBestPDG == 2112){
	      h_nSelectedTrueNeutronsDelayedNhits->Fill(Nhits);
	      h_nSelectedTrueNeutronsDelayedPDG->Fill(fBestPDG);
	      h_nSelectedTrueNeutronsDelayedCT->Fill(clusterTime);
	      h_nSelectedTrueNeutronsDelayedTVtxXY->Fill(fMCX, fMCY);
	      h_nSelectedTrueNeutronsDelayedTVtxYZ->Fill(fMCY, fMCZ);
	      h_nSelectedTrueNeutronsDelayedTVtxXZ->Fill(fMCX, fMCZ);
	      h_nSelectedTrueNeutronsDelayedNE->Fill(NeutrinoEnergy);
	      // std::cout << "Selected true neutrons DELAYED" << std::endl;
	      
	      //Combination of both delayed and prompt                                                                                                                                                       
	      h_nSelectedTrueNeutronsTVtxXZ->Fill(fMCX, fMCZ);
	      nSelectedTrueNeutronsDelayed++;
	      nSumingAllTrueNeutron++;
	    }
	    
	    ParticleCountsDelayed[fBestPDG]++;
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

  std::cout << "Sum of all Total True Neutrons:-" << nSumingTotalTrueNeutron <<std::endl;
  std::cout << "Sum of All selected True Neutrons:-" << nSumingAllTrueNeutron << std::endl;

  this->WriteHist();

  //Stores all particle information along with its PDG code for calculating contamination
  std::ofstream Delayed("ParticleCountsDelayed.csv", std::ios::trunc);
  Delayed << "PDG,Count\n";
  for (const auto& entry : ParticleCountsDelayed) {
    Delayed << entry.first << "," << entry.second << "\n"; // Write particle PDG and count to CSV                                                                                                             
  }
  Delayed.close();
  
  std::ofstream Prompt("ParticleCountsPrompt.csv", std::ios::trunc);
  Prompt << "PDG,Count\n";
  for (const auto& promptentry : ParticleCountsPrompt) {
    Prompt << promptentry.first << "," << promptentry.second << "\n"; // Write particle PDG and count to CSV                                                                                                                                                                                                                                                                                                            
  }
  Prompt.close();

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

  h_nTotalTrueNeutronsPromptTVtxXZ = new TH2F("h_nTotalTrueNeutronsPromptTVtxXZ", "h_nTotalTrueNeutronsPromptTVtxXZ", 40, -5, 5, 40, -5, 5);
  h_nTotalTrueNeutronsDelayedTVtxXZ = new TH2F("h_nTotalTrueNeutronsDelayedTVtxXZ", "h_nTotalTrueNeutronsDelayedTVtxXZ", 40, -5, 5, 40, -5, 5);
  h_nAllSelectedClustersPromptTVtxXZ = new TH2F("h_nAllSelectedClustersPromptTVtxXZ", "h_nAllSelectedClustersPromptTVtxXZ", 40, -5, 5, 40, -5, 5);
  h_nSelectedTrueNeutronsPromptTVtxXZ = new TH2F("h_nSelectedTrueNeutronsPromptTVtxXZ", "h_nSelectedTrueNeutronsPromptTVtxXZ", 40, -5, 5, 40, -5, 5);
  h_nAllSelectedClustersDelayedTVtxXZ = new TH2F("h_nAllSelectedClustersDelayedTVtxXZ", "h_nAllSelectedClustersDelayedTVtxXZ", 40, -5, 5, 40, -5, 5);
  h_nSelectedTrueNeutronsDelayedTVtxXZ = new TH2F("h_nSelectedTrueNeutronsDelayedTVtxXZ", "h_nSelectedTrueNeutronsDelayedTVtxXZ", 40, -5, 5, 40, -5, 5);

  h_nTotalTrueNeutronsPromptTVtxXY = new TH2F("h_nTotalTrueNeutronsPromptTVtxXY", "h_nTotalTrueNeutronsPromptTVtxXY", 40, -5, 5, 40, -5, 5);
  h_nTotalTrueNeutronsDelayedTVtxXY = new TH2F("h_nTotalTrueNeutronsDelayedTVtxXY", "h_nTotalTrueNeutronsDelayedTVtxXY", 40, -5, 5, 40, -5, 5);
  h_nAllSelectedClustersPromptTVtxXY = new TH2F("h_nAllSelectedClustersPromptTVtxXY", "h_nAllSelectedClustersPromptTVtxXY", 40, -5, 5, 40, -5, 5);
  h_nSelectedTrueNeutronsPromptTVtxXY = new TH2F("h_nSelectedTrueNeutronsPromptTVtxXY", "h_nSelectedTrueNeutronsPromptTVtxXY", 40, -5, 5, 40, -5, 5);
  h_nAllSelectedClustersDelayedTVtxXY = new TH2F("h_nAllSelectedClustersDelayedTVtxXY", "h_nAllSelectedClustersDelayedTVtxXY", 40, -5, 5, 40, -5, 5);
  h_nSelectedTrueNeutronsDelayedTVtxXY = new TH2F("h_nSelectedTrueNeutronsDelayedTVtxXY", "h_nSelectedTrueNeutronsDelayedTVtxXY", 40,-5, 5, 40, -5, 5);

  h_nTotalTrueNeutronsPromptTVtxYZ = new TH2F("h_nTotalTrueNeutronsPromptTVtxYZ", "h_nTotalTrueNeutronsPromptTVtxYZ", 40, -5, 5, 40, -5, 5);
  h_nTotalTrueNeutronsDelayedTVtxYZ = new TH2F("h_nTotalTrueNeutronsDelayedTVtxYZ", "h_nTotalTrueNeutronsDelayedTVtxYZ", 40, -5, 5, 40, -5, 5);
  h_nAllSelectedClustersPromptTVtxYZ = new TH2F("h_nAllSelectedClustersPromptTVtxYZ", "h_nAllSelectedClustersPromptTVtxYZ", 40, -5, 5, 40, -5, 5);
  h_nSelectedTrueNeutronsPromptTVtxYZ = new TH2F("h_nSelectedTrueNeutronsPromptTVtxYZ", "h_nSelectedTrueNeutronsPromptTVtxYZ", 40, -5, 5, 40, -5, 5);
  h_nAllSelectedClustersDelayedTVtxYZ = new TH2F("h_nAllSelectedClustersDelayedTVtxYZ", "h_nAllSelectedClustersDelayedTVtxYZ", 40, -5, 5, 40, -5, 5);
  h_nSelectedTrueNeutronsDelayedTVtxYZ = new TH2F("h_nSelectedTrueNeutronsDelayedTVtxYZ", "h_nSelectedTrueNeutronsDelayedTVtxYZ", 40, -5, 5, 40, -5, 5);
  
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

  h_allContaminationPrompt = new TH1F("h_allContaminationPrompt", "h_allContaminationPrompt", 150, 0, 150);
  h_allContaminationDelayed = new TH1F("h_allContaminationDelayed", "h_allContaminationDelayed", 150, 0, 150);

  //Combination of delayed and prompt
  h_nTotalTrueNeutronsTVtxXZ = new TH2F("h_nTotalTrueNeutronsAllTVtxXZ", "h_nTotalTrueNeutronsAllTVtxXZ", 40, -5, 5, 40, -5, 5);
  h_nSelectedTrueNeutronsTVtxXZ = new TH2F("h_nSelectedTrueNeutronsAllTVtxXZ", "h_nSelectedTrueNeutronsAllTVtxXZ", 40, -5, 5, 40, -5, 5);
  
  //Need to write histograms for the ClusterCharge and Charge Balance here and write it below the histograms
  h_nAllSelectedClusterPromptClusterCharge = new TH1F("h_nAllSelectedClusterPromptClusterCharge", "h_nAllSelectedClusterPromptClusterCharge", 50, 0, 300);
  h_nAllSelectedClusterPromptChargeBalance = new TH1F("h_nAllSelectedClusterPromptChargeBalance", "h_nAllSelectedClusterPromptChargeBalance", 50, 0, 1);
  h_nSelectedTrueNeutronsPromptClusterCharge = new TH1F("h_nSelectedTrueNeutronsPromptClusterCharge", "h_nSelectedTrueNeutronsPromptClusterCharge", 50, 0, 300);
  h_nSelectedTrueNeutronsPromptBalance = new TH1F("h_nSelectedTrueNeutronsPromptBalance", "h_nSelectedTrueNeutronsPromptBalance", 50, 0, 1);

  h_nAllSelectedClusterDelayedClusterCharge = new TH1F("h_nAllSelectedClusterDelayedClusterCharge", "h_nAllSelectedClusterDelayedClusterCharge", 50, 0, 300);
  h_nAllSelectedClusterDelayedChargeBalance = new TH1F("h_nAllSelectedClusterDelayedChargeBalance", "h_nAllSelectedClusterDelayedChargeBalance", 50, 0, 1);
  h_nSelectedTrueNeutronsDelayedClusterCharge = new TH1F("h_nSelectedTrueNeutronsDelayedClusterCharge", "h_nSelectedTrueNeutronsDelayedClusterCharge", 50, 0, 300);
  h_nSelectedTrueNeutronsDelayedBalance = new TH1F("h_nSelectedTrueNeutronsDelayedBalance", "h_nSelectedTrueNeutronsDelayedBalance", 50, 0, 1);

  h_nAllSelectedClusterPromptCBCC = new TH2F("h_nAllSelectedClusterPromptCBCC", "h_nAllSelectedClusterPromptCBCC", 50, 0, 300, 50, 0, 1);
  h_nSelectedTrueNeutronsPromptCBCC = new TH2F("h_nSelectedTrueNeutronsPromptCBCC", "h_nSelectedTrueNeutronsPromptCBCC", 50, 0, 300, 50, 0, 1);
  h_nAllSelectedClusterDelayedCBCC = new TH2F("h_nAllSelectedClusterDelayedCBCC", "h_nAllSelectedClusterDelayedCBCC", 50, 0, 300, 50, 0, 1);
  h_nSelectedTrueNeutronsDelayedCBCC = new TH2F("h_nSelectedTrueNeutronsDelayedCBCC", "h_nSelectedTrueNeutronsDelayedCBCC", 50, 0, 300, 50, 0, 1);
  

  gROOT->cd();

}

void SelectionEffnPurity::WriteHist()
{

  fOutFile->cd();
  TDirectory *dir_allhist = fOutFile->mkdir("Histograms");
  dir_allhist->cd();

  //combination of delayed and prompt hists
  h_nSelectedTrueNeutronsTVtxXZ->Write();
  h_nTotalTrueNeutronsTVtxXZ->Write();
  
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

  h_allContaminationPrompt->Write();
  h_allContaminationDelayed->Write();

  h_nAllSelectedClusterPromptClusterCharge->Write();
  h_nAllSelectedClusterPromptChargeBalance->Write();
  h_nSelectedTrueNeutronsPromptClusterCharge->Write();
  h_nSelectedTrueNeutronsPromptBalance->Write();

  h_nAllSelectedClusterDelayedClusterCharge->Write();
  h_nAllSelectedClusterDelayedChargeBalance->Write();
  h_nSelectedTrueNeutronsDelayedClusterCharge->Write();
  h_nSelectedTrueNeutronsDelayedBalance->Write();

  h_nAllSelectedClusterPromptCBCC->Write();
  h_nSelectedTrueNeutronsPromptCBCC->Write();

  h_nAllSelectedClusterDelayedCBCC->Write();
  h_nSelectedTrueNeutronsDelayedCBCC->Write();
  
  gROOT->cd();

}

bool SelectionEffnPurity::LoadFromStores()
{
  bool goodAnnieEvent = m_data->Stores.count("ANNIEEvent");
  if (!goodAnnieEvent) {
    logmessage = "1.1 SelectionEffnPurity:no ANNIEEvent store!";
    Log(logmessage, v_error, verbosity);
    return false;
  }
  
  bool goodClusterMap = m_data->Stores.at("ANNIEEvent")->Get(fClusterMapName, fClusterMap);
  if (!goodClusterMap) {
    logmessage = "1.2 SelectionEffnPurity: no " + fClusterMapName + " in the ANNIEEvent!";
    Log(logmessage, v_error, verbosity);
    return false;
  }
    
  bool goodMCParticles = m_data->Stores.at("ANNIEEvent")->Get("MCParticles", fMCParticles);
  if (!goodMCParticles) {
    logmessage = "1.3 SelectionEffnPurity:no MCParticles in the ANNIEEvent!";
    Log(logmessage, v_error, verbosity);
    return false;
  }
  
  bool goodMCParticleIndexMap = m_data->Stores.at("ANNIEEvent")->Get("TrackId_to_MCParticleIndex", fMCParticleIndexMap);
  if (!goodMCParticleIndexMap) {
    logmessage = "1.4 SelectionEffnPurity:no TrackId_to_MCParticleIndex in the ANNIEEvent!";
    Log(logmessage, v_error, verbosity);
    return false;
  }
  
  bool goodClusterToBestParticleID = m_data->Stores.at("ANNIEEvent")->Get("ClusterToBestParticleID", fClusterToBestParticleID);
  if (!goodClusterToBestParticleID) {
    logmessage = "1.5 SelectionEffnPurity:no ClusterToBestParticleID in the ANNIEEvent!";
    Log(logmessage, v_error, verbosity);
    return false;
  }
  
  bool goodClusterToBestParticlePDG = m_data->Stores.at("ANNIEEvent")->Get("ClusterToBestParticlePDG", fClusterToBestParticlePDG);
  if (!goodClusterToBestParticlePDG) {
    logmessage = "1.6 SelectionEffnPurity: no ClusterToBestParticlePDG in the ANNIEEvent!";
    Log(logmessage, v_error, verbosity);
    return false;
  }
  
  bool goodClusterEfficiency = m_data->Stores.at("ANNIEEvent")->Get("ClusterEfficiency", fClusterEfficiency);
  if (!goodClusterEfficiency) {
    logmessage = "1.7 SelectionEffnPurity:no ClusterEfficiency in the ANNIEEvent!";
    Log(logmessage, v_error, verbosity);
    return false;
  }
  
  bool goodClusterPurity = m_data->Stores.at("ANNIEEvent")->Get("ClusterPurity", fClusterPurity);
  if (!goodClusterPurity) {
    logmessage = "1.8 SelectionEffnPurity:no ClusterPurity in the ANNIEEvent!";
    Log(logmessage, v_error, verbosity);
    return false;
  }
  
  bool goodClusterTotalCharge = m_data->Stores.at("ANNIEEvent")->Get("ClusterTotalCharge", fClusterTotalCharge);
  if (!goodClusterTotalCharge) {
    logmessage = "1.9 SelectionEffnPurity:no ClusterTotalCharge in the ANNIEEvent!";
    Log(logmessage, v_error, verbosity);
    return false;
  }
  
  bool goodClusterChargeBalance =  m_data->Stores["ANNIEEvent"]->Get("ClusterChargeBalances", ClusterChargeBalances);
  if (!goodClusterChargeBalance){
    Log("SelectionEffnPurity tool: One of the charge cluster classifiers is not available", v_debug, verbosity);
    return false;
  }

  return true;
}

