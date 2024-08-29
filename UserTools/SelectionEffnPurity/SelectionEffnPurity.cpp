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

    // std::cout << "Vertex X:-" << fTrueVtxX << "; Vertex Y:-" << fTrueVtxY << "; Vertex Z:-" << fTrueVtxZ << std::endl;
    
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
	      h_nSelectedTrueNeutronsPromptTVtxXY->Fill(fTrueVtxX, fTrueVtxY);
	      h_nSelectedTrueNeutronsPromptTVtxYZ->Fill(fTrueVtxY, fTrueVtxZ);
	      h_nSelectedTrueNeutronsPromptTVtxXZ->Fill(fTrueVtxX, fTrueVtxZ);
	      h_nSelectedTrueNeutronsPromptNE->Fill(NeutrinoEnergy);
	      nSelectedTrueNeutronsPrompt++;
	    }

	    else if (fBestPDG == 2212){nPromptProton++;
	    h_allContaminationPrompt->Fill(1, nPromptProton);} //proton                                                
	    //	    h_allContaminationPrompt->Fill(1, nPromptProton);
	      
	    else if (fBestPDG == -2212){nPromptAntiProton++;} //Anti-proton                                             
	    else if (fBestPDG == 11){nPromptElectron++;}
	    else if (fBestPDG == -11){nPromptPositron++;}
	    else if (fBestPDG == 12){nPromptElectronNeutrino++;}
	    else if (fBestPDG == -12){nPromptAntiElectronNeutrino++;}
	    else if (fBestPDG == 22) {nPromptGamma++;}
	    else if (fBestPDG == 2112) {nPromptNeutron++;}
	    else if (fBestPDG == -2112){nPromptAntiNeutron++;}
	    else if (fBestPDG == -13){nPromptMuonPlus++;}
	    else if (fBestPDG == 13){nPromptMuonMinus++;}
	    else if (fBestPDG == 130){nPromptKaonlong++;}
	    else if (fBestPDG == 211){nPromptPionPlus++;}
	    else if (fBestPDG == -211){nPromptPionMinus++;}
	    else if (fBestPDG == 321){nPromptKaonPlus++;}
	    else if (fBestPDG == -321){nPromptKaonMinus++;}
	    else if (fBestPDG == 310){nPromptKaonshort++;}
	    else if (fBestPDG == 111){nPromptPion0++;}
	    else if (fBestPDG == 311){nPromptKaon0++;}
	    else if (fBestPDG == 14){nPromptMuonNeutrino++;}
	    else if (fBestPDG == -14){nPromptAntiMuonNeutrino++;}
	    else if (fBestPDG == -15){nPromptTauPlus++;}
	    else if (fBestPDG == 15){nPromptTauMinus++;}

	    else if (fBestPDG== 3122){nPromptLambda++;}
	    else if (fBestPDG== -3122){nPromptAntiLambda++;}
	    else if (fBestPDG== 3112 ){nPromptSigmaMinus++;}
	    else if (fBestPDG== 3222){nPromptSigmaPlus++;}
	    else if (fBestPDG== 3212){nPromptSigma0++;}
	    else if (fBestPDG== -311){nPromptAntiKaon0++;}
	    else if (fBestPDG== -3222){nPromptAntiSigmaMinus++;}
	    else if (fBestPDG== -3212){nPromptAntiSigma0++;}
	    else if (fBestPDG== -3112){nPromptAntiSigmaPlus++;}
	    else if (fBestPDG== 3322){nPromptXsi0++;}
	    else if (fBestPDG== -3322){nPromptAntiXsi0++;}
	    else if (fBestPDG== 3312){nPromptXsiMinus++;}
	    else if (fBestPDG== -3312){nPromptXsiPlus++;}
	    else if (fBestPDG== 3334){nPromptOmegaMinus++;}
	    else if (fBestPDG== -3334){nPromptOmegaPlus++;}
	    else if (fBestPDG== 100){nPromptOpticalPhoton++;}
	    else if (fBestPDG== 3328){nPromptAlpha++;}
	    else if (fBestPDG== 3329){nPromptDeuteron++;}
	    else if (fBestPDG== 3330){nPromptTriton++;}
	    else if (fBestPDG== 3351){nPromptLi7++;}
	    else if (fBestPDG== 3331){nPromptC10++;}
	    else if (fBestPDG== 3345){nPromptB11++;}
	    else if (fBestPDG== 3332){nPromptC12++;}
	    else if (fBestPDG== 3350){nPromptC13++;}
	    else if (fBestPDG== 3349){nPromptN13++;}
	    else if (fBestPDG== 3340){nPromptN14++;}
	    else if (fBestPDG== 3333){nPromptN15++;}
	    else if (fBestPDG== 3334){nPromptN16++;}
	    else if (fBestPDG== 3335){nPromptO16++;}
	    else if (fBestPDG== 3346){nPromptAl27++;}
	    else if (fBestPDG== 3341){nPromptFe54++;}
	    else if (fBestPDG== 3348){nPromptMn54++;}
	    else if (fBestPDG== 3342){nPromptMn55++;}
	    else if (fBestPDG== 3352){nPromptMn56++;}
	    else if (fBestPDG== 3343){nPromptFe56++;}
	    else if (fBestPDG== 3344){nPromptFe57++;}
	    else if (fBestPDG== 3347){nPromptFe58++;}
	    else if (fBestPDG== 3353){nPromptEu154++;}
	    else if (fBestPDG== 3336){nPromptGd158++;}
	    else if (fBestPDG== 3337){nPromptGd156++;}
	    else if (fBestPDG== 3338){nPromptGd157++;}
	    else if (fBestPDG== 3339){nPromptGd155++;}

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
	    else if (fBestPDG == 2212){nDelayedProton++;} //proton                                                                                                                  
            else if (fBestPDG == -2212){nDelayedAntiProton++;} //Anti-proton                                                                                                        
            else if (fBestPDG == 11){nDelayedElectron++;}
            else if (fBestPDG == -11){nDelayedPositron++;}
            else if (fBestPDG == 12){nDelayedElectronNeutrino++;}
            else if (fBestPDG == -12){nDelayedAntiElectronNeutrino++;}
            else if (fBestPDG == 22) {nDelayedGamma++;}
            else if (fBestPDG == 2112) {nDelayedNeutron++;}
            else if (fBestPDG == -2112){nDelayedAntiNeutron++;}
            else if (fBestPDG == -13){nDelayedMuonPlus++;}
            else if (fBestPDG == 13){nDelayedMuonMinus++;}
            else if (fBestPDG == 130){nDelayedKaonlong++;}
            else if (fBestPDG == 211){nDelayedPionPlus++;}
            else if (fBestPDG == -211){nDelayedPionMinus++;}
            else if (fBestPDG == 321){nDelayedKaonPlus++;}
            else if (fBestPDG == -321){nDelayedKaonMinus++;}
            else if (fBestPDG == 310){nDelayedKaonshort++;}
            else if (fBestPDG == 111){nDelayedPion0++;}
            else if (fBestPDG == 311){nDelayedKaon0++;}
            else if (fBestPDG == 14){nDelayedMuonNeutrino++;}
            else if (fBestPDG == -14){nDelayedAntiMuonNeutrino++;}
            else if (fBestPDG == -15){nDelayedTauPlus++;}
            else if (fBestPDG == 15){nDelayedTauMinus++;}

	    else if (fBestPDG== 3122){nDelayedLambda++;}
            else if (fBestPDG== -3122){nDelayedAntiLambda++;}
            else if (fBestPDG== 3112 ){nDelayedSigmaMinus++;}
            else if (fBestPDG== 3222){nDelayedSigmaPlus++;}
            else if (fBestPDG== 3212){nDelayedSigma0++;}
            else if (fBestPDG== -311){nDelayedAntiKaon0++;}
            else if (fBestPDG== -3222){nDelayedAntiSigmaMinus++;}
            else if (fBestPDG== -3212){nDelayedAntiSigma0++;}
            else if (fBestPDG== -3112){nDelayedAntiSigmaPlus++;}
            else if (fBestPDG== 3322){nDelayedXsi0++;}
            else if (fBestPDG== -3322){nDelayedAntiXsi0++;}
            else if (fBestPDG== 3312){nDelayedXsiMinus++;}
            else if (fBestPDG== -3312){nDelayedXsiPlus++;}
            else if (fBestPDG== 3334){nDelayedOmegaMinus++;}
            else if (fBestPDG== -3334){nDelayedOmegaPlus++;}
            else if (fBestPDG== 100){nDelayedOpticalPhoton++;}
            else if (fBestPDG== 3328){nDelayedAlpha++;}
            else if (fBestPDG== 3329){nDelayedDeuteron++;}
            else if (fBestPDG== 3330){nDelayedTriton++;}
            else if (fBestPDG== 3351){nDelayedLi7++;}
            else if (fBestPDG== 3331){nDelayedC10++;}
            else if (fBestPDG== 3345){nDelayedB11++;}
            else if (fBestPDG== 3332){nDelayedC12++;}
            else if (fBestPDG== 3350){nDelayedC13++;}
            else if (fBestPDG== 3349){nDelayedN13++;}
            else if (fBestPDG== 3340){nDelayedN14++;}
            else if (fBestPDG== 3333){nDelayedN15++;}
            else if (fBestPDG== 3334){nDelayedN16++;}
            else if (fBestPDG== 3335){nDelayedO16++;}
            else if (fBestPDG== 3346){nDelayedAl27++;}
            else if (fBestPDG== 3341){nDelayedFe54++;}
            else if (fBestPDG== 3348){nDelayedMn54++;}
            else if (fBestPDG== 3342){nDelayedMn55++;}
            else if (fBestPDG== 3352){nDelayedMn56++;}
            else if (fBestPDG== 3343){nDelayedFe56++;}
            else if (fBestPDG== 3344){nDelayedFe57++;}
            else if (fBestPDG== 3347){nDelayedFe58++;}
            else if (fBestPDG== 3353){nDelayedEu154++;}
            else if (fBestPDG== 3336){nDelayedGd158++;}
            else if (fBestPDG== 3337){nDelayedGd156++;}
            else if (fBestPDG== 3338){nDelayedGd157++;}
            else if (fBestPDG== 3339){nDelayedGd155++;}


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

  //Containimation outputs
  std::cout << "Contamination of nPromptProton:-" << nPromptProton << std::endl;
  std::cout << "Contamination of nPromptAntiProton:-" << nPromptAntiProton << std::endl;
  std::cout << "Contamination of nPromptElectron:-" << nPromptElectron << std::endl;
  std::cout << "Contamination of nPromptPositron:-" << nPromptPositron <<std::endl;
  std::cout << "Contamination of nPromptElectronNeutrino:-" << nPromptElectronNeutrino<<std::endl;
  std::cout << "Contamination of nPromptAntiElectronNeutrino:-" << nPromptAntiElectronNeutrino <<std::endl;
  std::cout << "Contamination of nPromptGamma:-" << nPromptGamma<<std::endl;
  std::cout << "Contamination of nPromptNeutron:-" << nPromptNeutron<<std::endl;
  std::cout << "Contamination of nPromptAntiNeutron:-" << nPromptAntiNeutron<<std::endl;
  std::cout << "Contamination of nPromptMuonPlus:-" << nPromptMuonPlus<<std::endl;
  std::cout << "Contamination of nPromptMuonMinus:-" << nPromptMuonMinus<<std::endl;
  std::cout << "Contamination of nPromptKaonlong:-" << nPromptKaonlong<<std::endl;
  std::cout << "Contamination of nPromptPionPlus:-" << nPromptPionPlus<<std::endl;
  std::cout << "Contamination of nPromptPionMinus:-" << nPromptPionMinus<<std::endl;
  std::cout << "Contamination of nPromptKaonPlus:-" << nPromptKaonPlus<<std::endl;
  std::cout << "Contamination of nPromptKaonMinus:-" << nPromptKaonMinus<<std::endl;
  std::cout << "Contamination of nPromptKaonshort:-" << nPromptKaonshort<<std::endl;
  std::cout << "Contamination of nPromptPion0:-" << nPromptPion0<<std::endl;
  std::cout << "Contamination of nPromptKaon0:-" << nPromptKaon0<<std::endl;
  std::cout << "Contamination of nPromptMuonNeutrino:-" << nPromptMuonNeutrino<<std::endl;
  std::cout << "Contamination of nPromptAntiMuonNeutrino:-" << nPromptAntiMuonNeutrino<<std::endl;
  std::cout << "Contamination of nPromptTauPlus:-" << nPromptTauPlus<<std::endl;
  std::cout << "Contamination of nPromptTauMinus:-" << nPromptTauMinus<<std::endl;

  std::cout << "Contamination of nDelayedProton:-" << nDelayedProton<<std::endl;
  std::cout << "Contamination of nDelayedAntiProton:-" << nDelayedAntiProton<<std::endl;
  std::cout << "Contamination of nDelayedElectron:-" << nDelayedElectron<<std::endl;
  std::cout << "Contamination of nDelayedPositron:-" << nDelayedPositron<<std::endl;
  std::cout << "Contamination of nDelayedElectronNeutrino:-" << nDelayedElectronNeutrino<<std::endl;
  std::cout << "Contamination of nDelayedAntiElectronNeutrino:-" << nDelayedAntiElectronNeutrino<<std::endl;
  std::cout << "Contamination of nDelayedGamma:-" << nDelayedGamma<<std::endl;
  std::cout << "Contamination of nDelayedNeutron:-" << nDelayedNeutron<<std::endl;
  std::cout << "Contamination of nDelayedAntiNeutron:-" << nDelayedAntiNeutron<<std::endl;
  std::cout << "Contamination of nDelayedMuonPlus:-" << nDelayedMuonPlus<<std::endl;
  std::cout << "Contamination of nDelayedMuonMinus:-" << nDelayedMuonMinus<<std::endl;
  std::cout << "Contamination of nDelayedKaonlong:-" << nDelayedKaonlong<<std::endl;
  std::cout << "Contamination of nDelayedPionPlus:-" << nDelayedPionPlus<<std::endl;
  std::cout << "Contamination of nDelayedPionMinus:-" << nDelayedPionMinus<<std::endl;
  std::cout << "Contamination of nDelayedKaonPlus:-" << nDelayedKaonPlus<<std::endl;
  std::cout << "Contamination of nDelayedKaonMinus:-" << nDelayedKaonMinus<<std::endl;
  std::cout << "Contamination of nDelayedKaonshort:-" << nDelayedKaonshort<<std::endl;
  std::cout << "Contamination of nDelayedPion0:-" << nDelayedPion0<<std::endl;
  std::cout << "Contamination of nDelayedKaon0:-" << nDelayedKaon0<<std::endl;
  std::cout << "Contamination of nDelayedMuonNeutrino:-" << nDelayedMuonNeutrino<<std::endl;
  std::cout << "Contamination of nDelayedAntiMuonNeutrino:-" << nDelayedAntiMuonNeutrino<<std::endl;
  std::cout << "Contamination of nDelayedTauPlus:-" << nDelayedTauPlus<<std::endl;
  std::cout << "Contamination of nDelayedTauMinus:-" << nDelayedTauMinus<<std::endl;

  
  fOutTree->Fill();
  fOutTreeContPrompt->Fill();
  fOutTreeContDelayed->Fill();
  fOutFile->cd();
  fOutTree->Write();
  fOutTreeContPrompt->Write();
  fOutTreeContDelayed->Write();
  this->WriteHist();
  fOutFile->Close();
  return true;
}

void SelectionEffnPurity::SetupTTree()
{
  
  fOutFile = new TFile("SelectionEffnPurity.root", "RECREATE");
  fOutTree = new TTree("tree", "tree");
  fOutTreeContPrompt = new TTree("ContaminationPrompt", "ContaminationPrompt");  
  fOutTreeContDelayed = new TTree("ContaminationDelayed", "ContaminationDelayed");
  
  fOutTree->Branch("nTotalTrueNeutronsWorld",            &nTotalTrueNeutronsWorld);
  fOutTree->Branch("nTotalTrueNeutrons",                 &nTotalTrueNeutrons);
  fOutTree->Branch("nTotalTrueNeutronsPrompt",           &nTotalTrueNeutronsPrompt);
  fOutTree->Branch("nTotalTrueNeutronsDelayed",          &nTotalTrueNeutronsDelayed);
  fOutTree->Branch("nAllSelectedClustersWorld",          &nAllSelectedClustersWorld);
  fOutTree->Branch("nAllSelectedClustersPrompt",         &nAllSelectedClustersPrompt);
  fOutTree->Branch("nSelectedTrueNeutronsPrompt",        &nSelectedTrueNeutronsPrompt);
  fOutTree->Branch("nAllSelectedClustersDelayed",        &nAllSelectedClustersDelayed);
  fOutTree->Branch("nSelectedTrueNeutronsDelayed",       &nSelectedTrueNeutronsDelayed);

  fOutTreeContPrompt->Branch("nPromptProton", &nPromptProton);
  fOutTreeContPrompt->Branch("nPromptAntiProton", &nPromptAntiProton);
  fOutTreeContPrompt->Branch("nPromptElectron", &nPromptElectron);
  fOutTreeContPrompt->Branch("nPromptPositron", &nPromptPositron );
  fOutTreeContPrompt->Branch("nPromptElectronNeutrino", &nPromptElectronNeutrino );
  fOutTreeContPrompt->Branch("nPromptAntiElectronNeutrino", &nPromptAntiElectronNeutrino );
  fOutTreeContPrompt->Branch("nPromptGamma", &nPromptGamma );
  fOutTreeContPrompt->Branch("nPromptNeutron", &nPromptNeutron );
  fOutTreeContPrompt->Branch("nPromptAntiNeutron", &nPromptAntiNeutron );
  fOutTreeContPrompt->Branch("nPromptMuonPlus", &nPromptMuonPlus );
  fOutTreeContPrompt->Branch("nPromptMuonMinus", &nPromptMuonMinus );
  fOutTreeContPrompt->Branch("nPromptKaonlong", &nPromptKaonlong );
  fOutTreeContPrompt->Branch("nPromptPionPlus", &nPromptPionPlus );
  fOutTreeContPrompt->Branch("nPromptPionMinus", &nPromptPionMinus );
  fOutTreeContPrompt->Branch("nPromptKaonPlus", &nPromptKaonPlus );
  fOutTreeContPrompt->Branch("nPromptKaonMinus", &nPromptKaonMinus );
  fOutTreeContPrompt->Branch("nPromptKaonshort", &nPromptKaonshort );
  fOutTreeContPrompt->Branch("nPromptPion0", &nPromptPion0 );
  fOutTreeContPrompt->Branch("nPromptKaon0", &nPromptKaon0 );
  fOutTreeContPrompt->Branch("nPromptMuonNeutrino", &nPromptMuonNeutrino );
  fOutTreeContPrompt->Branch("nPromptAntiMuonNeutrino", &nPromptAntiMuonNeutrino );
  fOutTreeContPrompt->Branch("nPromptTauPlus", &nPromptTauPlus );
  fOutTreeContPrompt->Branch("nPromptTauMinus", &nPromptTauMinus );


  fOutTreeContPrompt->Branch("nPromptLambda", &nPromptLambda );
  fOutTreeContPrompt->Branch("nPromptAntiLambda", &nPromptAntiLambda );
  fOutTreeContPrompt->Branch("nPromptSigmaMinus", &nPromptSigmaMinus );
  fOutTreeContPrompt->Branch("nPromptSigmaPlus", &nPromptSigmaPlus );
  fOutTreeContPrompt->Branch("nPromptSigma0", &nPromptSigma0 );
  fOutTreeContPrompt->Branch("nPromptAntiKaon0", &nPromptAntiKaon0 );
  fOutTreeContPrompt->Branch("nPromptAntiSigmaMinus", &nPromptAntiSigmaMinus );
  fOutTreeContPrompt->Branch("nPromptAntiSigma0", &nPromptAntiSigma0 );
  fOutTreeContPrompt->Branch("nPromptAntiSigmaPlus", &nPromptAntiSigmaPlus );
  fOutTreeContPrompt->Branch("nPromptXsi0", &nPromptXsi0 );
  fOutTreeContPrompt->Branch("nPromptAntiXsi0", &nPromptAntiXsi0 );
  fOutTreeContPrompt->Branch("nPromptXsiMinus", &nPromptXsiMinus );
  fOutTreeContPrompt->Branch("nPromptXsiPlus", &nPromptXsiPlus );
  fOutTreeContPrompt->Branch("nPromptOmegaMinus", &nPromptOmegaMinus );
  fOutTreeContPrompt->Branch("nPromptOmegaPlus", &nPromptOmegaPlus );
  fOutTreeContPrompt->Branch("nPromptOpticalPhoton", &nPromptOpticalPhoton );
  fOutTreeContPrompt->Branch("nPromptAlpha", &nPromptAlpha );
  fOutTreeContPrompt->Branch("nPromptDeuteron", &nPromptDeuteron );
  fOutTreeContPrompt->Branch("nPromptTriton", &nPromptTriton );
  fOutTreeContPrompt->Branch("nPromptLi7", &nPromptLi7 );
  fOutTreeContPrompt->Branch("nPromptC10", &nPromptC10 );
  fOutTreeContPrompt->Branch("nPromptB11", &nPromptB11 );
  fOutTreeContPrompt->Branch("nPromptC12", &nPromptC12 );
  fOutTreeContPrompt->Branch("nPromptC13", &nPromptC13 );
  fOutTreeContPrompt->Branch("nPromptN13", &nPromptN13 );
  fOutTreeContPrompt->Branch("nPromptN14", &nPromptN14 );
  fOutTreeContPrompt->Branch("nPromptN15", &nPromptN15 );
  fOutTreeContPrompt->Branch("nPromptN16", &nPromptN16 );
  fOutTreeContPrompt->Branch("nPromptO16", &nPromptO16 );
  fOutTreeContPrompt->Branch("nPromptAl27", &nPromptAl27 );
  fOutTreeContPrompt->Branch("nPromptFe54", &nPromptFe54 );
  fOutTreeContPrompt->Branch("nPromptMn54", &nPromptMn54 );
  fOutTreeContPrompt->Branch("nPromptMn55", &nPromptMn55 );
  fOutTreeContPrompt->Branch("nPromptMn56", &nPromptMn56 );
  fOutTreeContPrompt->Branch("nPromptFe56", &nPromptFe56 );
  fOutTreeContPrompt->Branch("nPromptFe57", &nPromptFe57 );
  fOutTreeContPrompt->Branch("nPromptFe58", &nPromptFe58 );
  fOutTreeContPrompt->Branch("nPromptEu154", &nPromptEu154 );
  fOutTreeContPrompt->Branch("nPromptGd158", &nPromptGd158 );
  fOutTreeContPrompt->Branch("nPromptGd156", &nPromptGd156 );
  fOutTreeContPrompt->Branch("nPromptGd157", &nPromptGd157 );
  fOutTreeContPrompt->Branch("nPromptGd155", &nPromptGd155 );

  
  fOutTreeContDelayed->Branch("nDelayedProton", &nDelayedProton );
  fOutTreeContDelayed->Branch("nDelayedAntiProton", &nDelayedAntiProton );
  fOutTreeContDelayed->Branch("nDelayedElectron", &nDelayedElectron );
  fOutTreeContDelayed->Branch("nDelayedPositron", &nDelayedPositron );
  fOutTreeContDelayed->Branch("nDelayedElectronNeutrino", &nDelayedElectronNeutrino );
  fOutTreeContDelayed->Branch("nDelayedAntiElectronNeutrino", &nDelayedAntiElectronNeutrino );
  fOutTreeContDelayed->Branch("nDelayedGamma", &nDelayedGamma );
  fOutTreeContDelayed->Branch("nDelayedNeutron", &nDelayedNeutron );
  fOutTreeContDelayed->Branch("nDelayedAntiNeutron", &nDelayedAntiNeutron );
  fOutTreeContDelayed->Branch("nDelayedMuonPlus", &nDelayedMuonPlus );
  fOutTreeContDelayed->Branch("nDelayedMuonMinus", &nDelayedMuonMinus );
  fOutTreeContDelayed->Branch("nDelayedKaonlong", &nDelayedKaonlong );
  fOutTreeContDelayed->Branch("nDelayedPionPlus", &nDelayedPionPlus );
  fOutTreeContDelayed->Branch("nDelayedPionMinus", &nDelayedPionMinus );
  fOutTreeContDelayed->Branch("nDelayedKaonPlus", &nDelayedKaonPlus );
  fOutTreeContDelayed->Branch("nDelayedKaonMinus", &nDelayedKaonMinus );
  fOutTreeContDelayed->Branch("nDelayedKaonshort", &nDelayedKaonshort );
  fOutTreeContDelayed->Branch("nDelayedPion0", &nDelayedPion0 );
  fOutTreeContDelayed->Branch("nDelayedKaon0", &nDelayedKaon0 );
  fOutTreeContDelayed->Branch("nDelayedMuonNeutrino", &nDelayedMuonNeutrino );
  fOutTreeContDelayed->Branch("nDelayedAntiMuonNeutrino", &nDelayedAntiMuonNeutrino );
  fOutTreeContDelayed->Branch("nDelayedTauPlus", &nDelayedTauPlus );
  fOutTreeContDelayed->Branch("nDelayedTauMinus", &nDelayedTauMinus );


  fOutTreeContDelayed->Branch("nDelayedLambda", &nDelayedLambda);
  fOutTreeContDelayed->Branch("nDelayedAntiLambda", &nDelayedAntiLambda);
  fOutTreeContDelayed->Branch("nDelayedSigmaMinus", &nDelayedSigmaMinus);
  fOutTreeContDelayed->Branch("nDelayedSigmaPlus", &nDelayedSigmaPlus);
  fOutTreeContDelayed->Branch("nDelayedSigma0", &nDelayedSigma0);
  fOutTreeContDelayed->Branch("nDelayedAntiKaon0", &nDelayedAntiKaon0);
  fOutTreeContDelayed->Branch("nDelayedAntiSigmaMinus", &nDelayedAntiSigmaMinus);
  fOutTreeContDelayed->Branch("nDelayedAntiSigma0", &nDelayedAntiSigma0);
  fOutTreeContDelayed->Branch("nDelayedAntiSigmaPlus", &nDelayedAntiSigmaPlus);
  fOutTreeContDelayed->Branch("nDelayedXsi0", &nDelayedXsi0);
  fOutTreeContDelayed->Branch("nDelayedAntiXsi0", &nDelayedAntiXsi0);
  fOutTreeContDelayed->Branch("nDelayedXsiMinus", &nDelayedXsiMinus);
  fOutTreeContDelayed->Branch("nDelayedXsiPlus", &nDelayedXsiPlus);
  fOutTreeContDelayed->Branch("nDelayedOmegaMinus", &nDelayedOmegaMinus);
  fOutTreeContDelayed->Branch("nDelayedOmegaPlus", &nDelayedOmegaPlus);
  fOutTreeContDelayed->Branch("nDelayedOpticalPhoton", &nDelayedOpticalPhoton);
  fOutTreeContDelayed->Branch("nDelayedAlpha", &nDelayedAlpha);
  fOutTreeContDelayed->Branch("nDelayedDeuteron", &nDelayedDeuteron);
  fOutTreeContDelayed->Branch("nDelayedTriton", &nDelayedTriton);
  fOutTreeContDelayed->Branch("nDelayedLi7", &nDelayedLi7);
  fOutTreeContDelayed->Branch("nDelayedC10", &nDelayedC10);
  fOutTreeContDelayed->Branch("nDelayedB11", &nDelayedB11);
  fOutTreeContDelayed->Branch("nDelayedC12", &nDelayedC12);
  fOutTreeContDelayed->Branch("nDelayedC13", &nDelayedC13);
  fOutTreeContDelayed->Branch("nDelayedN13", &nDelayedN13);
  fOutTreeContDelayed->Branch("nDelayedN14", &nDelayedN14);
  fOutTreeContDelayed->Branch("nDelayedN15", &nDelayedN15);
  fOutTreeContDelayed->Branch("nDelayedN16", &nDelayedN16);
  fOutTreeContDelayed->Branch("nDelayedO16", &nDelayedO16);
  fOutTreeContDelayed->Branch("nDelayedAl27", &nDelayedAl27);
  fOutTreeContDelayed->Branch("nDelayedFe54", &nDelayedFe54);
  fOutTreeContDelayed->Branch("nDelayedMn54", &nDelayedMn54);
  fOutTreeContDelayed->Branch("nDelayedMn55", &nDelayedMn55);
  fOutTreeContDelayed->Branch("nDelayedMn56", &nDelayedMn56);
  fOutTreeContDelayed->Branch("nDelayedFe56", &nDelayedFe56);
  fOutTreeContDelayed->Branch("nDelayedFe57", &nDelayedFe57);
  fOutTreeContDelayed->Branch("nDelayedFe58", &nDelayedFe58);
  fOutTreeContDelayed->Branch("nDelayedEu154", &nDelayedEu154);
  fOutTreeContDelayed->Branch("nDelayedGd158", &nDelayedGd158);
  fOutTreeContDelayed->Branch("nDelayedGd156", &nDelayedGd156);
  fOutTreeContDelayed->Branch("nDelayedGd157", &nDelayedGd157);
  fOutTreeContDelayed->Branch("nDelayedGd155", &nDelayedGd155);

  
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

  h_allContaminationPrompt->Write();
  h_allContaminationDelayed->Write();
  
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
