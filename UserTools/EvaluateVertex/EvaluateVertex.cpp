#include "EvaluateVertex.h"

EvaluateVertex::EvaluateVertex():Tool(){}


bool EvaluateVertex::Initialise(std::string configfile, DataModel &data)
{

  /////////////////// Useful header ///////////////////////
  if(configfile!="") m_variables.Initialise(configfile); // loading config file
  //m_variables.Print();

  m_data= &data; //assigning transient data pointer
  /////////////////////////////////////////////////////////////////

  // Load my config parameters
  bool gotVerbosity = m_variables.Get("verbosity",verbosity);
  if (!gotVerbosity) {
    verbosity = 0;
    Log("EvaluateVertex: \"verbosity\" not set in the config, defaulting to 0", v_error, verbosity);
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

  SetupTTree();
  
  return true;
}

//------------------------------------------------------------------------------
bool EvaluateVertex::Execute()
{
  if (!LoadFromStores())
    return false;

  // Loop over all the clusters
  for (auto clusterPair : *fClusterMap) {
    // Grab all the stuff for this cluster
    double clusterTime = clusterPair.first;
    Position vertex = fVertexMap->at(clusterTime);
    int bestPrtID = fClusterToBestParticleID->at(clusterTime);
    int bestPrtIdx = fMCParticleIndexMap->at(bestPrtID);

    fBestPDG = fClusterToBestParticlePDG->at(clusterTime);
    fEff = fClusterEfficiency->at(clusterTime);
    fPur = fClusterPurity->at(clusterTime);
    fTotalQ = fClusterTotalCharge->at(clusterTime);
    fNeutronQ = fClusterNeutronCharge->at(clusterTime);
    if (fNeutronQ > fTotalQ/2.) fMoreNeutronQ = 1;
    else fMoreNeutronQ = 0;

    fRecoVtxX = vertex.X();
    fRecoVtxY = vertex.Y();
    fRecoVtxZ = vertex.Z();

    // The vertex should be the start location of the particle
    MCParticle bestParticle = fMCParticles->at(bestPrtIdx);
    Position trueVertex = bestParticle.GetStartVertex();
    fTrueVtxX = trueVertex.X();
    fTrueVtxY = trueVertex.Y();
    fTrueVtxZ = trueVertex.Z();

    fDistX = fRecoVtxX - fTrueVtxX;
    fDistY = fRecoVtxY - fTrueVtxY;
    fDistZ = fRecoVtxZ - fTrueVtxZ;
    fDist  = sqrt(fDistX*fDistX + fDistY*fDistY + fDistZ*fDistZ);
    
    fOutTree->Fill();
  }
  
  return true;
}

//------------------------------------------------------------------------------
bool EvaluateVertex::Finalise()
{
  fOutFile->cd();
  fOutTree->Write();
  fOutFile->Close();
  return true;
}

//------------------------------------------------------------------------------
void EvaluateVertex::SetupTTree()
{
  fOutFile = new TFile("EvaluateVertex.root", "RECREATE");
  fOutTree = new TTree("tree", "tree");

  fOutTree->Branch("TrueVtxX",     &fTrueVtxX    );
  fOutTree->Branch("TrueVtxY",     &fTrueVtxY    );
  fOutTree->Branch("TrueVtxZ",     &fTrueVtxZ    );
  fOutTree->Branch("RecoVtxX",     &fRecoVtxX    );
  fOutTree->Branch("RecoVtxY",     &fRecoVtxY    );
  fOutTree->Branch("RecoVtxZ",     &fRecoVtxZ    );
  fOutTree->Branch("DistX",        &fDistX       );
  fOutTree->Branch("DistY",        &fDistY       );
  fOutTree->Branch("DistZ",        &fDistZ       );
  fOutTree->Branch("Dist",         &fDist        );
  fOutTree->Branch("BestPDG",      &fBestPDG     );
  fOutTree->Branch("MoreNeutronQ", &fMoreNeutronQ);
  fOutTree->Branch("Efficiency",   &fEff         );
  fOutTree->Branch("Purity",       &fPur         );
  fOutTree->Branch("TotalQ",       &fTotalQ      );
  fOutTree->Branch("NeutronQ",     &fNeutronQ    );


}

//------------------------------------------------------------------------------
bool EvaluateVertex::LoadFromStores()
{
  bool goodAnnieEvent = m_data->Stores.count("ANNIEEvent");
  if (!goodAnnieEvent) {
    logmessage = "EvaluateVertex: no ANNIEEvent store!";
    Log(logmessage, v_error, verbosity);
    return false;
  }

  bool goodClusterMap = m_data->Stores.at("ANNIEEvent")->Get(fClusterMapName, fClusterMap);
  if (!goodClusterMap) {
    logmessage = "EvaluateVertex: no " + fClusterMapName + " in the ANNIEEvent!";
    Log(logmessage, v_error, verbosity);
    return false;
  }

  bool goodVertexMap = m_data->Stores.at("ANNIEEvent")->Get(fVertexMapName, fVertexMap);
  if (!goodVertexMap) {
    logmessage = "EvaluateVertex: no " + fVertexMapName + " in the ANNIEEvent!";
    Log(logmessage, v_error, verbosity);
    return false;
  }

  bool goodMCParticles = m_data->Stores.at("ANNIEEvent")->Get("MCParticles", fMCParticles);
  if (!goodMCParticles) {
    logmessage = "EvaluateVertex: no MCParticles in the ANNIEEvent!";
    Log(logmessage, v_error, verbosity);
    return false;
  }

  bool goodMCParticleIndexMap = m_data->Stores.at("ANNIEEvent")->Get("TrackId_to_MCParticleIndex", fMCParticleIndexMap);
  if (!goodMCParticleIndexMap) {
    logmessage = "EvaluateVertex: no TrackId_to_MCParticleIndex in the ANNIEEvent!";
    Log(logmessage, v_error, verbosity);
    return false;
  }

  bool goodClusterToBestParticleID = m_data->Stores.at("ANNIEEvent")->Get("ClusterToBestParticleID", fClusterToBestParticleID);
  if (!goodClusterToBestParticleID) {
    logmessage = "EvaluateVertex: no ClusterToBestParticleID in the ANNIEEvent!";
    Log(logmessage, v_error, verbosity);
    return false;
  }

  bool goodClusterToBestParticlePDG = m_data->Stores.at("ANNIEEvent")->Get("ClusterToBestParticlePDG", fClusterToBestParticlePDG);
  if (!goodClusterToBestParticlePDG) {
    logmessage = "EvaluateVertex: no ClusterToBestParticlePDG in the ANNIEEvent!";
    Log(logmessage, v_error, verbosity);
    return false;
  }

  bool goodClusterEfficiency = m_data->Stores.at("ANNIEEvent")->Get("ClusterEfficiency", fClusterEfficiency);
  if (!goodClusterEfficiency) {
    logmessage = "EvaluateVertex: no ClusterEfficiency in the ANNIEEvent!";
    Log(logmessage, v_error, verbosity);
    return false;
  }

  bool goodClusterPurity = m_data->Stores.at("ANNIEEvent")->Get("ClusterPurity", fClusterPurity);
  if (!goodClusterPurity) {
    logmessage = "EvaluateVertex: no ClusterPurity in the ANNIEEvent!";
    Log(logmessage, v_error, verbosity);
    return false;
  }

  bool goodClusterTotalCharge = m_data->Stores.at("ANNIEEvent")->Get("ClusterTotalCharge", fClusterTotalCharge);
  if (!goodClusterTotalCharge) {
    logmessage = "EvaluateVertex: no ClusterTotalCharge in the ANNIEEvent!";
    Log(logmessage, v_error, verbosity);
    return false;
  }

  bool goodClusterNeutronCharge = m_data->Stores.at("ANNIEEvent")->Get("ClusterNeutronCharge", fClusterNeutronCharge);
  if (!goodClusterNeutronCharge) {
    logmessage = "EvaluateVertex: no ClusterNeutronCharge in the ANNIEEvent!";
    Log(logmessage, v_error, verbosity);
    return false;
  }

  return true;
}
