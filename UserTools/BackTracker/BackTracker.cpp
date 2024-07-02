#include "BackTracker.h"

BackTracker::BackTracker():Tool(){}

// To sort
struct sort_by_charge {
    bool operator()(const std::pair<int,int> &left, const std::pair<int,int> &right) {
        return left.second < right.second;
    }
};


bool BackTracker::Initialise(std::string configfile, DataModel &data){

  /////////////////// Useful header ///////////////////////
  if(configfile!="") m_variables.Initialise(configfile); // loading config file
  //m_variables.Print();

  m_data= &data; //assigning transient data pointer
  /////////////////////////////////////////////////////////////////


  ClusterToBestParticleID  = new std::map<double, int>;
  ClusterToBestParticlePDG = new std::map<double, int>;
  ClusterEfficiency        = new std::map<double, double>;
  ClusterPurity            = new std::map<double, double>;
  ClusterTotalCharge       = new std::map<double, double>;
  ClusterNeutronCharge     = new std::map<double, double>;

  
  return true;
}

//------------------------------------------------------------------------------
bool BackTracker::Execute()
{
  if (!LoadFromStores())
    return false;

  ClusterToBestParticleID ->clear();
  ClusterToBestParticlePDG->clear();
  ClusterEfficiency       ->clear();
  ClusterPurity           ->clear();
  ClusterTotalCharge      ->clear();
  ClusterNeutronCharge    ->clear();

  
  // Loop over the clusters and do the things
  for (std::pair<double, std::vector<MCHit>>&& apair : *ClusterMapMC) {
    int prtId = -5;
    int prtPdg = -5;
    double eff = -5;
    double pur = -5;
    double totalCharge = 0;
    double neutronCharge = 0;

    MatchMCParticle(apair.second, prtId, prtPdg, eff, pur, totalCharge, neutronCharge);

    ClusterToBestParticleID ->emplace(apair.first, prtId);
    ClusterToBestParticlePDG->emplace(apair.first, prtPdg);
    ClusterEfficiency       ->emplace(apair.first, eff);
    ClusterPurity           ->emplace(apair.first, pur);
    ClusterTotalCharge      ->emplace(apair.first, totalCharge);
    ClusterNeutronCharge    ->emplace(apair.first, neutronCharge);
  }

  m_data->Stores.at("ANNIEEvent")->Set("ClusterToBestParticleID",  ClusterToBestParticleID );
  m_data->Stores.at("ANNIEEvent")->Set("ClusterToBestParticlePDG", ClusterToBestParticlePDG);
  m_data->Stores.at("ANNIEEvent")->Set("ClusterEfficiency",        ClusterEfficiency       );
  m_data->Stores.at("ANNIEEvent")->Set("ClusterPurity",            ClusterPurity           );
  m_data->Stores.at("ANNIEEvent")->Set("ClusterTotalCharge",       ClusterTotalCharge      );
  m_data->Stores.at("ANNIEEvent")->Set("ClusterNeutronCharge",     ClusterNeutronCharge    );

  return true;
}

//------------------------------------------------------------------------------
bool BackTracker::Finalise()
{

  return true;
}

//------------------------------------------------------------------------------
bool BackTracker::LoadFromStores()
{
  // Grab the stuff we need from the stores
  bool goodMCClusters = m_data->CStore.Get("ClusterMapMC", ClusterMapMC);
  if (!goodMCClusters) {
    std::cerr<<"BackTracker: no ClusterMapMC in the CStore!"<<endl;
    return false;
  }

  bool goodAnnieEvent = m_data->Stores.count("ANNIEEvent");
  if (!goodAnnieEvent) {
    std::cerr<<"BackTracker: no ANNIEEvent store!"<<endl;
    return false;
  }
    
  bool goodParticleTankTubeMap = m_data->Stores.at("ANNIEEvent")->Get("ParticleId_to_TankTubeIds", ParticleToTankTube);
  if (!goodParticleTankTubeMap) {
    std::cerr<<"BackTracker: no ParticleId_to_TankTubeIds in the ANNIEEvent!"<<endl;
    return false;
  }

  bool goodParticleTankChargeMap = m_data->Stores.at("ANNIEEvent")->Get("ParticleId_to_TankCharge", ParticleToTankCharge);
  if (!goodParticleTankChargeMap) {
    std::cerr<<"BackTracker: no ParticleId_to_TankCharge in the ANNIEEvent!"<<endl;
    return false;
  }
  
  bool goodMCParticles = m_data->Stores.at("ANNIEEvent")->Get("MCParticles", MCParticles);
  if (!goodMCParticles) {
    std::cerr<<"BackTracker: no MCParticles in the ANNIEEvent!"<<endl;
    return false;
  }

  bool goodMCParticleIndexMap = m_data->Stores.at("ANNIEEvent")->Get("TrackId_to_MCParticleIndex", MCParticles);
  if (!goodMCParticleIndexMap) {
    std::cerr<<"BackTracker: no TrackId_to_MCParticleIndex in the ANNIEEvent!"<<endl;
    return false;
  }

  return true;
}

//------------------------------------------------------------------------------
void BackTracker::MatchMCParticle(std::vector<MCHit> const &mchits, int &prtId, int &prtPdg, double &eff, double &pur, double &totalCharge, double &neutronCharge)
{
  // Loop over the hits and get all of their parents and the energy that each one contributed
  //  be sure to bunch up all neutronic contributions
  std::map<int, double> map_ParticleId_TotalClusterCharge;
  totalCharge = 0;
  neutronCharge = 0;  
  for (auto mchit : mchits) {
    unsigned long tubeId = mchit.GetTubeId();
    for (auto particleId : *mchit.GetParents()) {
      double depositedCharge = ParticleToTankTube->at(particleId).at(tubeId);
      totalCharge += depositedCharge;

      if (map_ParticleId_TotalClusterCharge.count(particleId) == 0) 
	map_ParticleId_TotalClusterCharge.emplace(particleId, depositedCharge);
      else
	map_ParticleId_TotalClusterCharge[particleId] += depositedCharge;

      auto tempParticle = MCParticles->at(MCParticleIndexMap->at(particleId));
      if (tempParticle.GetParentPdg() == 2112) 
	neutronCharge += depositedCharge;
    }
  }

  // Loop over the particleIds to find the primary contributer to the cluster
  double maxCharge = 0;
  for (auto apair : map_ParticleId_TotalClusterCharge) {
    if (apair.second > maxCharge) {
      maxCharge = apair.second;
      prtId = apair.first;
    }
  }


  // Check that we have some charge, if not then something is wrong so pass back all -5
  if (totalCharge > 0) {
    eff = maxCharge/ParticleToTankCharge->at(prtId);
    pur = maxCharge/totalCharge;
    prtPdg = (MCParticles->at(MCParticleIndexMap->at(prtId))).GetPdgCode();
  } else {
    prtId = -5;
    eff = -5;
    pur = -5;
    totalCharge = -5;
    neutronCharge = -5;
  }
}
