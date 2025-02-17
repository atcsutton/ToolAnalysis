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

  // Load my config parameters
  bool gotVerbosity = m_variables.Get("verbosity",verbosity);
  if (!gotVerbosity) {
    verbosity = 0;
    logmessage = "BackTracker::Initialize: \"verbosity\" not set in the config, defaulting to 0";
    Log(logmessage, v_error, verbosity);
  }

  bool gotMCWaveforms = m_variables.Get("MCWaveforms", fMCWaveforms);
  if (!gotMCWaveforms) {
    fMCWaveforms = false;
    logmessage = "BackTracker::Initialize: \"MCWaveforms\" not set in the config, defaulting to false";
    Log(logmessage, v_error, verbosity);
  }


  // Set up the pointers we're going to save. No need to 
  // delete them at Finalize, the store will handle it
  fClusterToBestParticleID  = new std::map<double, int>;
  fClusterToBestParticlePDG = new std::map<double, int>;
  fClusterEfficiency        = new std::map<double, double>;
  fClusterPurity            = new std::map<double, double>;
  fClusterTotalCharge       = new std::map<double, double>;
  
  return true;
}

//------------------------------------------------------------------------------
bool BackTracker::Execute()
{
  int load_status = LoadFromStores();
  if (load_status == 0) return false;
  if (load_status == 2) return true;
  
  
  fClusterToBestParticleID ->clear();
  fClusterToBestParticlePDG->clear();
  fClusterEfficiency       ->clear();
  fClusterPurity           ->clear();
  fClusterTotalCharge      ->clear();

  fParticleToTankTotalCharge.clear();
  
  SumParticleTankCharge();

  if (fMCWaveforms) { // using clusters of Hits made from simulated PMT pulses

    // Produce the map from channel ID to pulse time to MCHit index
    //std::cout << "BT: Calling MapPulsesToParentIdxs" << std::endl;
    bool gotPulseMap = MapPulsesToParentIdxs();
    if (!gotPulseMap) {
      logmessage = "BackTracker: No good pulse map.";
      Log(logmessage, v_error, verbosity);
      return false;
    }
    
    // Loop over the clusters
    for (std::pair<double, std::vector<Hit>>&& apair : *fClusterMap) {
      // Create a vector of MCHits associated with the vector of Hits
      std::vector<MCHit> mcHits;
      for (auto hit : apair.second) {
	int channel_key = hit.GetTubeId();
	double hitTime = hit.GetTime();

	// Catches if something goes wrong
	// First make sure that we have the PMT in the outer map
	// Then make sure we have the hit time in the inner map
	if (fMapChannelToPulseTimeToMCHitIdx.find(channel_key) == fMapChannelToPulseTimeToMCHitIdx.end()) {
	  std::cout << "BackTracker: No hit on this PMT: " << channel_key << std::endl;
	  return false;
	}
	if (fMapChannelToPulseTimeToMCHitIdx.at(channel_key).find(hitTime) == fMapChannelToPulseTimeToMCHitIdx.at(channel_key).end()) {
	  std::cout << "BackTracker: No hit on this PMT: " << channel_key << " at this time: " << hitTime << std::endl;
	  return false;
	}

	std::vector<int> mcHitIdxVec = fMapChannelToPulseTimeToMCHitIdx[channel_key][hitTime];
	for (auto mcHitIdx : mcHitIdxVec) 
	  mcHits.push_back((fMCHitsMap->at(channel_key)).at(mcHitIdx));
      }// end loop over cluster hits

      int prtId = -5;
      int prtPdg = -5;
      double eff = -5;
      double pur = -5;
      double totalCharge = 0;

      MatchMCParticle(mcHits, prtId, prtPdg, eff, pur, totalCharge);

      fClusterToBestParticleID ->emplace(apair.first, prtId);
      fClusterToBestParticlePDG->emplace(apair.first, prtPdg);
      fClusterEfficiency       ->emplace(apair.first, eff);
      fClusterPurity           ->emplace(apair.first, pur);
      fClusterTotalCharge      ->emplace(apair.first, totalCharge);

      m_data->Stores.at("ANNIEEvent")->Set("MapChannelToPulseTimeToMCHitIdx",  fMapChannelToPulseTimeToMCHitIdx );
    }// end loop over the cluster map
  } else { // using clusters of MCHits
    // Loop over the MC clusters and do the things
    for (std::pair<double, std::vector<MCHit>>&& apair : *fClusterMapMC) {
      int prtId = -5;
      int prtPdg = -5;
      double eff = -5;
      double pur = -5;
      double totalCharge = 0;

      MatchMCParticle(apair.second, prtId, prtPdg, eff, pur, totalCharge);

      fClusterToBestParticleID ->emplace(apair.first, prtId);
      fClusterToBestParticlePDG->emplace(apair.first, prtPdg);
      fClusterEfficiency       ->emplace(apair.first, eff);
      fClusterPurity           ->emplace(apair.first, pur);
      fClusterTotalCharge      ->emplace(apair.first, totalCharge);
    }// end loop over the MC cluster map
  }// end if/else fMCWaveforms

  m_data->Stores.at("ANNIEEvent")->Set("ClusterToBestParticleID",  fClusterToBestParticleID );
  m_data->Stores.at("ANNIEEvent")->Set("ClusterToBestParticlePDG", fClusterToBestParticlePDG);
  m_data->Stores.at("ANNIEEvent")->Set("ClusterEfficiency",        fClusterEfficiency       );
  m_data->Stores.at("ANNIEEvent")->Set("ClusterPurity",            fClusterPurity           );
  m_data->Stores.at("ANNIEEvent")->Set("ClusterTotalCharge",       fClusterTotalCharge      );

  return true;
}

//------------------------------------------------------------------------------
bool BackTracker::Finalise()
{

  return true;
}

//------------------------------------------------------------------------------
void BackTracker::SumParticleTankCharge()
{
  for (auto apair : *fMCHitsMap) {
    std::vector<MCHit> mcHits = apair.second;
    for (uint mcHitIdx = 0; mcHitIdx < mcHits.size(); ++mcHitIdx) {

      // technically a MCHit could have multiple parents, but they don't appear to in practice
      // skip any cases we come across
      std::vector<int> parentIdxs = *(mcHits[mcHitIdx].GetParents());
      if (parentIdxs.size() != 1) continue;
      
      int particleId = -5;
      for (auto bpair : *fMCParticleIndexMap) {
	if (bpair.second == parentIdxs[0]) particleId = bpair.first;
      }
      if (particleId == -5) continue;
	
      double depositedCharge = mcHits[mcHitIdx].GetCharge();      
      if (!fParticleToTankTotalCharge.count(particleId)) 
	fParticleToTankTotalCharge.emplace(particleId, depositedCharge);
      else 
	fParticleToTankTotalCharge.at(particleId) += depositedCharge;
    }    
  }
}

//------------------------------------------------------------------------------
void BackTracker::MatchMCParticle(std::vector<MCHit> const &mchits, int &prtId, int &prtPdg, double &eff, double &pur, double &totalCharge)
{
  // Loop over the hits and get all of their parents and the energy that each one contributed
  //  be sure to bunch up all neutronic contributions
  std::map<int, double> mapParticleToTotalClusterCharge;
  totalCharge = 0;

  for (auto mchit : mchits) {    
    std::vector<int> parentIdxs = *(mchit.GetParents());
    if (parentIdxs.size() != 1) {
      logmessage = "BackTracker::MatchMCParticle: this MCHit has ";
      logmessage += std::to_string(parentIdxs.size()) + " parents!";
      Log(logmessage, v_debug, verbosity);
      continue;
    }
    
    int particleId = -5;
    for (auto apair : *fMCParticleIndexMap) {
      if (apair.second == parentIdxs[0]) particleId = apair.first;
    }
    if (particleId == -5) continue;
    
    double depositedCharge = mchit.GetCharge();
    totalCharge += depositedCharge;
    
    if (mapParticleToTotalClusterCharge.count(particleId) == 0) 
      mapParticleToTotalClusterCharge.emplace(particleId, depositedCharge);
    else
      mapParticleToTotalClusterCharge[particleId] += depositedCharge;    
  }       

  // Loop over the particleIds to find the primary contributer to the cluster
  double maxCharge = 0;
  for (auto apair : mapParticleToTotalClusterCharge) {
    if (apair.second > maxCharge) {
      maxCharge = apair.second;
      prtId = apair.first;
    }
  }

  // Check that we have some charge, if not then something is wrong so pass back all -5
  if (totalCharge > 0) {
    eff = maxCharge/fParticleToTankTotalCharge.at(prtId);
    pur = maxCharge/totalCharge;
    prtPdg = (fMCParticles->at(fMCParticleIndexMap->at(prtId))).GetPdgCode();
  } else {
    prtId = -5;
    eff = -5;
    pur = -5;
    totalCharge = -5;
  }

  logmessage = "BackTracker::MatchMCParticle: best particleId is : ";
  logmessage += std::to_string(prtId) + " which has PDG: " + std::to_string(prtPdg);
  Log(logmessage, v_message, verbosity);

}

//------------------------------------------------------------------------------
bool BackTracker::MapPulsesToParentIdxs()
{
  // Clear out the map
  fMapChannelToPulseTimeToMCHitIdx.clear();

  // Grab the pulses
  std::map<unsigned long, std::vector< std::vector<ADCPulse>> > adcPulseMap;
  bool goodADCPulses = m_data->Stores.at("ANNIEEvent")->Get("RecoADCHits", adcPulseMap);
  if (!goodADCPulses) {
    logmessage = "BackTracker: no RecoADCHits in the ANNIEEvent!";
    Log(logmessage, v_error, verbosity);
    return false;
  }

  // Loop over the ADCPulses and find MCHits that fall within the start and stop times
  // also record the pulse time to match the MCHits to the reco Hits
  // Pulses are indexed by PMT id then stacked in a two-deep vector

  for (auto apair: adcPulseMap) {
    int channel_key = apair.first;
      
    std::map<double, std::vector<int>> mapHitTimeToParents;
    bool goodPulses = false;
    for (auto pulseVec : apair.second) {
      for (auto pulse : pulseVec) {
	goodPulses = true;
	double pulseTime = pulse.peak_time();

	// Record the hit index if it occurred within the pulse window
	// If there is only one hit the no need to check the times
	std::vector<MCHit> mcHits = fMCHitsMap->at(channel_key);
	if (mcHits.size() == 1)
	  mapHitTimeToParents[pulseTime].push_back(0);
	else {
	  for (uint mcHitIdx = 0; mcHitIdx < mcHits.size(); ++mcHitIdx) {
	    double hitTime = mcHits[mcHitIdx].GetTime();
	    // The hit finding has to contend with noise so allow for a 10 ns 
	    // slew in the pulse start time (I know it seems large)
	    if ( hitTime + 10 >= pulse.start_time() && hitTime <= pulse.stop_time())
	      mapHitTimeToParents[pulseTime].push_back(mcHitIdx);
	  }// end loop over MCHits
	}
      }// end loop over inner pulse vector
    }// end loop over outer pulse vector

    if (mapHitTimeToParents.size() == 0 && goodPulses) {
      logmessage = "BackTracker::MapPulsesToParentIdxs: No MCHits match with this pulse! PMT channel: ";
      logmessage += std::to_string(channel_key);
      Log(logmessage, v_error, verbosity);
      return false;
    }

    fMapChannelToPulseTimeToMCHitIdx.emplace(channel_key, std::move(mapHitTimeToParents));
  }// end loop over pulse map

  return true;
}

//------------------------------------------------------------------------------
int BackTracker::LoadFromStores()
{
  // Grab the stuff we need from the stores

  bool goodAnnieEvent = m_data->Stores.count("ANNIEEvent");
  if (!goodAnnieEvent) {
    logmessage = "BackTracker: no ANNIEEvent store!";
    Log(logmessage, v_error, verbosity);
    return 0;
  }

  bool skip = false;
  bool goodSkipStatus = m_data->Stores.at("ANNIEEvent")->Get("SkipExecute", skip);
  if (goodSkipStatus && skip) {
    logmessage = "BackTracker: An upstream tool told me to skip this event.";
    Log(logmessage, v_warning, verbosity);
    return 2;
  }

  bool goodMCHits = m_data->Stores.at("ANNIEEvent")->Get("MCHits", fMCHitsMap);
  if (!goodMCHits) {
    logmessage = "BackTracker: no MCHits in the ANNIEEvent!";
    Log(logmessage, v_error, verbosity);
    return 0;
  }

  if (fMCWaveforms) {
    bool goodClusters = m_data->CStore.Get("ClusterMap", fClusterMap);
    if (!goodClusters) {
      logmessage = "BackTracker: no ClusterMap in the CStore!";
      Log(logmessage, v_error, verbosity);
      return 0;
    }

    bool goodRecoHits = m_data->Stores.at("ANNIEEvent")->Get("Hits", fRecoHitsMap);
    if (!goodRecoHits) {
      logmessage = "BackTracker: no Hits in the ANNIEEvent!";
      Log(logmessage, v_error, verbosity);
      return 0;
    }
    
  } else {
    bool goodMCClusters = m_data->CStore.Get("ClusterMapMC", fClusterMapMC);
    if (!goodMCClusters) {
      logmessage = "BackTracker: no ClusterMapMC in the CStore!";
      Log(logmessage, v_error, verbosity);
      return 0;
    }
  }// end if/else fMCWaveforms
  
  bool goodMCParticles = m_data->Stores.at("ANNIEEvent")->Get("MCParticles", fMCParticles);
  if (!goodMCParticles) {
    logmessage = "BackTracker: no MCParticles in the ANNIEEvent!";
    Log(logmessage, v_error, verbosity);
    return 0;
  }

  bool goodMCParticleIndexMap = m_data->Stores.at("ANNIEEvent")->Get("TrackId_to_MCParticleIndex", fMCParticleIndexMap);
  if (!goodMCParticleIndexMap) {
    logmessage = "BackTracker: no TrackId_to_MCParticleIndex in the ANNIEEvent!";
    Log(logmessage, v_error, verbosity);
    return 0;
  }

  return 1;
}

