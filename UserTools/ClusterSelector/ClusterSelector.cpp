#include "ClusterSelector.h"

ClusterSelector::ClusterSelector():Tool(){}


bool ClusterSelector::Initialise(std::string configfile, DataModel &data)
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
    Log("ClusterSelector: \"verbosity\" not set in the config, defaulting to 0", v_error, verbosity);
  }

  bool gotClusterMapName = m_variables.Get("ClusterMapName", fClusterMapName);
  if (!gotClusterMapName) {
    Log("ClusterSelector: \"ClusterMapName\" not set in the config! Aborting!", v_error, verbosity);
    return false;
  }

  bool gotUseMCHits = m_variables.Get("UseMCHits", fUseMCHits);
  if (!gotUseMCHits) {
    Log("ClusterSelector: \"UseMCHits\" not set in the config! Aborting!", v_error, verbosity);
    return false;
  } else {
    bool gotPDGCut = m_variables.Get("PDGCut", fPDGCut);
    if (!gotPDGCut) fPDGCut = -9999;

    bool gotFindTrueCapts = m_variables.Get("FindTrueCapts", fFindTrueCapts);
    if (!gotFindTrueCapts) fFindTrueCapts = 0;
  }

  bool gotMaxClusterCharge = m_variables.Get("MaxClusterCharge", fMaxClusterCharge);
  if (!gotMaxClusterCharge) {
    Log("ClusterSelector: \"MaxClusterCharge\" not set in the config! Aborting!", v_error, verbosity);
    return false;
  }

  bool gotMaxClusterCB = m_variables.Get("MaxClusterCB", fMaxClusterCB);
  if (!gotMaxClusterCB) {
    Log("ClusterSelector: \"MaxClusterCB\" not set in the config! Aborting!", v_error, verbosity);
    return false;
  }

  bool gotCBvQIntercept = m_variables.Get("CBvQIntercept", fCBvQIntercept);
  if (!gotCBvQIntercept) {
    Log("ClusterSelector: \"CBvQIntercept\" not set in the config! Aborting!", v_error, verbosity);
    return false;
  }

  bool gotCBvQSlopeInverse = m_variables.Get("CBvQSlopeInverse", fCBvQSlopeInverse);
  if (!gotCBvQSlopeInverse) {
    Log("ClusterSelector: \"CBvQSlopeInverse\" not set in the config! Aborting!", v_error, verbosity);
    return false;
  }

  // Set up the pointers we're going to save. No need to 
  // delete them at Finalize, the store will handle it
  fClusterMapOut   = new std::map<double, std::vector<Hit>>;
  fClusterMapOutMC = new std::map<double, std::vector<MCHit>>;

  return true;
}

//------------------------------------------------------------------------------
bool ClusterSelector::Execute()
{
  // Grab cluster classifier outputs from the ANNIE Event
  std::map<double, double> clusterCBMap;
  bool gotClusterChargeBalance = m_data->Stores.at("ANNIEEvent")->Get("ClusterChargeBalances", clusterCBMap);
  if(!gotClusterChargeBalance){
    logmessage = "ClusterSelector: No ClusterChargeBalances found! Aborting!";
    Log(logmessage, v_error, verbosity);
    return false;
  }

  if (fUseMCHits) {
    
    bool gotClusters = m_data->CStore.Get("ClusterMapMC", fClusterMapMC);
    if(!gotClusters){
      logmessage = "ClusterSelector: No ClusterMapMC found! Aborting!";
      Log(logmessage, v_error, verbosity);
      return false;
    }

    // If we're using the PDG then load what we need from backtracker
    if (fPDGCut != -9999) {
      bool goodClusterToBestParticlePDG = m_data->Stores.at("ANNIEEvent")->Get("ClusterToBestParticlePDG", fClusterToBestParticlePDG);
      if (!goodClusterToBestParticlePDG) {
	logmessage = "EvaluateVertex: no ClusterToBestParticlePDG in the ANNIEEvent!";
	Log(logmessage, v_error, verbosity);
	return false;
      }
    }

    // If we're just considering true captures then load what we need from backtracker
    if (fFindTrueCapts) {
      bool goodClusterToBestParticleID = m_data->Stores.at("ANNIEEvent")->Get("ClusterToBestParticleID", fClusterToBestParticleID);
      if (!goodClusterToBestParticleID) {
	logmessage = "EvaluateVertex: no ClusterToBestParticleID in the ANNIEEvent!";
	Log(logmessage, v_error, verbosity);
	return false;
      }

    }

    fClusterMapOutMC->clear();
    FindClustersMC(clusterCBMap);
    
    m_data->Stores.at("ANNIEEvent")->Set(fClusterMapName,  fClusterMapOutMC);
    
  } else {
    
    bool gotClusters = m_data->CStore.Get("ClusterMap", fClusterMap);
    if(!gotClusters){
      logmessage = "ClusterSelector: No ClusterMap found! Aborting!";
      Log(logmessage, v_error, verbosity);
      return false;
    }

    fClusterMapOut->clear();
    FindClusters(clusterCBMap);
    m_data->Stores.at("ANNIEEvent")->Set(fClusterMapName, fClusterMapOut);
    
  }

  return true;
}

//------------------------------------------------------------------------------
bool ClusterSelector::Finalise()
{

  return true;
}

//------------------------------------------------------------------------------
void ClusterSelector::FindClusters(const std::map<double, double> &clusterCBMap)
{
  // Loop over input clusters
  for (auto cluster : *fClusterMap) {
    // Grab the charge balance
    double chargeBalance = clusterCBMap.at(cluster.first);

    // Loop over all hits to sum up the charge
    double chargeTotal = 0;
    for (auto hit : cluster.second)
      chargeTotal += hit.GetCharge();

    bool passCut = true;

    if ( passCut && chargeTotal > fMaxClusterCharge )
      passCut = false;
    
    if ( passCut && chargeBalance > fMaxClusterCB )
      passCut = false;

    if ( passCut && chargeBalance > (fCBvQIntercept + chargeTotal/fCBvQSlopeInverse) )
      passCut = false;
    
    if ( passCut )
      fClusterMapOut->emplace(cluster.first, cluster.second);
  }
}

//------------------------------------------------------------------------------
void ClusterSelector::FindClustersMC(const std::map<double, double> &clusterCBMap)
{
  if (fFindTrueCapts) {
    FindTrueCaptures();
    return;
  }
  
  
  // Loop over input clusters
  for (auto cluster : *fClusterMapMC) {
    // Grab the charge balance
    double chargeBalance = clusterCBMap.at(cluster.first);

    // Loop over all hits to sum up the charge
    double chargeTotal = 0;
    for (auto hit : cluster.second)
      chargeTotal += hit.GetCharge();

    bool passCut = true;

    if ( passCut && chargeTotal > fMaxClusterCharge )
      passCut = false;
    
    if ( passCut && chargeBalance > fMaxClusterCB )
      passCut = false;

    if ( passCut && chargeBalance > (fCBvQIntercept + chargeTotal/fCBvQSlopeInverse) )
      passCut = false;

    if ( passCut && fPDGCut != -9999 && fClusterToBestParticlePDG->at(cluster.first) != fPDGCut )      
      passCut = false;

    if ( passCut ) {
      fClusterMapOutMC->emplace(cluster.first, cluster.second);

    }
  }
}

//------------------------------------------------------------------------------
void ClusterSelector::FindTrueCaptures()
{
  // Map of neutron capture info
  std::map<std::string,std::vector<double>> fMCNeutCap;
  bool goodMCNeutCap = m_data->Stores.at("ANNIEEvent")->Get("MCNeutCap", fMCNeutCap);
  if (!goodMCNeutCap) {
    logmessage = "EvaluateVertex: no MCNeutCap in the ANNIEEvent!";
    Log(logmessage, v_error, verbosity);
    return;
  }

  // Bail if no neutron captures in this event
  if (!fMCNeutCap.size()) 
    return;

  // Grab the parent neutron's track ID and look for it in backtracker
  std::vector<double> cptPrtIDs = fMCNeutCap["CaptParent"];

  for (auto clusterTimeToParticleID : *fClusterToBestParticleID) {
    double clustertime = clusterTimeToParticleID.first;
    int id = clusterTimeToParticleID.second;
    
    for (auto cptPrtID : cptPrtIDs) {
      if (cptPrtID == id) {
  	fClusterMapOutMC->emplace(clustertime, fClusterMapMC->at(clustertime));
      }
    }
  }
}
