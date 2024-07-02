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

    if ( passCut )
      fClusterMapOutMC->emplace(cluster.first, cluster.second);
  }
}
