#include "TimeGridVertex.h"

TimeGridVertex::TimeGridVertex():Tool(){}

//------------------------------------------------------------------------------
bool TimeGridVertex::Initialise(std::string configfile, DataModel &data)
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
    Log("TimeGridVertex: \"verbosity\" not set in the config, defaulting to 0", v_error, verbosity);
  }

  bool gotClusterMapName = m_variables.Get("ClusterMapName", fClusterMapName);
  if (!gotClusterMapName) {
    Log("TimeGridVertex: \"ClusterMapName\" not set in the config! Aborting!", v_error, verbosity);
    return false;
  }

  bool gotUseMCHits = m_variables.Get("UseMCHits", fUseMCHits);
  if (!gotUseMCHits) {
    Log("TimeGridVertex: \"UseMCHits\" not set in the config! Aborting!", v_error, verbosity);
    return false;
  }

  bool gotInitialSpacing = m_variables.Get("InitialSpacing", fInitialSpacing);
  if (!gotInitialSpacing) {
    fInitialSpacing = 0.5;
    Log("TimeGridVertex: \"InitialSpacing\" not set in the config! Using 0.5 m.", v_warning, verbosity);
  }

  bool gotMinSpacing = m_variables.Get("MinSpacing", fMinSpacing);
  if (!gotMinSpacing) {
    fMinSpacing = 0.005;
    Log("TimeGridVertex: \"fMinSpacing\" not set in the config! Using 0.005 (0.5 cm).", v_warning, verbosity); // 
  }


  // Load the geometry service
  bool gotGeometry = m_data->Stores.at("ANNIEEvent")->Header->Get("AnnieGeometry",fGeom);
  if(!gotGeometry){
    Log("TimeGridVertex: Error retrieving Geometry from ANNIEEvent! Aborting!", v_error, verbosity); 
    return false; 
  }

  // Calculate the number of vertex nodes based on the initial spacing
  fNX = int(fGeom->GetPMTEnclosedRadius()*2 / fInitialSpacing);
  fNY = int(fGeom->GetPMTEnclosedHalfheight()*2 / fInitialSpacing);
  fNZ = int(fGeom->GetPMTEnclosedRadius()*2 / fInitialSpacing);

  std::cout << "TimeGridVertex: Number of nodes (X, Y, Z): (" << fNX << ", " << fNY << ", " << fNZ << ")" << std::endl;

  fClusterToVertexMap = new std::map<double, Position>;
  
  return true;
}

//------------------------------------------------------------------------------
bool TimeGridVertex::Execute()
{
  fClusterToVertexMap->clear();
    
  if (fUseMCHits) {
    bool gotClusters = m_data->Stores.at("ANNIEEvent")->Get(fClusterMapName, fClusterMapMC);
    if (!gotClusters) {
      logmessage = "TimeGridVertex: no \"" + fClusterMapName + "\" in the ANNIEEvent!";
      Log(logmessage, v_error, verbosity);
      return false;
    }

    RunLoopMC();
  } else {
    bool gotClusters = m_data->Stores.at("ANNIEEvent")->Get(fClusterMapName, fClusterMap);
    if (!gotClusters) {
      logmessage = "TimeGridVertex: no \"" + fClusterMapName + "\" in the ANNIEEvent!";
      Log(logmessage, v_error, verbosity);
      return false;
    }

    RunLoop();
  }
  
  m_data->Stores.at("ANNIEEvent")->Set("ClusterToVertex",  fClusterToVertexMap);
  
  return true;
}

//------------------------------------------------------------------------------
bool TimeGridVertex::Finalise()
{

  return true;
}

//------------------------------------------------------------------------------
std::vector<Position> TimeGridVertex::GenerateVetices(const Position &center, double spacing)
{
  std::vector<Position> vertices;

  // Side lengths in each dimension
  double lX = fNX * spacing;
  double lY = fNY * spacing;
  double lZ = fNZ * spacing;


  double x[fNX] = {0};
  for (int idX = 0; idX < fNX; ++idX) {
    if (idX == 0) x[idX] = center.X() - lX/2.;
    else          x[idX] = x[idX-1] + spacing;
  }

  double y[fNY] = {0};
  for (int idY = 0; idY < fNY; ++idY) {
    if (idY == 0) y[idY] = center.Y() - lY/2.;
    else          y[idY] = y[idY-1] + spacing;
  }

  double z[fNZ] = {0};
  for (int idZ = 0; idZ < fNZ; ++idZ) {
    if (idZ == 0) z[idZ] = center.Z() - lZ/2.;
    else          z[idZ] = z[idZ-1] + spacing;
  }

  for (int idX = 0; idX < fNX; ++idX) {
    for (int idY = 0; idY < fNY; ++idY) {
      for (int idZ = 0; idZ < fNZ; ++idZ) {
	
	vertices.push_back(Position(x[idX], y[idY], z[idZ]));
      }
    }
  }

  // Cut out any vertices that are not inside the tank
  for (auto it = vertices.begin(); it != vertices.end(); ) {
    if (!fGeom->GetTankContained(*it))
      it = vertices.erase(it);
    else
      ++it;
  }
    
  return vertices;
}

//------------------------------------------------------------------------------
std::pair<int, double> TimeGridVertex::FindBestVertex(const std::vector<Hit> &hits, const std::vector<Position> &vertices)
{
  double minRMS = 9999;
  int minVtxIdx = -1;
  
  for (uint vtxIdx = 0; vtxIdx < vertices.size(); ++vtxIdx) {
    Position vertex = vertices[vtxIdx];
        
    double tsum = 0;
    double tsum_sq = 0;
    for (auto hit : hits) {
      Detector *det = fGeom->ChannelToDetector(hit.GetTubeId());
      double t0 = hit.GetTime();
    
      Position dist = det->GetDetectorPosition() - vertex;
      double t = t0 - dist.Mag()/fSoL;

      tsum += t;
      tsum_sq += t * t;
    }
    double tmean = tsum / hits.size();
    double tRMS = tsum_sq / hits.size() - tmean * tmean;

    if (tRMS < minRMS) {
      minRMS = tRMS;
      minVtxIdx = vtxIdx;
    }
  }

  return std::make_pair(minVtxIdx, sqrt(minRMS));
}

//------------------------------------------------------------------------------
std::pair<int, double> TimeGridVertex::FindBestVertex(const std::vector<MCHit> &hits, const std::vector<Position> &vertices)
{
  double minRMS = 9999;
  int minVtxIdx = -1;
  
  for (uint vtxIdx = 0; vtxIdx < vertices.size(); ++vtxIdx) {
    Position vertex = vertices[vtxIdx];
        
    double tsum = 0;
    double tsum_sq = 0;
    for (auto hit : hits) {
      Detector *det = fGeom->ChannelToDetector(hit.GetTubeId());
      double t0 = hit.GetTime();
    
      Position dist = det->GetDetectorPosition() - vertex;
      double t = t0 - dist.Mag()/fSoL;

      tsum += t;
      tsum_sq += t * t;
    }
    double tmean = tsum / hits.size();
    double tRMS = tsum_sq / hits.size() - tmean * tmean;

    if (tRMS < minRMS) {
      minRMS = tRMS;
      minVtxIdx = vtxIdx;
    }
  }

  return std::make_pair(minVtxIdx, sqrt(minRMS));
}

//------------------------------------------------------------------------------
void TimeGridVertex::RunLoop()
{
  for (auto clusterpair : *fClusterMap) {
    // Build the initial test vertices around the tank center
    Position bestVertex = fGeom->GetTankCentre();
    double spacing = fInitialSpacing;
    
    double tRMS = 9999;
    int counter = 0;
    bool done = false;      
    while (!done) {
      std::vector<Position> vertices = GenerateVetices(bestVertex, spacing);
      std::pair<int, double> bestIdxRMS = FindBestVertex(clusterpair.second, vertices);
    
      // Grab the new center position, decrease the spacing by half, and check if we should quit
      bestVertex = vertices[bestIdxRMS.first];
      tRMS = bestIdxRMS.second;
      spacing = spacing / 2.;

      // When we hit the minimum spacing we want to perform one final go
      if (spacing < fMinSpacing) {
	spacing = fMinSpacing;

	std::vector<Position> vertices = GenerateVetices(bestVertex, spacing);
	std::pair<int, double> bestIdxRMS = FindBestVertex(clusterpair.second, vertices);
      
	bestVertex = vertices[bestIdxRMS.first];
	tRMS = bestIdxRMS.second;
	done = true;
      }
    
      ++counter;
      if (counter > 100) break;
    }

    logmessage = ("TimeGridVertex: Found best vertex at (X, Y, Z): (" +
		  std::to_string(bestVertex.X()) + ", " +
		  std::to_string(bestVertex.Y()) + ", " +
		  std::to_string(bestVertex.Z()) +
		  ") with RMS: " + std::to_string(tRMS));
    Log(logmessage, v_debug, verbosity);

    
    fClusterToVertexMap->emplace(clusterpair.first, bestVertex);
  }
}

//------------------------------------------------------------------------------
void TimeGridVertex::RunLoopMC()
{
  for (auto clusterpair : *fClusterMapMC) {
    // Build the initial test vertices around the tank center
    Position bestVertex = fGeom->GetTankCentre();
    double spacing = fInitialSpacing;
    
    double tRMS = 9999;
    int counter = 0;
    bool done = false;      
    while (!done) {
      std::vector<Position> vertices = GenerateVetices(bestVertex, spacing);
      std::pair<int, double> bestIdxRMS = FindBestVertex(clusterpair.second, vertices);
    
      // Grab the new center position, decrease the spacing by half, and check if we should quit
      bestVertex = vertices[bestIdxRMS.first];
      tRMS = bestIdxRMS.second;
      spacing = spacing / 2.;

      // When we hit the minimum spacing we want to perform one final go
      if (spacing < fMinSpacing) {
	spacing = fMinSpacing;

	std::vector<Position> vertices = GenerateVetices(bestVertex, spacing);
	std::pair<int, double> bestIdxRMS = FindBestVertex(clusterpair.second, vertices);
      
	bestVertex = vertices[bestIdxRMS.first];
	tRMS = bestIdxRMS.second;
	done = true;
      }
    
      ++counter;
      if (counter > 100) break;
    }

    logmessage = ("TimeGridVertex: Found best vertex at (X, Y, Z): (" +
		  std::to_string(bestVertex.X()) + ", " +
		  std::to_string(bestVertex.Y()) + ", " +
		  std::to_string(bestVertex.Z()) +
		  ") with RMS: " + std::to_string(tRMS));
    Log(logmessage, v_debug, verbosity);
        
    fClusterToVertexMap->emplace(clusterpair.first, bestVertex);
  }
}
