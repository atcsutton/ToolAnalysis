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

  bool gotNumNodes = m_variables.Get("NumNodes", fNumNodes);
  if (!gotNumNodes) {
    fNumNodes = 5;
    Log("TimeGridVertex: \"fNumNodes\" not set in the config! Using 5.", v_warning, verbosity); // 
  }

  
  bool gotDebugTree = m_variables.Get("DebugTree", fDebugTree);
  if (!gotDebugTree) {
    fDebugTree = 0;
    Log("TimeGridVertex: \"fDebugTree\" not set in the config! Using 0 (false).", v_warning, verbosity); // 
  } else {
    SetupDebugTree();
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
  //  fNX = fNumNodes; fNY = fNumNodes; fNZ = fNumNodes; 

  // Set up the pointer we're going to save. No need to 
  // delete it at Finalize, the store will handle it
  fVertexMap = new std::map<double, Position>;

  // Initialize some things
  fEventNum = 0;
  
  return true;
}

//------------------------------------------------------------------------------
bool TimeGridVertex::Execute()
{
  fVertexMap->clear();
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
  
  m_data->Stores.at("ANNIEEvent")->Set("TimeGridVertexMap",  fVertexMap);

  ++fEventNum;
  
  return true;
}

//------------------------------------------------------------------------------
bool TimeGridVertex::Finalise()
{
  if (fDebugTree) {
    fOutFile->cd();
    fVtxTree->Write();
    fHitTree->Write();
    fOutFile->Close();
  }
  
  return true;
}

//------------------------------------------------------------------------------
std::vector<Position> TimeGridVertex::GenerateVetices(const Position &center, double spacing)
{
  std::vector<Position> vertices;

  // if (fLoopCount == 0) {
  //     fNX = int(fGeom->GetPMTEnclosedRadius()*2 / fInitialSpacing);
  //     fNY = int(fGeom->GetPMTEnclosedHalfheight()*2 / fInitialSpacing);
  //     fNZ = int(fGeom->GetPMTEnclosedRadius()*2 / fInitialSpacing);
  // } else {
  //   fNX *= 1./0.75;
  //   fNY *= 1./0.75;
  //   fNZ *= 1./0.75;
  // }

  logmessage = "TimeGridVertex: Number of nodes (X, Y, Z): (" + std::to_string(fNX);
  logmessage += ", " + std::to_string(fNY) + ", " + std::to_string(fNZ) + ")";
  logmessage += "  Spacing: " + std::to_string(spacing);
  Log(logmessage, v_message, verbosity);


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

	// Check if this position is contained in the tank
	Position temppos(x[idX], y[idY], z[idZ]);
	double radialpos = sqrt( pow(temppos.X(), 2.) +
				 pow(temppos.Z() - fGeom->GetTankCentre().Z(), 2.) );
	double vertpos = abs(temppos.Y() - fGeom->GetTankCentre().Y());
	bool pmtcontained = ( (radialpos < fGeom->GetPMTEnclosedRadius()) &&
			      (vertpos < fGeom->GetPMTEnclosedHalfheight()) );

	if (!pmtcontained) continue;

	// if (!fGeom->GetTankContained(temppos)) continue;

	vertices.push_back(temppos);
      }
    }
  }

  logmessage = "TimeGridVertex: Number of vertices: " + std::to_string(vertices.size());
  Log(logmessage, v_message, verbosity);

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

    if (fDebugTree) {
      fVtxX = vertex.X();
      fVtxY = vertex.Y();
      fVtxZ = vertex.Z();
      fRMSt = tRMS;
      
      fVtxTree->Fill();
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
      double ti = hit.GetTime();
    
      Position dist = det->GetDetectorPosition() - vertex;
      
      // Predicted emission time from the vertex with this hit
      double t0 = ti - dist.Mag()/fSoL;
	    
      tsum += t0;
      tsum_sq += t0 * t0;
    }
    double tmean = tsum / hits.size();
    double tRMS = tsum_sq / hits.size() - tmean * tmean;

    if (tRMS < minRMS) {
      minRMS = tRMS;
      minVtxIdx = vtxIdx;
    }

    if (fDebugTree) {
      fVtxX = vertex.X();
      fVtxY = vertex.Y();
      fVtxZ = vertex.Z();
      fRMSt = tRMS;
      
      fVtxTree->Fill();
    }

  }

  // Debugging
  std::cout << "Best RMS: " << minRMS
  	    << " for vtx (x, y, z), " << vertices[minVtxIdx].X() << ", " << vertices[minVtxIdx].Y() << ", " << vertices[minVtxIdx].Z() << ""
  	    << std::endl;
  for (auto hit : hits) {
    Detector *det = fGeom->ChannelToDetector(hit.GetTubeId());
    double ti = hit.GetTime();

    Position detpos = det->GetDetectorPosition();
    Position dist = det->GetDetectorPosition() - vertices[minVtxIdx];
    double t0 = ti - dist.Mag()/fSoL;
    std::cout << " pos (x, y, z), " << detpos.X() << ", " << detpos.Y() << ", " << detpos.Z() << "" << std::endl;
    std::cout << " times (ti, tof, t0), " << ti << ", " << dist.Mag()/fSoL << ", " << t0 << "" << std::endl;
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
    fLoopCount = 0;
    bool done = false;      
    while (!done) {
      std::vector<Position> vertices = GenerateVetices(bestVertex, spacing);
      std::pair<int, double> bestIdxRMS = FindBestVertex(clusterpair.second, vertices);
    
      // Grab the new center position, decrease the spacing by half, and check if we should quit
      bestVertex = vertices[bestIdxRMS.first];
      tRMS = bestIdxRMS.second;
      spacing = spacing * 0.75;

      // When we hit the minimum spacing we want to perform one final go
      if (spacing < fMinSpacing) {
	spacing = fMinSpacing;

	std::vector<Position> vertices = GenerateVetices(bestVertex, spacing);
	std::pair<int, double> bestIdxRMS = FindBestVertex(clusterpair.second, vertices);
      
	bestVertex = vertices[bestIdxRMS.first];
	tRMS = bestIdxRMS.second;
	done = true;
      }
    
      ++fLoopCount;
      if (fLoopCount > 100) break;
    }

    logmessage = ("TimeGridVertex: Found best vertex at (X, Y, Z): (" +
		  std::to_string(bestVertex.X()) + ", " +
		  std::to_string(bestVertex.Y()) + ", " +
		  std::to_string(bestVertex.Z()) +
		  ") with RMS: " + std::to_string(tRMS));
    Log(logmessage, v_debug, verbosity);

    
    fVertexMap->emplace(clusterpair.first, bestVertex);
  }
}

//------------------------------------------------------------------------------
void TimeGridVertex::RunLoopMC()
{
  for (auto clusterpair : *fClusterMapMC) {
    // Build the initial test vertices around the tank center
    Position bestVertex = fGeom->GetTankCentre();
    double spacing = fInitialSpacing;

    // Filter hits
    auto filtered = FilterHitsMC(clusterpair.second);

    std::cout << "Num hits unfiltered: " << clusterpair.second.size() << ", filtered: " << filtered.size() << std::endl;
    
    // We need at least 4 hits
    if (filtered.size() < 4) continue;
   
    double tRMS = 9999;
    fLoopCount = 0;
    bool done = false;      
    while (!done) {
      std::vector<Position> vertices = GenerateVetices(bestVertex, spacing);
      std::pair<int, double> bestIdxRMS = FindBestVertex(filtered, vertices);
    
      // Grab the new center position, decrease the spacing by half, and check if we should quit
      bestVertex = vertices[bestIdxRMS.first];
      tRMS = bestIdxRMS.second;

      if (spacing == fMinSpacing) break;
      spacing = spacing * 0.75;

      // When we hit the minimum spacing we want to perform one final go
      if (spacing < fMinSpacing) {
	spacing = fMinSpacing;

	std::vector<Position> vertices = GenerateVetices(bestVertex, spacing);
	std::pair<int, double> bestIdxRMS = FindBestVertex(filtered, vertices);
      
	bestVertex = vertices[bestIdxRMS.first];
	tRMS = bestIdxRMS.second;
	done = true;
      }
    
      ++fLoopCount;
      if (fLoopCount > 100) break;
    }

    if (fDebugTree) {      
      for (auto hit : filtered) {
    	Detector *det = fGeom->ChannelToDetector(hit.GetTubeId());    	
    	Position detpos = det->GetDetectorPosition();

    	fHitX = detpos.X();
    	fHitY = detpos.Y();
    	fHitZ = detpos.Z();
	fHitT0 = hit.GetTime();
    	fHitTree->Fill();
      }
    }

    logmessage = ("TimeGridVertex: Found best vertex at (X, Y, Z): (" +
		  std::to_string(bestVertex.X()) + ", " +
		  std::to_string(bestVertex.Y()) + ", " +
		  std::to_string(bestVertex.Z()) +
		  ") with RMS: " + std::to_string(tRMS));
    Log(logmessage, v_debug, verbosity);
        
    fVertexMap->emplace(clusterpair.first, bestVertex);
  }
}

//------------------------------------------------------------------------------
std::vector<MCHit> TimeGridVertex::FilterHitsMC(std::vector<MCHit> hits)
{
  std::vector<MCHit> filtered;
  
  // Sort the hits by time
  std::sort(hits.begin(), hits.end(),
	    [](const MCHit &a, const MCHit &b) {return a.GetTime() < b.GetTime();});

  // Keep the first four and find the mean time
  double mean = 0;
  for (uint idx = 0; idx < 4; ++idx) {
    mean += hits[idx].GetTime();
    filtered.push_back(hits[idx]);
  }
  mean = mean/4.;

  // Look through the rest of the hits,
  // check that they are within 10 ns (for reflections),
  // and check that they could come from the same point
  for (uint idx = 4; idx < hits.size(); ++idx) {
    if ((hits[idx].GetTime() - mean) > 10)
      continue;

    bool keep = true;
    for (uint jdx = 0; jdx < 4; ++jdx) {
      Detector *det_i = fGeom->ChannelToDetector(hits[idx].GetTubeId());
      Position pos_i = det_i->GetDetectorPosition();

      Detector *det_j = fGeom->ChannelToDetector(filtered[jdx].GetTubeId());
      Position pos_j = det_j->GetDetectorPosition();

      double deltaT = abs(filtered[jdx].GetTime() - hits[idx].GetTime());
      double tof = (pos_i - pos_j).Mag() / fSoL;

      if (deltaT > tof) {
	keep = false;
	break;
      }
    }

    if (!keep) continue;
    filtered.push_back(hits[idx]);
  }

  return filtered;
}

//------------------------------------------------------------------------------
void TimeGridVertex::SetupDebugTree()
{
  fOutFile = new TFile("TimeGridVertex.debug.root", "RECREATE");

  fVtxTree = new TTree("vtxtree", "vtxtree");
  fVtxTree->Branch("EventNumber", &fEventNum );
  fVtxTree->Branch("CountNumber", &fLoopCount);
  fVtxTree->Branch("VtxX",        &fVtxX     );
  fVtxTree->Branch("VtxY",        &fVtxY     );
  fVtxTree->Branch("VtxZ",        &fVtxZ     );
  fVtxTree->Branch("TimeRMS",     &fRMSt     );

  fHitTree = new TTree("hittree", "hittree");
  fHitTree->Branch("EventNumber", &fEventNum );
  fHitTree->Branch("CountNumber", &fLoopCount);
  fHitTree->Branch("HitX",        &fHitX     );
  fHitTree->Branch("HitY",        &fHitY     );
  fHitTree->Branch("HitZ",        &fHitZ     );
  fHitTree->Branch("HitT0",       &fHitT0    );

}
