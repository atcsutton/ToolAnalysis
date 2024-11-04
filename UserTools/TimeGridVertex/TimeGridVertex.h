#ifndef TimeGridVertex_H
#define TimeGridVertex_H

#include <string>
#include <iostream>

#include "Tool.h"
#include "Position.h"
#include "Hit.h"

#include "TFile.h"
#include "TTree.h"


/**
 * \class TimeGridVertex
 *
 * This is a blank template for a Tool used by the script to generate a new custom tool. Please fill out the description and author information.
*
* $Author: B.Richards $
* $Date: 2019/05/28 10:44:00 $
* Contact: b.richards@qmul.ac.uk
*/
class TimeGridVertex: public Tool {


 public:

  TimeGridVertex(); ///< Simple constructor
  bool Initialise(std::string configfile,DataModel &data); ///< Initialise Function for setting up Tool resources. @param configfile The path and name of the dynamic configuration file to read in. @param data A reference to the transient data class used to pass information between Tools.
  bool Execute(); ///< Execute function used to perform Tool purpose.
  bool Finalise(); ///< Finalise function used to clean up resources.

  // Generate a rectangular prism with n{XYZ} nodes in each dimension with the specified spacing
  std::vector<Position> GenerateVetices(const Position &center, double spacing);

  // Perform the actual tRMS calculation for each test point. Returns the pair <vertex Idx, tRMS> of the best match.
  std::pair<int, double> FindBestVertex(const std::vector<Hit> &hits, const std::vector<Position> &vertices); 
  std::pair<int, double> FindBestVertex(const std::vector<MCHit> &hits, const std::vector<Position> &vertices);

  // Run the while loop to which actually decreases the grid spacing to hone in on the vertex
  void RunLoop();
  void RunLoopMC();

  // Filter hits based on time separation and causality
  std::vector<MCHit> FilterHitsMC(std::vector<MCHit> hits);
    
  void SetupDebugTree();

 private:

  // Configuration parameters
  std::string fClusterMapName;
  bool fUseMCHits;
  double fInitialSpacing;
  double fMinSpacing;
  int fNumNodes;
  bool fDebugTree;
  
  // The ANNIE geometry service
  Geometry *fGeom = nullptr;

  // The clusters we'll load from the ANNIEEvent
  std::map<double, std::vector<Hit>>   *fClusterMap   = nullptr;
  std::map<double, std::vector<MCHit>> *fClusterMapMC = nullptr;

  // Vertices we'll save to the ANNIEEvent
  std::map<double, Position> *fVertexMap = nullptr;
      
  // The number of nodes in each dimension. These are calculated on the fly from the initial spacing
  int fNX;
  int fNY;
  int fNZ;

  // Output ROOT file for debug ttree
  TFile *fOutFile;
  TTree *fVtxTree;
  TTree *fHitTree;
  int fEventNum, fLoopCount;
  double fVtxX, fVtxY, fVtxZ;
  double fRMSt;
  double fHitX, fHitY, fHitZ;
  double fHitT0; 
  double fEff, fPur, fTotalQ, fNeutronQ;


  // Speed of light in water, m/ns
  const double fSoL = 0.299792458 * 3/4; 

  /// \brief verbosity levels: if 'verbosity' < this level, the message type will be logged.
  int verbosity;
  int v_error=0;
  int v_warning=1;
  int v_message=2;
  int v_debug=3;
  std::string logmessage;
};


#endif
