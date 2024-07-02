#ifndef ClusterSelector_H
#define ClusterSelector_H

#include <string>
#include <iostream>

#include "Tool.h"


/**
 * \class ClusterSelector
 *
 * This is a blank template for a Tool used by the script to generate a new custom tool. Please fill out the description and author information.
*
* $Author: B.Richards $
* $Date: 2019/05/28 10:44:00 $
* Contact: b.richards@qmul.ac.uk
*/
class ClusterSelector: public Tool {


 public:

  ClusterSelector(); ///< Simple constructor
  bool Initialise(std::string configfile,DataModel &data); ///< Initialise Function for setting up Tool resources. @param configfile The path and name of the dynamic configuration file to read in. @param data A reference to the transient data class used to pass information between Tools.
  bool Execute(); ///< Execute function used to perform Tool purpose.
  bool Finalise(); ///< Finalise function used to clean up resources.

  void FindClusters(const std::map<double, double> &clusterCBMap);
  void FindClustersMC(const std::map<double, double> &clusterCBMap);
		      

 private:

  // Configuration variables
  double fMaxClusterCharge;
  double fMaxClusterCB;
  double fCBvQIntercept;
  double fCBvQSlopeInverse;
  bool fUseMCHits;
  std::string fClusterMapName;

  // The clusters we'll load from the CStore
  std::map<double, std::vector<Hit>>   *fClusterMap   = nullptr;
  std::map<double, std::vector<MCHit>> *fClusterMapMC = nullptr;

  // The clusters we'll save to the ANNIEEvent
  std::map<double, std::vector<Hit>>   *fClusterMapOut   = nullptr;
  std::map<double, std::vector<MCHit>> *fClusterMapOutMC = nullptr;

  /// \brief verbosity levels: if 'verbosity' < this level, the message type will be logged.
  int verbosity;
  int v_error=0;
  int v_warning=1;
  int v_message=2;
  int v_debug=3;
  std::string logmessage;

};


#endif
