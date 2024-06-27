#ifndef BackTracker_H
#define BackTracker_H

#include <string>
#include <iostream>

#include "Tool.h"
#include "Hit.h"
#include "Particle.h"


/**
 * \class BackTracker
 *
 * A tool to link reco info to the paticle(s) that generated the light 
*
* $Author: A.Sutton $
* $Date: 2024/06/16 $
* Contact: atcsutton@gmail.com
*/
class BackTracker: public Tool {


 public:

  BackTracker(); ///< Simple constructor
  bool Initialise(std::string configfile,DataModel &data); ///< Initialise Function for setting up Tool resources. @param configfile The path and name of the dynamic configuration file to read in. @param data A reference to the transient data class used to pass information between Tools.
  bool Execute(); ///< Execute function used to perform Tool purpose.
  bool Finalise(); ///< Finalise function used to clean up resources.

  bool LoadFromStores(); ///< Does all the loading so I can move it away from the Execute function
  void MatchMCParticle(std::vector<MCHit> const &mchits, int &prtId, int &prtPdg, double &eff, double &pur, double &totalCharge, double &neutronCharge); ///< The meat and potatoes
  
 private:

  // Things we need to pull out of the store
  std::map<double, std::vector<MCHit>>*           ClusterMapMC = nullptr;         ///< Clusters that we will be linking MCParticles to
  std::map<int, std::map<unsigned long, double>>* ParticleToTankTube = nullptr;   ///< Map between the particle Id and the charge deposited on a given tube
  std::map<int, double>*                          ParticleToTankCharge = nullptr; ///< Map between the particle Id and the total charge that it deposited on all tubes
  std::vector<MCParticle>*                        MCParticles = nullptr;          ///< The true particles from the event
  std::map<int, int>*                             MCParticleIndexMap = nullptr;   ///< Map between the particle Id and it's position in MCParticles vector
};


#endif
