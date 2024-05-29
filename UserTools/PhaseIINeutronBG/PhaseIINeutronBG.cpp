#include "PhaseIINeutronBG.h"
#include <iostream>
#include <fstream>

PhaseIINeutronBG::PhaseIINeutronBG():Tool(){}

bool PhaseIINeutronBG::Initialise(std::string configfile, DataModel &data){

  
  ///////////////////////// Useful header /////////////////////////
  if (configfile != "") m_variables.Initialise(configfile); // loading config file
  //m_variables.Print();

  m_data= &data; // assigning transient data pointer
  /////////////////////////////////////////////////////////////////

  // default variable values
  verbosity = 1;
  outpfile_prefix = "PhaseIINBG_default_output";
  min_clusterPE = 5;
  max_clusterPE = 100;
  BeamOn = 1;
  pethresh = 5.;
  PMTthresh = 1;  
  goodnoncc = 0;  
  goodnoncctrig = 0;

  // load user-defined variables from config file 
  m_variables.Get("verbosity", verbosity);
  m_variables.Get("OutputFilePrefix", outpfile_prefix);
  m_variables.Get("ClusterPEMin", min_clusterPE);
  m_variables.Get("ClusterPEMax", max_clusterPE);
  m_variables.Get("BeamOn", BeamOn);
  if (verbosity > 2)
  {
    std::cout << " PhaseIINeutronBG tool: config variables loaded. Moving on.." << std::endl;
  }

  std::cout << "BeamOn is set to " << BeamOn << std::endl;

  // must load geo in Initialise or else it gets overwritten
  bool get_ok = m_data->Stores["ANNIEEvent"]->Header->Get("AnnieGeometry", fGeo);
  if (!get_ok)
  {
    std::cout << " PhaseIINeutronBG tool: UH OH! Could not load AnnieGeometry!" << std::endl;
    return false;
  }

  // Things have changed slightly. ClusterClassifiers now puts some stuff into the RecoEvent store
  // So we have to make sure that exists
  if (!m_data->Stores.count("RecoEvent"))
    m_data->Stores["RecoEvent"] = new BoostStore(false, 2);


  // initialize ROOT stuff
  std::string root_outpfile_ext = ".root";
  std::string root_outpfile_name = outpfile_prefix + root_outpfile_ext;

  if (verbosity > 2) { std::cout << " PhaseIINeutronBG tool: root file will be saved as " << root_outpfile_name << std::endl; }

  p2nbg_root_outp = new TFile(root_outpfile_name.c_str(), "RECREATE");

  this->InitTree();
  this->InitHist();

  m_data->CStore.Get("ChannelNumToTankPMTSPEChargeMap", map_chankey2spe);

  std::cout << " PhaseIINeutronBG tool: Initalization complete." << std::endl;
  if (verbosity > 3)
  {
    std::cout << "   Moving on to world destru--I mean, neutron background analysis........" << std::endl;
  }

  return true;
}


bool PhaseIINeutronBG::Execute(){
  
  if (verbosity > 5) { std::cout << "MUCH VERBOSE! SO WORDS! MANY TALKS!" << std::endl; }

  //  load event info  //
  int annieeventexists = m_data->Stores.count("ANNIEEvent");
  if (!annieeventexists)
  {
    std::cerr << "PhaseIINeutronBG tool: No ANNIEEvent store!" << std::endl;
    return false;
  }

  nevents++;
  goodnoncctrig=0;
  bool get_ok;
  get_ok = m_data->Stores["ANNIEEvent"]->Get("RunNumber", fRunNumber);
  if (!get_ok) { Log("PhaseIINeutronBG tool: No RunNumber object in ANNIEEvent! Abort!", v_error, verbosity); return true; }
  get_ok = m_data->Stores["ANNIEEvent"]->Get("SubrunNumber", fSubrunNumber);
  if (!get_ok) { Log("PhaseIINeutronBG tool: No SubrunNumber object in ANNIEEvent! Abort!", v_error, verbosity); return true; }
  get_ok = m_data->Stores["ANNIEEvent"]->Get("RunType", fRunType);
  if (!get_ok) { Log("PhaseIINeutronBG tool: No RunType object in ANNIEEvent! Abort!", v_error, verbosity); return true; }
  get_ok = m_data->Stores["ANNIEEvent"]->Get("RunStartTime", fStartTime);
  if (!get_ok) { Log("PhaseIINeutronBG tool: No RunStartTime object in ANNIEEvent! Abort!", v_error, verbosity); return true; }
  get_ok = m_data->Stores["ANNIEEvent"]->Get("EventNumber", fEventNumber);
  if (!get_ok) { Log("PhaseIINeutronBG tool: No EventNumber object in ANNIEEvent! Abort!", v_error, verbosity); return true; }
  get_ok = m_data->Stores["ANNIEEvent"]->Get("EventTimeTank", fEventTimeTank);
  if (!get_ok) { Log("PhaseIINeutronBG tool: No EventTimeTank object in ANNIEEvent! Abort!", v_error, verbosity); return true; }
  get_ok = m_data->Stores["ANNIEEvent"]->Get("BeamStatus", beamstat);
  if (!get_ok) { Log("PhaseIINeutronBG tool: No BeamStatus object in ANNIEEvent! Abort!", v_error, verbosity); return true;}

  
  bool isBeam = false;
  bool isExtTrig = false;
  bool isLED = false;
  get_ok = m_data->Stores["ANNIEEvent"]->Get("TriggerWord", trigword);
  if (!get_ok) { Log("PhaseIINeutronBG tool: No TriggerWord object in ANNIEEvent! Abort!", v_error, verbosity); return true; }
  get_ok = m_data->Stores["ANNIEEvent"]->Get("TriggerExtended", trigext);
  if (!get_ok) { Log("PhaseIINeutronBG tool: No TriggerExtended object in ANNIEEvent! Abort!", v_error, verbosity); return true; }
  if (trigword == 5) { isBeam = true; } //beam trigger
  if (trigext == 2) { isExtTrig = true; } //noncc (forced extended window)
  if (trigword == 31) { isLED = true; } // LED triggerword, used as beam off
  bool isnoncc = isBeam && isExtTrig; 
  bool noprompt = true;
  fpot = beamstat.pot(); // protons on target (delivered by beam)
  totalpot += fpot;
  beamflag = -1;

  
  if (isLED == true) //beam off, readout window is 2us
    {
       beamoff_lt += 2;
       beamoff_pot += fpot;
    }

  if (isBeam && trigext == 0) //beam on without extended trigger, readout window is 2us
    {
      woext_lt += 2;
      woext_pot +=fpot;
    }

  if (isBeam && trigext == 1) //Beam on with CC extended trigger, readout window is 70us
    {
      cc_lt += 70;
      cc_pot += fpot;
    }

  if (isBeam && trigext == 2) //Beam on with NC extended trigger, readout window is 70us
    {
      nc_lt += 70;
      nc_pot += fpot;
    }



  if (beamstat.condition() == BeamCondition::NonBeamMinibuffer) {beamflag=1;} // beam is off
  else if (beamstat.condition() == BeamCondition::Ok) {beamflag=0;} // beam is on
  else if (beamstat.condition() == BeamCondition::Bad) {beamflag=2;} // beam is funky
  else {beamflag=3;}
  
  h_beamFlag->Fill(beamflag);

  //  check for Veto hits  //
  fVetoHit = 0;
  bool hasVeto = false;

  get_ok = m_data->Stores["ANNIEEvent"]->Get("TDCData", tdcdata);
  if (!get_ok) { Log("PhaseIINeutronBG tool: No TDCData object in ANNIEEvent! Abort!", v_error, verbosity); return true; }

  if (tdcdata->size() > 0)
  {
    Log("PhaseIINeutronBG tool: Looping over FMV/MRD hits... looking for Veto activity", v_debug, verbosity);
    for (auto&& anmrdpmt : (*tdcdata))
    {
      unsigned long chankey = anmrdpmt.first;
      Detector *thedetector = fGeo->ChannelToDetector(chankey);
      unsigned long detkey = thedetector->GetDetectorID();
      if (thedetector->GetDetectorElement() == "Veto") { fVetoHit = 1; hasVeto = true;}
    }
  }

  //  check for MRD tracks  //
  bool hasMRDTracks = false;
  int n_tracks_evt = -99;
  fnMRDTracks = -99;
  get_ok = m_data->Stores["MRDTracks"]->Get("NumMrdTracks", n_tracks_evt);
  if (!get_ok) { Log("PhaseIINeutronBG tool: No NumMrdTracks object in ANNIEEvent! Abort!", v_error, verbosity); return true; }
  fnMRDTracks = n_tracks_evt;
  
  if (n_tracks_evt != 0)
  {
    hasMRDTracks = true; 
  }


  //  looking at tank hits  //
  std::map<unsigned long, std::vector<Hit>> *tank_hits = nullptr;

  get_ok = m_data->Stores["ANNIEEvent"]->Get("Hits", tank_hits);
  if (!get_ok) { Log("PhaseIINeutronBG tool: No Hits object in ANNIEEvent! Abort!", v_error, verbosity); return true; }

  bool t_hitPE_ok = true;
  fHitPE.clear();
  fHitPMT.clear();
  for (std::pair<unsigned long, std::vector<Hit>>&& apair : *tank_hits)
  { 
    unsigned long chankey = apair.first;
    fHitPMT.push_back(chankey);
    std::map<int, double>::iterator it = map_chankey2spe.find(chankey);
    if (it != map_chankey2spe.end())
    {
      std::vector<Hit>& thisPMTHits = apair.second;
      for (Hit &ahit : thisPMTHits)
      {
        double hit_charge = ahit.GetCharge();
        double t_hitPE = hit_charge / map_chankey2spe.at(chankey);
        fHitPE.push_back(t_hitPE);
	double hit_time = ahit.GetTime();
       
	// Select events - beam on or off?
        if (BeamOn == 1 && t_hitPE > 5. && isExtTrig && hit_time<2000) { t_hitPE_ok = false;}   //only want nonCC events without prompt cluster,
        else if (BeamOn == 0 && beamflag != 1) {t_hitPE_ok = false; } //only want good beam off events *change to BeamOn==0 && !isLED ??????? TODO 
      }
    }
  }


  //  working w/ tank clusters  //
  get_ok = m_data->CStore.Get("ClusterMap", m_all_clusters);
  if (!get_ok) { Log("PhaseIINeutronBG tool: No ClusterMap object in ANNIEEvent! Abort!", v_error, verbosity); return true; }

  int n_cluster = 0; //NOT the number of clusters in evt; issa label
  int n_prompt_cluster = 0;
  bool hasPromptCluster = false;

  for (std::pair<double, std::vector<Hit>>&& cluster_pair : *m_all_clusters) //looping over all clusters in tank
  {  
    // reset variables
    bool isPrompt = false;
    bool isNonCC = false;
    bool pc_hitPE_ok = true;
    double cluster_charge = 0.;
    double cluster_time = cluster_pair.first;
    double cluster_PE = 0;
    std::vector<Hit> cluster_hits = cluster_pair.second;
    fClusterTime = -9999.;
    fClusterHitPE.clear();
    fClusterPMT.clear();
   
    if (cluster_time < 2000.) // this is a prompt cluster
    {
      isPrompt = true;  // this is for individual cluster; resets later
      hasPromptCluster = true;  // this is for the entire evt
      n_prompt_cluster += 1; 
    }

   
   
    double pmtcount = 0;
    // calculate cluster charge
    for (int i = 0; i < cluster_hits.size(); i++) //for each hit in this cluster
    {  
      bool lonely_PMT = false;  // this will not happen bc min hits to form clusters = 5
      if (cluster_hits.size() < 2) { lonely_PMT = true; std::cout << " [[ DEBUG ]] FOUND A LONELY PMT!" << std::endl; }

      int hit_ID = cluster_hits.at(i).GetTubeId(); // which PMT did this hit come from
      std::map<int, double>::iterator it = map_chankey2spe.find(hit_ID);
      fClusterPMT.push_back(hit_ID);
      if (it != map_chankey2spe.end()) // now find the PE for this hit given the PMT that saw the hit
      {
        double hit_charge = cluster_hits.at(i).GetCharge(); 
	if (hit_charge < 0) {continue;}
        double c_hitPE = hit_charge / map_chankey2spe.at(hit_ID); 

        fClusterHitPE.push_back(c_hitPE);
        cluster_charge += hit_charge;
        cluster_PE += c_hitPE;
        if (BeamOn == 1 && isPrompt && c_hitPE > pethresh) {pmtcount++;} //count number of PMTs that see PE > pethresh (do we consider prompt activity to be 1 PMT seeing >= 5pe, 2 PMTs seeing >= 3pe, etc)

      if (c_hitPE < 0) {pc_hitPE_ok = false;} //sometimes when the beam is funky it'll show up as a negative POT value
      }


      else
      {
        if (verbosity > 2) { std::cout << "PhaseIINeutronBG tool: FOUND A HIT FOR CHANKEY " << hit_ID << " BUT NO CONVERSION TO PE AVAILABLE. SKIPPING PE..." << std::endl; }
      }
    }


    if (BeamOn == 1 && pmtcount >= PMTthresh && isPrompt)
        {
	  pc_hitPE_ok = false;
	  noprompt = false;
	} //if prompt cluster with has >= PMTthresh PMTs with > pethresh, don't count

    if (verbosity > 3 && !pc_hitPE_ok) { std::cout << "PhaseIINeutronBG tool: Found a PROMPT cluster with >= " << PMTthresh << " PMTs with > " << pethresh << std::endl;}

    // get charge balance
    bool good_class = this->LoadTankClusterClassifiers(cluster_time);
    if (!good_class) { Log("PhaseIINeutronBG tool: NO cluster classifiers..", v_debug, verbosity); }


    fClusterTime = cluster_pair.first;
    fClusterNumber = n_cluster;
    fClusterCharge = cluster_charge;
    fClusterPE = cluster_PE;    
    fClusterHits = cluster_hits.size();

    // ALL cluster events
    h_clusterCharge->Fill(cluster_charge);
    h_clusterTime->Fill(cluster_time);
    h_clusterPE->Fill(cluster_PE);


    // basically saying "isNonCC" is true if trigger is non-cc AND no prompt activity that would indicate signal data
    goodnoncc = 0;

    if (isExtTrig && pc_hitPE_ok) 
    { 
      isNonCC = true; 
      goodnoncc = 1;
    } 

    // Updated version with filling up hist of beam off and beam on at the same time.

    if ((isBeam && BeamOn==1)) 
    {
      h_clusterTime_beam->Fill(cluster_time);
      h_clusterPE_beam->Fill(cluster_PE);
    }
    if ( beamflag==1 && isLED)
    {
	h_clusterTime_beamoff->Fill(cluster_time);
	h_clusterPE_beamoff->Fill(cluster_PE);
    }


    if (isBeam  && isNonCC && BeamOn==1)
    {
      h_clusterTime_nonCCbeam->Fill(cluster_time);
      h_clusterPE_nonCCbeam->Fill(cluster_PE);
    }
   

    if (cluster_time < 2000.0) // All prompt
      {
	h_clusterTime_prompt->Fill(cluster_time);
	h_clusterPE_prompt->Fill(cluster_PE);

	for (auto i = fClusterPMT.begin(); i != fClusterPMT.end(); ++i)
	  {
	    h_clusterTime_prompt_hitPMT->Fill(*i); // Filling PMT hit histogram for debugging
	  }
        

       	if (beamflag == 1 && isLED)
	  {
	    h_clusterTime_nonCCbeamoff_prompt->Fill(cluster_time); //Related to beam off readout in 2us
	    h_clusterPE_nonCCbeamoff_prompt->Fill(cluster_PE);
	  }

	if (isBeam && isNonCC && BeamOn == 1)
	  {
	    h_clusterTime_nonCCbeam_prompt->Fill(cluster_time);
	    h_clusterPE_nonCCbeam_prompt->Fill(cluster_PE);

	    if (!hasVeto)
	      {
		h_clusterTime_nonCCbeam_prompt_noVeto->Fill(cluster_time);
		h_clusterPE_nonCCbeam_prompt_noVeto->Fill(cluster_PE);

		if (!hasMRDTracks)
		  {
		    h_clusterTime_nonCCbeam_prompt_noVetoMRD->Fill(cluster_time);
		    h_clusterPE_nonCCbeam_prompt_noVetoMRD->Fill(cluster_PE);

		    if (BeamOn == 0 && isLED && fpot == 0 && fClusterCharge < 120 && fClusterChargeBalance < 0.4)
		      { 
			// Beam off uses LED trigger which is just 2us,
			// hence filling these histograms for cluster time < 2us if looking for beam off data
			h_clusterTime_background_neutrons->Fill(cluster_time);
			h_clusterPE_background_neutrons->Fill(cluster_PE);
		      }
		  }
	      }
	  }
      }
    else if (cluster_time > 2000.0) // All delayed
      {
	h_clusterTime_delayed->Fill(cluster_time); // Filling histograms for delayed clusters
	h_clusterPE_delayed->Fill(cluster_PE);

	for (auto i = fClusterPMT.begin(); i != fClusterPMT.end(); ++i)
	  {
	    h_clusterTime_delayed_hitPMT->Fill(*i); // Filling PMT hit histogram for debugging
	  }

       	if (beamflag == 1 && isLED)
	  {
	    h_clusterTime_nonCCbeamoff_delayed->Fill(cluster_time);
	    h_clusterPE_nonCCbeamoff_delayed->Fill(cluster_PE);
	  }

	if (isBeam && isNonCC && BeamOn == 1)
	  {
	    h_clusterTime_nonCCbeam_delayed->Fill(cluster_time);
	    h_clusterPE_nonCCbeam_delayed->Fill(cluster_PE);

	    if (!hasVeto)
	      {
		h_clusterTime_nonCCbeam_delayed_noVeto->Fill(cluster_time);
		h_clusterPE_nonCCbeam_delayed_noVeto->Fill(cluster_PE);

		if (!hasMRDTracks)
		  {
		    h_clusterTime_nonCCbeam_delayed_noVetoMRD->Fill(cluster_time);
		    h_clusterPE_nonCCbeam_delayed_noVetoMRD->Fill(cluster_PE);
                
		    if (BeamOn == 1 && cluster_time > 10000.0 && fClusterCharge < 120 && fClusterChargeBalance < 0.4)
		      {
			std::cout << "MADE IT HERE" << std::endl; // Debugging message
			h_clusterTime_background_neutrons->Fill(cluster_time);
			h_clusterPE_background_neutrons->Fill(cluster_PE);
		      }
		  }
	      }
	  }
      }
     
   
    if (cluster_time > 2000. && cluster_time < 4000.) //looking at different possible regions of interest, commented out bc wasn't needed but I left it in case it would be useful in the future. can probs just be deleted.
    {
	if (isBeam && BeamOn==1) {h_clusterCB_michel->Fill(fClusterChargeBalance);}
	if (beamflag==1 && BeamOn == 0 && isLED) {h_clusterCB_beamoff_michel->Fill(fClusterChargeBalance);}
        if ((isBeam && BeamOn==1) && isNonCC && !hasVeto && !hasMRDTracks) {h_clusterCB_michel_goodNonCC->Fill(fClusterChargeBalance);}
        if ((beamflag==1 && BeamOn == 0 && isLED) && isNonCC && !hasVeto && !hasMRDTracks) {h_clusterCB_michel_goodNonCC_beamoff->Fill(fClusterChargeBalance);}

    }

    if (cluster_time >= 4000. && cluster_time < 12000. && fpot==0){
	if (isBeam && BeamOn==1) {h_clusterCB_afterpulse->Fill(fClusterChargeBalance);}
	if ((isBeam && BeamOn==1) && isNonCC && !hasVeto && !hasMRDTracks) {h_clusterCB_afterpulse_goodNonCC->Fill(fClusterChargeBalance);}
	if (beamflag==1 && BeamOn == 0) {h_clusterCB_afterpulse_beamoff->Fill(fClusterChargeBalance);}
        if ((beamflag==1 && BeamOn == 0) && isNonCC && !hasVeto && !hasMRDTracks) {h_clusterCB_afterpulse_goodNonCC_beamoff->Fill(fClusterChargeBalance);}
    }

    if (cluster_time >= 12000. && cluster_time < 50000. && fpot==0){
	if (isBeam && BeamOn==1) {h_clusterCB_neutrons->Fill(fClusterChargeBalance);}
	if ((isBeam && BeamOn==1) && isNonCC && !hasVeto && !hasMRDTracks) {h_clusterCB_neutrons_goodNonCC->Fill(fClusterChargeBalance);}
	if (beamflag==1 && BeamOn == 0) {h_clusterCB_neutrons_beamoff->Fill(fClusterChargeBalance);}
        if ((beamflag==1 && BeamOn == 0) && isNonCC && !hasVeto && !hasMRDTracks) {h_clusterCB_neutrons_goodNonCC_beamoff->Fill(fClusterChargeBalance);}

    }

    if (cluster_time >= 50000. && fpot==0)
    {
	if (isBeam && BeamOn==1) {h_clusterCB_muons->Fill(fClusterChargeBalance);}
	if ((isBeam && BeamOn==1) && isNonCC && !hasVeto && !hasMRDTracks) {h_clusterCB_muons_goodNonCC->Fill(fClusterChargeBalance);}
	if (beamflag==1 && BeamOn == 0) {h_clusterCB_muons_beamoff->Fill(fClusterChargeBalance);}
        if ((beamflag==1 && BeamOn == 0) && isNonCC && !hasVeto && !hasMRDTracks) {h_clusterCB_muons_goodNonCC_beamoff->Fill(fClusterChargeBalance);}
    }

    t_TankCluster->Fill();
    n_cluster += 1;

  }
  if (isnoncc && noprompt) {goodnoncctrig==1;}
  t_Trigger->Fill(); // Want to look at number of triggers, not number of clusters for an event with that trigger. I'm not sure if this is actually working like I thought it would tho 
  return true;
}


bool PhaseIINeutronBG::Finalise(){
  std::cout << "Total POT: " << totalpot << std::endl;
  std::cout << "beamoff_lt: " << beamoff_lt << std::endl;
  std::cout << "woext_lt: " << woext_lt << std::endl;
  std::cout << "cc_lt: " << cc_lt << std::endl;
  std::cout << "nc_lt: " << nc_lt << std::endl;
  std::cout << "beamoff_pot: " << beamoff_pot << std::endl;
  std::cout << "woext_pot: " << woext_pot << std::endl;
  std::cout << "cc_pot: " << cc_pot << std::endl;
  std::cout << "nc_pot: " << nc_pot << std::endl;


  
  p2nbg_root_outp->cd();
  t_TankCluster->Write("", TObject::kOverwrite);
  t_Trigger->Write("", TObject::kOverwrite);
  //Live time and POT histograms for acculated values of respective event selection.
  h_totalpot->Fill(1,totalpot);
  h_beamoff_lt->Fill(1, beamoff_lt);
  h_beamoff_pot->Fill(1,beamoff_pot);
  h_woext_pot->Fill(1, woext_pot);
  h_woext_lt->Fill(1, woext_lt);
  h_cc_lt->Fill(1,cc_lt);
  h_cc_pot->Fill(1,cc_pot);
  h_nc_lt->Fill(1, nc_lt);
  h_nc_pot->Fill(1, nc_pot);
  this->WriteHist();
  p2nbg_root_outp->Close();
  delete p2nbg_root_outp;
  std::cout << "PhaseIINeutronBG tool exitting" << std::endl;
  if (verbosity > 5)
  {
    std::cout << " ...I'm going to go now...." << std::endl;
  }

  return true;
}

void PhaseIINeutronBG::InitTree()
{
  /***********************************
   ********** Define TTree ***********
   ***********************************/
  
  p2nbg_root_outp->cd();

  t_TankCluster = new TTree("phaseIITankClusterTree", "ANNIE Phase II Tank Cluster Tree");
  t_Trigger = new TTree("phaseIITriggerTree", "ANNIE Phase II Trigger Tree");
  t_TankCluster->Branch("runNumber", &fRunNumber, "runNumber/I");
  t_TankCluster->Branch("subrunNumber", &fSubrunNumber, "subrunNumber/I");
  t_TankCluster->Branch("runType", &fRunType, "runType/I");
  t_TankCluster->Branch("startTime", &fStartTime, "startTime/l");
  t_TankCluster->Branch("POT",&fpot,"POT/D"); 
  t_TankCluster->Branch("eventNumber", &fEventNumber, "eventNumber/I");
  t_TankCluster->Branch("eventTimeTank", &fEventTimeTank, "eventTimeTank/l");
  t_TankCluster->Branch("clusterNumber", &fClusterNumber, "clusterNumber/I");
  t_TankCluster->Branch("clusterCharge", &fClusterCharge, "clusterCharge/D");
  t_TankCluster->Branch("clusterTime", &fClusterTime, "clusterTime/D");
  t_TankCluster->Branch("clusterPE", &fClusterPE, "clusterPE/D");
  t_TankCluster->Branch("clusterMaxPE", &fClusterMaxPE, "clusterMaxPE/D");
  t_TankCluster->Branch("clusterHits", &fClusterHits, "clusterHits/I");
  t_TankCluster->Branch("clusterChargeBalance", &fClusterChargeBalance, "clusterChargeBalance/D");
  t_TankCluster->Branch("triggerWord", &trigword, "triggerWord/I");
  t_TankCluster->Branch("triggerWordExtended", &trigext, "triggerWordExtended/I");
  t_TankCluster->Branch("goodNonCC", &goodnoncc, "goodNonCC/I");
  t_TankCluster->Branch("vetoHit", &fVetoHit, "vetoHit/I");
  t_TankCluster->Branch("nMRDTracks", &fnMRDTracks, "nMRDTracks/I");
  t_TankCluster->Branch("hitPE", &fHitPE);
  t_TankCluster->Branch("hitPMT",&fHitPMT);
  t_TankCluster->Branch("cluster_hitPE", &fClusterHitPE);

  t_Trigger->Branch("goodNonCCTrig",&goodnoncctrig,"goodNonCCTrig/I");
  t_Trigger->Branch("triggerWord",&trigword, "triggerWord/I");
  t_Trigger->Branch("triggerWordExtended",&trigext, "triggerWordExtended/I");
  gROOT->cd();
}

void PhaseIINeutronBG::InitHist() // initialize all the histograms declared in the header file
{
  /***********************************
   ******** Define Histograms ********
   ***********************************/

  p2nbg_root_outp->cd();

  h_clusterCharge = new TH1F("h_clusterCharge", "clusterCharge", 1000, 0, 1);
  h_beamFlag = new TH1F("h_beamFlag", "0=good, 1=off, 2=bad, 3=missing, -1=fill error", 5, -1, 4);
  h_clusterTime_beam = new TH1F("h_clusterTime_beam", "clusterTime for BEAM evts", 250, 0, 71000);
  h_clusterTime_beamoff = new TH1F("h_clusterTime_beamoff", "clusterTime for BEAM off evts", 250, 0, 71000);

  h_clusterTime = new TH1F("h_clusterTime", "clusterTime", 250, 0, 71000);
  
  h_beamoff_lt = new TH1F("h_beamoff_lt", "beam off live time", 3,0,3);
  h_woext_lt = new TH1F("h_woext_lt", "without extended trigger live time",3,0,3);
  h_cc_lt = new TH1F("h_cc_lt", "cc live time",3,0,3);
  h_nc_lt = new TH1F("h_nc_lt", "nc live time",3,0,3);

  h_totalpot = new TH1F("h_totalpot", "total POT", 3,0,3);


  h_beamoff_pot = new TH1F("h_beamoff_pot", "beam off POT", 3,0,3);
  h_woext_pot = new TH1F("h_woext_pot", "without extended trigger POT",3,0,3);
  h_cc_pot = new TH1F("h_cc_pot", "cc POT",3,0,3);
  h_nc_pot = new TH1F("h_nc_pot", "nc POT",3,0,3);
  // h_clusterTime_beam = new TH1F("h_clusterTime_beam", "clusterTime for BEAM evts", 250, 0, 71000);
  h_clusterTime_prompt = new TH1F("h_clusterTime_prompt", "clusterTime for ALL PROMPT evts", 200, 0, 2100);
  h_clusterTime_delayed = new TH1F("h_clusterTime_delayed", "clusterTime for ALL DELAYED evts", 250, 0, 71000);
  h_clusterTime_nonCCbeam = new TH1F("h_clusterTime_nonCCbeam", "clusterTime for NON-CC BEAM evts", 250, 0, 71000);
  h_clusterTime_nonCCbeamoff = new TH1F("h_clusterTime_nonCCbeamoff", "clusterTime for NON-CC BEAM-off evts", 250, 0, 71000);
  
  h_clusterTime_nonCCbeam_prompt = new TH1F("h_clusterTime_nonCCbeam_prompt", "clusterTime for NON-CC BEAM evts - prompt", 200, 0, 2100);
  h_clusterTime_nonCCbeamoff_prompt = new TH1F("h_clusterTime_nonCCbeamoff_prompt", "clusterTime for NON-CC BEAM off evts - prompt", 200, 0, 2100);

  h_clusterTime_nonCCbeam_prompt_noVeto = new TH1F("h_clusterTime_nonCCbeam_prompt_noVeto", "clusterTime for NON-CC BEAM evts w NO VETO - prompt", 200, 0, 2100);
  h_clusterTime_nonCCbeam_prompt_noVetoMRD = new TH1F("h_clusterTime_nonCCbeam_prompt_noVetoMRD", "clusterTime for NON-CC BEAM evts w NO VETO OR MRD - prompt", 200, 0, 2100);
  h_clusterTime_nonCCbeam_delayed = new TH1F("h_clusterTime_nonCCbeam_delayed", "clusterTime for NON-CC BEAM evts - delayed", 250, 0, 71000);
  h_clusterTime_nonCCbeamoff_delayed = new TH1F("h_clusterTime_nonCCbeamoff_delayed", "clusterTime for NON-CC BEAM off evts - delayed", 250, 0, 71000);

  h_clusterTime_nonCCbeam_delayed_noVeto = new TH1F("h_clusterTime_nonCCbeam_delayed_noVeto", "clusterTime for NON-CC BEAM evts w NO VETO - delayed", 250, 0, 71000);
  h_clusterTime_nonCCbeam_delayed_noVetoMRD = new TH1F("h_clusterTime_nonCCbeam_delayed_noVetoMRD", "clusterTime for NON-CC BEAM evts w NO VETO OR MRD - delayed",250,0,71000);
  h_clusterTime_background_neutrons = new TH1F("h_clusterTime_background_neutrons", "Background neutrons cluster time",250,0,71000);
  h_clusterTime_hitPMT = new TH1F("h_clusterTime_hitPMT","Tube IDs for cluster hits from 6000-8000ns",200,300,500);
  h_clusterTime_delayed_hitPMT = new TH1F("h_clusterTime_delayed_hitPMT","Tube IDs for all cluster hits in delayed window", 200,300,500);
  h_clusterTime_prompt_hitPMT = new TH1F("h_clusterTime_prompt_hitPMT","Tube IDs for all cluster hits in prompt window",200,300,500);
  h_clusterPE = new TH1F("h_clusterPE", "clusterPE", 300, 0, 300);
  h_clusterPE_beam = new TH1F("h_clusterPE_beam", "clusterPE for BEAM evts", 300, 0, 300);
  h_clusterPE_beamoff = new TH1F("h_clusterPE_beamoff", "clusterPE for BEAM off evts", 300, 0, 300);

  h_clusterPE_prompt = new TH1F("h_clusterPE_prompt", "clusterPE for ALL PROMPT evts", 300, 0, 300);
  h_clusterPE_delayed = new TH1F("h_clusterPE_delayed", "clusterPE for ALL DELAYED evts", 300, 0, 300);
  h_clusterPE_nonCCbeam = new TH1F("h_clusterPE_nonCCbeam", "clusterPE for NON-CC BEAM evts", 300, 0, 300);
  h_clusterPE_nonCCbeamoff = new TH1F("h_clusterPE_nonCCbeamoff", "clusterPE for NON-CC BEAM-off evts", 300, 0, 300);

  h_clusterPE_nonCCbeam_prompt = new TH1F("h_clusterPE_nonCCbeam_prompt", "clusterPE for NON-CC BEAM evts - prompt", 300, 0, 300);
  h_clusterPE_nonCCbeamoff_prompt = new TH1F("h_clusterPE_nonCCbeamoff_prompt", "clusterPE for NON-CC BEAM off evts - prompt", 300, 0, 300);

  h_clusterPE_nonCCbeam_prompt_noVeto = new TH1F("h_clusterPE_nonCCbeam_prompt_noVeto", "clusterPE for NON-CC BEAM evts w NO VETO - prompt", 300, 0, 300);
  h_clusterPE_nonCCbeam_prompt_noVetoMRD = new TH1F("h_clusterPE_nonCCbeam_prompt_noVetoMRD","clusterPE for NON-CC BEAM evts w NO VETO OR MRD - prompt", 300, 0, 300);
  h_clusterPE_nonCCbeam_delayed = new TH1F("h_clusterPE_nonCCbeam_delayed", "clusterPE for NON-CC BEAM evts - delayed", 300, 0, 300);
  h_clusterPE_nonCCbeamoff_delayed = new TH1F("h_clusterPE_nonCCbeamoff_delayed", "clusterPE for NON-CC BEAM off evts - delayed", 300, 0, 300);

  h_clusterPE_nonCCbeam_delayed_noVeto = new TH1F("h_clusterPE_nonCCbeam_delayed_noVeto", "clusterPE for NON-CC BEAM evts w NO VETO - delayed", 300, 0, 300);
  h_clusterPE_nonCCbeam_delayed_noVetoMRD = new TH1F("h_clusterPE_nonCCbeam_delayed_noVetoMRD", "clusterPE for NON-CC BEAM evts w NO VETO OR MRD - delayed", 300, 0, 300);
  h_clusterPE_background_neutrons = new TH1F("h_clusterPE_background_neutrons","cluster PE for background neutrons",300,0,300);
  h_clusterCB_michel = new TH1F("h_clusterCB_michel","",100,0,1);
  h_clusterCB_michel_goodNonCC = new TH1F("h_clusterCB_michel_goodNonCC","",100,0,1);
  h_clusterCB_afterpulse = new TH1F("h_clusterCB_afterpulse","",100,0,1);
  h_clusterCB_afterpulse_goodNonCC = new TH1F("h_clusterCB_afterpulse_goodNonCC","",100,0,1);
  h_clusterCB_neutrons = new TH1F("h_clusterCB_neutrons","",100,0,1);
  h_clusterCB_neutrons_goodNonCC = new TH1F("h_clusterCB_neutrons_goodNonCC","",100,0,1);
  h_clusterCB_muons = new TH1F("h_clusterCB_muons","",100,0,1);
  h_clusterCB_muons_goodNonCC = new TH1F("h_clusterCB_muons_goodNonCC","",100,0,1);
  
  h_clusterCB_beamoff_michel = new TH1F("h_clusterCB_beamoff_michel","",100,0,1);
  h_clusterCB_michel_goodNonCC_beamoff = new TH1F("h_clusterCB_michel_goodNonCC_beamoff","",100,0,1);
  h_clusterCB_afterpulse_beamoff = new TH1F("h_clusterCB_afterpulse_beamoff","",100,0,1);
  h_clusterCB_afterpulse_goodNonCC_beamoff = new TH1F("h_clusterCB_afterpulse_goodNonCC_beamoff","",100,0,1);
  h_clusterCB_neutrons_beamoff = new TH1F("h_clusterCB_neutrons_beamoff","",100,0,1);
  h_clusterCB_neutrons_goodNonCC_beamoff = new TH1F("h_clusterCB_neutrons_goodNonCC_beamoff","",100,0,1);
  h_clusterCB_muons_beamoff = new TH1F("h_clusterCB_muons_beamoff","",100,0,1);
  h_clusterCB_muons_goodNonCC_beamoff = new TH1F("h_clusterCB_muons_goodNonCC_beamoff","",100,0,1);

  gROOT->cd();
}

void PhaseIINeutronBG::WriteHist()
{
  /***********************************
   ******** Write Histograms *********
   ***********************************/

  p2nbg_root_outp->cd();

  TDirectory *dir_allhist = p2nbg_root_outp->mkdir("Histograms");
  dir_allhist->cd();

  h_clusterCharge->Write();
  h_beamFlag->Write();
  h_beamoff_lt->Write();
  h_woext_lt->Write();
  h_cc_lt->Write();
  h_nc_lt->Write();

  h_totalpot->Write();
  h_beamoff_pot->Write();
  h_woext_pot->Write();
  h_cc_pot->Write();
  h_nc_pot->Write();
  
  h_clusterTime->Write();
  h_clusterTime_beamoff->Write();

  h_clusterTime_beam->Write();
  h_clusterTime_prompt->Write();
  h_clusterTime_delayed->Write();
  h_clusterTime_nonCCbeam->Write();
  h_clusterTime_nonCCbeamoff->Write();


  h_clusterTime_nonCCbeam_prompt->Write();
  h_clusterTime_nonCCbeamoff_prompt->Write();

  h_clusterTime_nonCCbeam_prompt_noVeto->Write();
  h_clusterTime_nonCCbeam_prompt_noVetoMRD->Write();
  h_clusterTime_nonCCbeam_delayed->Write();
  h_clusterTime_nonCCbeamoff_delayed->Write();

  h_clusterTime_nonCCbeam_delayed_noVeto->Write();
  h_clusterTime_nonCCbeam_delayed_noVetoMRD->Write();
  h_clusterTime_background_neutrons->Write();
  h_clusterTime_hitPMT->Write();
  h_clusterTime_delayed_hitPMT->Write();
  h_clusterTime_prompt_hitPMT->Write();
  h_clusterPE->Write();
  h_clusterPE_beamoff->Write();

  h_clusterPE_beam->Write();
  h_clusterPE_prompt->Write();
  h_clusterPE_delayed->Write();
  h_clusterPE_nonCCbeam->Write();
  h_clusterPE_nonCCbeamoff->Write();

  h_clusterPE_nonCCbeam_prompt->Write();
  h_clusterPE_nonCCbeamoff_prompt->Write();

  h_clusterPE_nonCCbeam_prompt_noVeto->Write();
  h_clusterPE_nonCCbeam_prompt_noVetoMRD->Write();
  h_clusterPE_nonCCbeam_delayed->Write();
  h_clusterPE_nonCCbeamoff_delayed->Write();

  h_clusterPE_nonCCbeam_delayed_noVeto->Write();
  h_clusterPE_nonCCbeam_delayed_noVetoMRD->Write();
  h_clusterPE_background_neutrons->Write();
  h_clusterCB_michel->Write();
  h_clusterCB_michel_goodNonCC->Write();
  h_clusterCB_afterpulse->Write();
  h_clusterCB_afterpulse_goodNonCC->Write();
  h_clusterCB_neutrons->Write();
  h_clusterCB_neutrons_goodNonCC->Write();
  h_clusterCB_muons->Write();
  h_clusterCB_muons_goodNonCC->Write();

  h_clusterCB_beamoff_michel->Write();
  h_clusterCB_michel_goodNonCC_beamoff->Write();
  h_clusterCB_afterpulse_beamoff->Write();
  h_clusterCB_afterpulse_goodNonCC_beamoff->Write();
  h_clusterCB_neutrons_beamoff->Write();
  h_clusterCB_neutrons_goodNonCC_beamoff->Write();
  h_clusterCB_muons_beamoff->Write();
  h_clusterCB_muons_goodNonCC_beamoff->Write();  

  gROOT->cd();
}

bool PhaseIINeutronBG::LoadTankClusterClassifiers(double cluster_time)
{
  bool got_ccp = m_data->Stores["ANNIEEvent"]->Get("ClusterChargePoints", cluster_CP);
  bool got_ccb = m_data->Stores["ANNIEEvent"]->Get("ClusterChargeBalances", cluster_CB);
  bool got_cmpe = m_data->Stores["ANNIEEvent"]->Get("ClusterMaxPEs", cluster_maxPEs);
  bool good_class = got_ccp && got_ccb && got_cmpe;
  if (!good_class) { Log("PhaseIINeutronBG tool: One of the charge cluster classifiers is not available", v_debug, verbosity); }
  else
  {
    Log("PhaseIINeutronBG tool: Setting fCluster variables to classifier parameters", v_debug, verbosity);
    fClusterMaxPE = cluster_maxPEs.at(cluster_time);
    fClusterChargeBalance = cluster_CB.at(cluster_time);
  }

  return good_class;
}
