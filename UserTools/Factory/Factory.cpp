#include "Factory.h"

Tool* Factory(std::string tool) {
    Tool* ret = 0;
    if (tool=="PlotWaveforms") ret=new PlotWaveforms;
    if (tool=="LoadGeometry") ret=new LoadGeometry;
    if (tool=="LoadWCSim") ret=new LoadWCSim;
    if (tool=="PMTWaveformSim") ret=new PMTWaveformSim;
    if (tool=="PhaseIIADCHitFinder") ret=new PhaseIIADCHitFinder;
    if (tool=="ClusterFinder") ret=new ClusterFinder;
    if (tool=="BackTracker") ret=new BackTracker;
    if (tool=="ClusterClassifiers") ret=new ClusterClassifiers;
    if (tool=="ClusterSelector") ret=new ClusterSelector;
    if (tool=="SelectionEffnPurity") ret=new SelectionEffnPurity;

    return ret;
}
