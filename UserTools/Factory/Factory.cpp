#include "Factory.h"

Tool* Factory(std::string tool) {
    Tool* ret = 0;
    if (tool=="PlotWaveforms") ret=new PlotWaveforms;
    if (tool=="LoadGeometry") ret=new LoadGeometry;
    if (tool=="LoadWCSim") ret=new LoadWCSim;
    if (tool=="PMTWaveformSim") ret=new PMTWaveformSim;
    if (tool=="PhaseIIADCHitFinder") ret=new PhaseIIADCHitFinder;

    return ret;
}
