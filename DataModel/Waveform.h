/* vim:set noexpandtab tabstop=4 wrap */
#ifndef WAVEFORMCLASS_H
#define WAVEFORMCLASS_H

#include <SerialisableObject.h>
#include <iostream>

using namespace std;

template <class T>
class Waveform : public SerialisableObject{

	friend class boost::serialization::access;

	public:
	Waveform() : fStartTime(), fSamples(std::vector<T>{}) {serialise=true;}
	Waveform(double tsin, std::vector<T> samplesin) : fStartTime(tsin), fSamples(samplesin){serialise=true;}
	virtual ~Waveform(){}

	inline double GetStartTime() const {return fStartTime;}
	inline std::vector<T>* GetSamples() {return &fSamples;}
	inline const std::vector<T>& Samples() const { return fSamples; }
	inline T GetSample(int i) const {return fSamples.at(i);}

	inline void SetStartTime(double ts) {fStartTime=ts;}
	inline void SetSamples(std::vector<T> sam) {fSamples=sam;}
	inline void PushSample(T asample) {fSamples.push_back(asample);}
	inline void ClearSamples() {fSamples.clear();}

	bool Print() {
		int verbose=0;
		cout<<"StartTime : "<<fStartTime<<endl;
		cout<<"NSamples : "<<fSamples.size()<<endl;
		if(verbose){
			cout<<"Samples : {";
			for(int samplei=0; samplei<fSamples.size(); samplei++){
				cout<<fSamples.at(samplei);
				if((samplei+1)!=fSamples.size()) cout<<", ";
			}
			cout<<"}"<<endl;
		}

		return true;
	}

	protected:
	double fStartTime;
	std::vector<T> fSamples;

	template<class Archive> void serialize(Archive & ar, const unsigned int version){
		if(serialise){
			ar & fStartTime;
			ar & fSamples;
		}
	}
};

template <class T>
class MCWaveform : public Waveform<T> {

        friend class boost::serialization::access;
	
        public:
        MCWaveform() : Waveform<T>(), fParents(std::vector<std::vector<int>>{}) {}
        MCWaveform(double tsin, std::vector<T> samplesin, std::vector<std::vector<int>> theparents) : Waveform<T>(tsin, samplesin), fParents(theparents) {}
	virtual ~MCWaveform(){};

	inline const std::vector<std::vector<int>>* GetParents() { return &fParents; }
	inline const std::vector<std::vector<int>>& Parents()    { return fParents; }
	inline std::vector<int> ParentsAtSample(int sample)      { return fParents[sample]; }

	inline void SetParents(std::vector<std::vector<int>> parentsin){ fParents = parentsin; }

	inline Waveform<T> GetBaseWaveform(){ return *this; }

	// Override the base print to add in parent info
	bool Print() {
	  int verbose=0;
	  cout << "StartTime : " << this->fStartTime << endl;
	  cout << "NSamples : " << this->fSamples.size() << endl;
	  if(verbose){
	    cout << "Samples (parents) : {";
	    for(int samplei=0; samplei < this->fSamples.size(); samplei++){
	      cout << this->fSamples.at(samplei);

	      if (this->fParents.size() == this->fSamples.size()) {
		cout << "(";
		for (auto parent : this->ParentsAtSample(samplei)) {
		  cout << parent;
		  if ((samplei+1) != this->ParentsAtSample(samplei).size())
		    cout << ", ";
		  else
		    cout << ")";
		}
	      }
	      
	      if ((samplei+1) != this ->fSamples.size()) cout << ", ";
	    }
	    cout<<"}"<<endl;
	  }
	  
	  return true;
	}

	protected:
	// The outer vector has the same index as the samples.
	// The inner vector contains the MCParticle index of the parents contibuting to that sample
	std::vector<std::vector<int>> fParents;

};

#endif
