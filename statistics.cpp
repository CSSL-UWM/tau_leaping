#include "statistics.h"
#include <cassert>
#include <numeric>
#include <fstream>

using std::vector;
using std::map;
using std::pair;
using std::accumulate;
using std::make_pair;
using std::ofstream;
using std::endl;

void AllStat2StepStat(vector<stat_t> const &allStat, StepStat_t &stepStat, double timeStep){
	int reactantsToTrack=allStat[0].begin()->second.size();
	//stepStat.resize(reactantsToTrack);
	int realizationCount=allStat.size();

	double currTime=0;
	/*for(int i=0; i<reactantsToTrack; ++i){
		vector<vector<double>> &currTimeStat=stepStat[0];
		currTimeStat.resize(reactantsToTrack);
		for(int reactant=0; reactant<reactantsToTrack; ++reactant){
			currTimeStat[reactant].resize(realizationCount, 0);
			//currTimeStat[reactant].push_back(0);
			//std::fill(currTimeStat[reactant].begin(),
		}
	}*/
	//currTime+=timeStep;
	for(;;currTime+=timeStep){
		vector<vector<double>> &currTimeStat=stepStat[currTime];
		currTimeStat.resize(reactantsToTrack);

		for(int reactant=0; reactant<reactantsToTrack; ++reactant){
			currTimeStat[reactant].resize(realizationCount);
		}
		for(int realization=0; realization<realizationCount; ++realization){
			stat_t const &realizationAllStat=allStat[realization];
			auto it=realizationAllStat.upper_bound(currTime);
			if(it==realizationAllStat.end()){
				stepStat.erase(currTime);
				return;
			}
			auto prevIt=it;
			--prevIt;

			double time=it->first;
			double prevTime=prevIt->first;
			double shift=(currTime-prevTime)/(time-prevTime);
			assert(0<=shift&&shift<=1);

			for(int reactant=0; reactant<reactantsToTrack; ++reactant){
				double population=prevIt->second[reactant]+(it->second[reactant]-prevIt->second[reactant])*shift;
				currTimeStat[reactant][realization]=population;
			}
		}
	}
	

	/*int minEndTime=
	for(int i=0; i<realizationCount; ++i){
	}*/
}


void StepStat2StdDev(StepStat_t const &stepStat, StdDevStat_t &meanStdDev){
	int reactantsToTrack=stepStat.cbegin()->second.size();
	for(auto it=stepStat.cbegin(); it!=stepStat.cend(); ++it){
		vector<vector<double>> const &timeSlice=it->second;
		vector<pair<double, double>> &stdDevTimeSlice=meanStdDev[it->first];
		stdDevTimeSlice.resize(reactantsToTrack);
		for(int reactant=0; reactant<reactantsToTrack; ++reactant){
			vector<double> const &data=timeSlice[reactant];
			double avg=accumulate(data.begin(), data.end(),0.);
			avg/=double(data.size());

			double stdDev=0.;
			for(size_t i=0; i<data.size(); ++i){
				stdDev+=(avg-data[i])*(avg-data[i]);
			}
			stdDev=sqrt(stdDev/double(data.size()-1));
			stdDevTimeSlice[reactant]=make_pair(avg, stdDev);
		}
	}

}

void StatToFile(std::string fileName, StdDevStat_t const &meanStdDev)
{
    ofstream ofileMean("mean_"+fileName);
    assert(ofileMean.is_open());
    for(auto it=meanStdDev.begin(); it!=meanStdDev.end(); ++it){
    	ofileMean<<it->first<<" ";
    	//std::copy(it->second.begin(), it->second.end(), std::ostream_iterator<double>(ofile));
    	for(size_t i=0; i<it->second.size();++i){
    		ofileMean<<it->second[i].first<<" ";
    	}
    	ofileMean<<endl;
    }

	ofstream ofileStdDev("stddev_"+fileName);
    assert(ofileStdDev.is_open());
    for(auto it=meanStdDev.begin(); it!=meanStdDev.end(); ++it){
    	ofileStdDev<<it->first<<" ";
    	for(size_t i=0; i<it->second.size();++i){
    		ofileStdDev<<it->second[i].second<<" ";
    	}
    	ofileStdDev<<endl;
    }
}