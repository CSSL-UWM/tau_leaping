#include "CSCHostMatrixData.h"
#include "StochKitReader.h"
#include <thrust/sort.h>
#include <iomanip>

using std::vector;
using std::unordered_map;
using std::string;
using std::auto_ptr;
using thrust::host_vector;
using std::setw;



//we have information about the reaction set encoded 
//in reactionProps and reactantsPopulations containers
//And we need to produce CSCMatrix out of them
//Reactions (channels) are in columns within CSCHostMatrixData...
void ReadXMLToMatrixStruct(vector<ReactionProp> const &reactionProps, 
						   unordered_map<string, std::pair<int, int> > const &reactantsPopulations, 
						   CSCHostMatrixData &mData,
						   vector<string> &reactantNames,
						   host_vector<int> &reactantsPopulation, 
						   host_vector<float> &reactionConstants)
{
	typedef unordered_map<string, std::pair<int, int> > reactantsPopulations_t;
	mData.colPtr.resize(reactionProps.size()+1,0);

	reactantNames.resize(reactantsPopulations.size());
	reactantsPopulation.resize(reactantsPopulations.size());
	for(reactantsPopulations_t::const_iterator cit=reactantsPopulations.begin();
		cit!=reactantsPopulations.end(); ++cit)
	{	
		 reactantNames[cit->second.first]=cit->first;
		 reactantsPopulation[cit->second.first]=cit->second.second;
	}

	host_vector<int> rInd;

	int col=0;
	for(vector<ReactionProp>::const_iterator cit=reactionProps.begin();
		cit!=reactionProps.end(); ++cit)
	{
		int colLen=cit->reactants.size()+cit->products.size();//how many products+reactants are involved in this reaction
		mData.colPtr[col+1]=mData.colPtr[col]+colLen;

		rInd.clear();

		reactionConstants[col]=cit->rate;

		//unordered_map<string, char>::const_iterator
		for(auto it=cit->reactants.begin();
			it!=cit->reactants.end(); ++it)
		{
			reactantsPopulations_t::const_iterator rit=reactantsPopulations.find(it->first);
			assert(rit!=reactantsPopulations.end());
			//int i=reactantsPopulations.find(it->first).first;
			assert(rit->second.first<(int)reactantNames.size()&&rit->second.first>=0);
			rInd.push_back(rit->second.first);
			//rInd.push_back(reactantsPopulations[it->first].first);
		}

		for(auto it=cit->products.begin();
			it!=cit->products.end(); ++it)
		{
			reactantsPopulations_t::const_iterator rit=reactantsPopulations.find(it->first);
			assert(rit->second.first<reactantNames.size()&&rit->second.first>=0);
			//rInd.push_back(reactantsPopulations[it->first].first);
			rInd.push_back(rit->second.first);
		}

		thrust::sort(rInd.begin(), rInd.end());

		unordered_map<string, bool> foundInReactants;
		//ReactionProp const &rp=cit;
		for(size_t i=0; i<rInd.size(); ++i)
		{
			int ind=rInd[i];
			const string &rName=reactantNames[rInd[i]];
			auto it=cit->reactants.find(rName);
			char v=0;
			if(it==cit->reactants.end()||foundInReactants[rName])
			{
				it=cit->products.find(rName);
				if(it!=cit->products.end())
				{
					v=it->second;
				}else
				{
					assert(false);
				}
			}else
			{
			//	assert(cit->products.find(rName)==cit->products.end());
				v=-it->second;
				foundInReactants[rName]=true;
			}

			

			mData.rates.push_back(v);
			mData.rowInd.push_back(rInd[i]);
		}

		++col;
	}
}


std::ostream &operator<<(std::ostream &os, CSCHostMatrixData const &data)
{
	int channelsCount=data.colPtr.size()-1;
	for(int i=0; i<channelsCount; ++i)
	{
		os<<"reaction "<<setw(3)<<i+1<<";  ";
		int bi=data.colPtr[i];
		int ei=data.colPtr[i+1];
		for(int ind=bi; ind<ei; ++ind)
		{
			os<<"S"<<setw(2)<<data.rowInd[ind]+1<<":"<<setw(2)<<(int)data.rates[ind]<<"  ";
		}
		os<<"\n";
	}
	return os;
}