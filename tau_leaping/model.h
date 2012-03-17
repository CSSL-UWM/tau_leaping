#ifndef __MODEL_H__
#define __MODEL_H__


#include "defines.h"
#include <fstream>
#include <string>
#include <map>

#include "ModelConst.h"
#include "ModelState.h"
#include "statistics.h"

class CSRMatrixPacked; //forward declaration
class CSCMatrix; //forward declaration


//class stores all data we need for computation and 
//provides appropriate functions
class Model
{
private:

	ModelConst m_ConstData;
	ModelState m_StateData;
	
	//Number of reactions
	const int m_channelsCount;

	//Number of reactants
	const int m_reactantsCount;

	//an empirical variable for choosing direct/tau-leaping branch
	static const int m_N=10;

	//current time
	float m_tau;

	
public:
	Model(int channels_M, int reactants_N);
	~Model();

	int ChannelsCount()const{return m_channelsCount;}
	int ReactantsCount()const{return m_reactantsCount;}
	//void FillRandom();

	void InitSeeds();
	
	float PerformCalc(float endTau, float eps, std::vector<char> const &r2Track,  
		stat_t &statistics, int maxIterations);
	static __host__ Model* ReadXML(const char *fn, 
		std::vector<std::string> &reactantNames);

	Model *Clone()const;
	
private:
	
	//compute propensities
	void CalcPropensities();

	//invoke direct method (once)
	float DirectMethodStep();

	void TransformMatrix(CSCMatrix const *reactionsDescCSC);

	//initialize seeds for RNs 
	__host__ void CreateRandomSeedsUint4();

	//log statistics into a file
	void LogStatistics(std::ofstream &logF, std::vector<char> const &r2Track, int nIter, stat_t &statistics)const;
	
	//perform one step of the algorithm, returns time step
	float Step(float eps);

	//detects tau-leaping time step
	float TauLeapingTimeStep(float eps, float a0);

	//detects a crtitical reaction to be fired
	void GetReaction2Fire(vec_float_t const &rPropPartSum, float &tauCrit, int &k2Fire)const;

	//set a number to be fired for a reaction
#ifdef _DEBUG
	thrust::pair<vec_int_t::const_iterator, vec_int_t::const_iterator>
#else
	void 
#endif
		SetNumberOfFiring(bool fireCritical, int k2Fire, float tau);

	//fire reactions
	void FireReactions();

private:
	//disable copying
	Model(const Model &);
	Model operator =(const Model &);

};

#endif