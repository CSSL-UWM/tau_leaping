#pragma once

#include <thrust/device_vector.h>

struct ModelState
{
	//the amount of each reactants at this very moment of time
	thrust::device_vector<int> m_reactantsPopulation;

	//reactions' propensities
	thrust::device_vector<float> m_reactionsPropensity;

	//partial sums of propensities
	thrust::device_vector<float> m_rPropPartSum;

	//partial sums of propensities for non-critical reactions
	thrust::device_vector<float> m_rPropNonCritPartSum;//TODO: it seems we need critical instead...

	//array with trues for critical reactions and falses for non-critical ones
	thrust::device_vector<bool> m_criticalClassifier;
	
	//vec_int_t m_sortedIndices;

	//seeds for producing RNs
	thrust::device_vector<uint4> m_seeds;

	thrust::device_vector<int> m_k;//how many times we fire every reaction
public:
	ModelState(int channels_M, int reactants_N);
	~ModelState(void);
};
