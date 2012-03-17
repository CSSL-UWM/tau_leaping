#include "ModelState.h"

ModelState::ModelState(int channels_M, int reactants_N):
	m_reactantsPopulation(reactants_N),
	m_reactionsPropensity(channels_M),
	m_criticalClassifier(channels_M),
	m_rPropPartSum(channels_M),
    m_rPropNonCritPartSum(channels_M),
	m_k(channels_M),
	m_seeds(channels_M*2)
{
}

ModelState::~ModelState(void)
{
}
