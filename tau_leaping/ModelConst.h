#pragma once

#include <memory>
#include "defines.h"

class CSRMatrixPacked; //forward declaration
class CSCMatrix; //forward declaration

struct ModelConst
{
	ModelConst(int channels_M);

	//matrix that keeps reactions' description
	std::auto_ptr<CSRMatrixPacked> m_reactionsDesc;

	//reactions' constants
	vec_float_t m_reactionsConstants;

	//reactions' orders [1,2,3] 
	vec_char_t m_ReactionOrder;

	//indices of first and second reactants in reactions
	vec_int_t m_s0, m_s1;

	void Print()const;

	void Init(CSCMatrix const *reactionsDescCSC);

	void CopyTo(ModelConst &modelConst)const;
	
};