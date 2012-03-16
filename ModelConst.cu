#include "ModelConst.h"
#include "CSRMatrixPacked.h"
#include "CSCMatrix.h"
#include "CSRMatrix.h"
#include "functions.h"

#include <cassert>

ModelConst::ModelConst(int channels_M):
		m_ReactionOrder(channels_M, -1),
		m_reactionsConstants(channels_M),
		m_s0(channels_M),
	    m_s1(channels_M)
{
		
}

void ModelConst::Init(CSCMatrix const *reactionsDescCSC)
{
	reactionsDescCSC->DetectReactionOrdersAndFillReactantIndices(m_ReactionOrder, m_s0, m_s1);
	m_reactionsDesc=reactionsDescCSC->ToCSR()->ToCSRPacked();

}

void ModelConst::Print()const
{
	printf("reactions' ordres\n");
	print_vector<char, int>(m_ReactionOrder);
	printf("\n");

	printf("reaction constatns: \n");
	print_vector<float, float> (m_reactionsConstants);
	printf("\n");

	printf("s0:\n");
	print_vector<int, int>(m_s0);
	printf("\ns1:\n");
	print_vector<int, int>(m_s1);
	printf("\n");
}

void ModelConst::CopyTo(ModelConst &modelConst)const
{
	CSRMatrixPacked const *ptr=m_reactionsDesc.get();
	assert(ptr);
	modelConst.m_reactionsDesc= std::auto_ptr<CSRMatrixPacked>(new CSRMatrixPacked(*ptr));

	modelConst.m_reactionsConstants=m_reactionsConstants;
	modelConst.m_ReactionOrder=m_ReactionOrder;
	modelConst.m_s0=m_s0;
	modelConst.m_s1=m_s1;
}