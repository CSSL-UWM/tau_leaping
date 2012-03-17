#include "model.h"
#include "StochKitReader.h"
#include <new>
#include <thrust/host_vector.h>
#include <thrust/sort.h>
#include "CSCMatrix.h"
#include "CSCHostMatrixData.h"

using std::cout;
using std::endl;
using std::string;
using std::vector;
using std::unordered_map;

using thrust::host_vector;

Model* Model::ReadXML(const char *fn, std::vector<std::string> &reactantNames)
{
	cout<<"reading xml StochKit model... ";
	StochKitReader skr;
	vector<ReactionProp> reactionProps;
	typedef unordered_map<string, std::pair<int, int> > reactantsPopulations_t;
	reactantsPopulations_t reactantsPopulations;
	skr.ReadXML(fn, reactionProps, reactantsPopulations);

	CSCHostMatrixData mData;
	host_vector<float> reactionConstants(reactionProps.size());
	host_vector<int> reactantsPopulation;
	ReadXMLToMatrixStruct(reactionProps, reactantsPopulations, 
		mData, reactantNames, reactantsPopulation, reactionConstants);	

	/*{
		CSCMatrix(reactantsPopulations.size(),
			mData).print();
	}*/

	cout<<"done"<<endl;
	cout<<"Initializing data structure... ";

	CHECK_CUDA_ERRORS();

	Model* model = new Model(reactionProps.size(), reactantsPopulations.size());

	CHECK_CUDA_ERRORS();

	CSCMatrix *reactionsDescCSC=new CSCMatrix(reactantsPopulations.size(),	mData);

	assert(reactionProps.size()==mData.colPtr.size()-1);

	model->m_ConstData.m_reactionsConstants=reactionConstants;

	model->m_StateData.m_reactantsPopulation=reactantsPopulation;

	model->TransformMatrix(reactionsDescCSC);

	delete reactionsDescCSC;

	cout<<"done"<<endl;
	return model;
}
