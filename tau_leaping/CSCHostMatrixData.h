#pragma once
#include "CSCHostMatrixData.cuh"

//we need to reformat data obtained from an XML input file
//into a Compressed Sparse Column (CSC) format

#include <vector>
#include <string>

#ifdef _MSC_VER
#pragma warning (push)
#pragma warning (disable: 4100 4512 4127)
#endif

#include <thrust/host_vector.h>

#ifdef _MSC_VER
#pragma warning (pop)
#endif


#include <unordered_map>

#ifdef _WIN32
namespace std{
	using namespace tr1;
}
#endif


struct ReactionProp;

void ReadXMLToMatrixStruct(std::vector<ReactionProp> const &reactionProps, 
						   std::unordered_map<std::string, std::pair<int, int> > const &reactantsPopulations, 
						   CSCHostMatrixData &mData,
						   std::vector<std::string> &reactantNames,
						   thrust::host_vector<int> &reactantsPopulation, 
						   thrust::host_vector<float> &reactionConstants);

std::ostream &operator<<(std::ostream &os, CSCHostMatrixData  const &data);