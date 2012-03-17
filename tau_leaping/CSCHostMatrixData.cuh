#pragma once


#include <thrust/host_vector.h>

//this structure is for supporting CSCMatrix on the host side
struct CSCHostMatrixData
{
	thrust::host_vector<char> rates;
	thrust::host_vector<int> colPtr; 
	thrust::host_vector<int> rowInd;
};
