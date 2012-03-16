#ifndef _CSC_MATRIX_H_
#define _CSC_MATRIX_H_

//Compressed Sparse Column (CSC) matrix structure.
//These data are used in the initial stage of the algorithm
//to import data from an XML input file. Then it is transformed 
//into CSC->CSR->CSRPacked

#include "defines.h"
extern cusparseHandle_t g_handle;

//forward declaration
class CSRMatrix;
struct CSCHostMatrixData;

class CSCMatrix
{
	friend class CSRMatrix;
private:
	vec_float_t m_values;//values
	vec_int_t m_colPtr;//begins of columns
	vec_int_t m_rowInd;//row indices
	const int m_rows, m_cols;//matrix' dimension
	cusparseMatDescr_t m_descr;//auxilar structure for cusparse lib

private:
	//initialize cusparse library if it was not initialized before
	void CheckCusparseForDeviceVersion();

protected:

	CSCMatrix(int rows,
		thrust::device_vector<float> const &values,
		thrust::device_vector<int> const &colPtr,
		thrust::device_vector<int> const &rowInd);
public:
	//initialize from host data
	CSCMatrix(int rows,
		CSCHostMatrixData const &data);

	//make it out of CSR matrix
	static
	std::auto_ptr<CSCMatrix> FromCSR(const CSRMatrix &csr);

	//multiply by vector
	void mv_manual(float *res, float const *x);

	//print to a console
	void print()const;

	//number of non-zero values in the matrix
	int size()const;

	//multiple by vector (atomicAdding involved)
	void mv(float *res, float const *x)const;

	//number of Rows
	int Rows()const;

	//fill out arrays with reaction orders 
	//and reactants' (not products') indices
	//TODO: move out of here
	void DetectReactionOrdersAndFillReactantIndices(vec_char_t &reactionOrder,
		vec_int_t &s0, vec_int_t &s)const;

	//create and fill with random values
	static
	std::auto_ptr<CSCMatrix> FillRand(int rows, int cols);

	//convert to CSR matrix
	std::auto_ptr<CSRMatrix> ToCSR()const;
};
#endif
