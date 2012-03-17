#ifndef _CSR_MATRIX_H_
#define _CSR_MATRIX_H_


//Compressed Sparse Raw (CSR) matrix structure.
//This data are not used directly int the algorithm, but
//used in a toolchain for making Packed CSR matrix


#include "defines.h"
extern cusparseHandle_t g_handle;
class CSCMatrix;
class CSCMatrixPacked;
class CSRMatrixPacked;

class CSRMatrix
{
	friend class CSCMatrix;
private:
	vec_float_t m_values;//matrix' values
	vec_int_t m_rowPtr;//rows' begin indices
	vec_int_t m_colInd;//column indices
	const int m_cols, m_rows;//matrix' dimension
	cusparseMatDescr_t m_descr;//special structure for cusparce library

private:
	void CheckCusparseForDeviceVersion();

	//disable copying
	CSRMatrix operator=(CSRMatrix const&)const;
	CSRMatrix(CSRMatrix const&);

protected:
	//initialize it from host structures
	CSRMatrix(int rows, int cols,
		thrust::host_vector<float> const &values,
		thrust::host_vector<int> const &rowPtr, 
		thrust::host_vector<int> const &colInd);

	
public:

	//make a copy of a matrix
	std::auto_ptr<CSRMatrix> Clone()const;

	//modify values by multiplying them by scalar
	void mc(float x);

	//modify values by squaring them
	void square();

	//multiply by vector with the help of cusparse lib
	void mv(float *res, float const *x)const;

	//manually multiply by vector
	void mv_manual_row(float *res, float const *x)const;

	//print structure/values to a console
	void print()const;

	//number of Rows
	int Rows()const;

	//how many non-zero values a stored into it
	int size()const;

	//sorting rows in a lexicographical order
	void sort();

	//convert it into a CSC matrix (cusparse lib used)
	std::auto_ptr<CSCMatrix> ToCSC()const;

	//not implemented 
	//std::auto_ptr<CSCMatrixPacked> ToCSCPacked()const;

	//convert it into a CSR packed matrix
	std::auto_ptr<CSRMatrixPacked> ToCSRPacked()const;
};

#endif