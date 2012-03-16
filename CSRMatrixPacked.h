#ifndef _CSR_MATRIX_PACKED_H_
#define _CSR_MATRIX_PACKED_H_

//CSR packed matrix. Values and column inidices are packed
//together into one integer array.
#include "defines.h"

class CSRMatrixPacked
{
	vec_int_t m_values_colInd;//packed values/column indices
	vec_int_t m_rowPtr;//where rows begin
	const int m_cols, m_rows;//matrix dimension

public:
	CSRMatrixPacked(int rows, int cols,
		vec_float_t const &values,
		vec_int_t const &rowPtr,
		vec_int_t const &colInd);

	CSRMatrixPacked(CSRMatrixPacked const &other);

	CSRMatrixPacked(int rows, int cols,
		thrust::host_vector<int> const &values_colInd,
		thrust::host_vector<int> const &rowPtr);

	//print to a console
	void print()const;

	//internal checking. Is not used
	bool check()const;

	//Number of rows
	int Rows()const;

	//multiply by vector
	void mv(float *res, float const *x, bool const *isCritical)const;

	//compute a time step for tau-leaping
	float CompTau(bool const *crit,
		int const *x, char const *type, float const *a,
		float eps, float a0)const;

	//fire reactions
	void FireReactions(vec_int_t const &k, vec_int_t &rc)const;

	//Length of a combined values/column indices array
	__host__
	int ValColIndLen()const{return m_values_colInd.size();}

	//Length of a row array
	__host__
	int RowPtrLen()const{return m_rowPtr.size();}

	//Raw device pointer to a combined values/column indices array
	__host__
	int const *DeviceValColIndPtr()const;//{return thrust::raw_pointer_cast<int const>(m_values_colInd.data());}

	//Raw device pointer to a row array
	__host__
	int const *DeviceRowPtr()const;//{return thrust::raw_pointer_cast<int const>(m_rowPtr.data());}
private:
	//supress copying
	CSRMatrixPacked operator =(CSRMatrixPacked const &)const;
};
#endif
