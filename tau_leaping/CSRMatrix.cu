#include "CSRMatrix.h"
#include "CSCMatrix.h"
#include "CSRMatrixPacked.h"
//#include "CSCMatrixPacked.h"


#include <thrust/host_vector.h>
#include <thrust/iterator/constant_iterator.h>
#include <thrust/iterator/counting_iterator.h>
#include <thrust/sort.h>
#include <cassert>
#include <set>
#include <cmath>
#include <thrust/sequence.h>
#include <thrust/scan.h>

#include <functional>

using namespace thrust;
using std::cout;
using std::endl;

//cusparseHandle_t g_handle=NULL;
/*
CSRMatrix::CSRMatrix(int nNz, int nRows, int nCols):m_m(nRows), m_n(nCols), m_descr(NULL)
{
	CheckCusparseForDeviceVersion();

	m_values.resize(nNz);
	m_colInd.resize(nNz);
	m_rowPtr.resize(nCols+1);
}
*/
CSRMatrix::CSRMatrix(int rows, int cols,
		thrust::host_vector<float> const &values,
		thrust::host_vector<int> const &rowPtr,
		thrust::host_vector<int> const &colInd):
m_rows(rows), m_cols(cols), m_values(values), m_rowPtr(rowPtr), m_colInd(colInd)
{
	CheckCusparseForDeviceVersion();
}

void CSRMatrix::CheckCusparseForDeviceVersion()
{
#ifdef DEVICE
	if(g_handle==0)
		cusparseCreate(&g_handle);

	cusparseStatus_t status= cusparseCreateMatDescr(&m_descr);
	if (status != CUSPARSE_STATUS_SUCCESS)
	{
		fprintf(stderr, "Matrix descriptor initialization failed\n");
		exit(EXIT_FAILURE);
	}
	status=cusparseSetMatType(m_descr, CUSPARSE_MATRIX_TYPE_GENERAL);
	assert(status== CUSPARSE_STATUS_SUCCESS);
	status=cusparseSetMatIndexBase(m_descr, CUSPARSE_INDEX_BASE_ZERO);
	assert(status== CUSPARSE_STATUS_SUCCESS);
	status=cusparseSetMatDiagType(m_descr, CUSPARSE_DIAG_TYPE_NON_UNIT);
	assert(status== CUSPARSE_STATUS_SUCCESS);
#endif
}

/*
CSRMatrix::CSRMatrix(int nRows, int nCols):m_m(nRows), m_n(nCols), m_descr(NULL)
{
	if(g_handle==0)
		cusparseCreate(&g_handle);

	cusparseStatus_t status= cusparseCreateMatDescr(&m_descr);
	if (status != CUSPARSE_STATUS_SUCCESS)
	{
		fprintf(stderr, "Matrix descriptor initialization failed\n");
		exit(EXIT_FAILURE);
	}
	status=cusparseSetMatType(m_descr, CUSPARSE_MATRIX_TYPE_GENERAL);
	assert(status== CUSPARSE_STATUS_SUCCESS);
	status=cusparseSetMatIndexBase(m_descr, CUSPARSE_INDEX_BASE_ZERO);
	assert(status== CUSPARSE_STATUS_SUCCESS);
	status=cusparseSetMatDiagType(m_descr, CUSPARSE_DIAG_TYPE_NON_UNIT);
	assert(status== CUSPARSE_STATUS_SUCCESS);


	m_rowPtr.resize(nCols+1);
}*/
/*
std::auto_ptr<CSRMatrix> CSRMatrix::FillRand(int m, int n, float avgElsPerColumn)
{*/

	/*float vals[]={-2.f, -1.f, 0.f, 1.f, 2.f};

	//const int nNz=m_values.size();

	//assert(nNz==m_nNz);
	//assert(m_n+1==mp1);

	printf("Filling up the random matrix %dx%d\n",m, n);

	const int mp1=m+1;

	host_vector<int> cscRowInd, cscColPtr(m+1);
	host_vector<float> cscValues;
	int colPtr=0;
	cscColPtr[0]=0;
	for(int j=0; j<m; ++j)
	{
		const int rN=static_cast<int>(rand()/static_cast<float>(RAND_MAX)*4+1);
		assert(rN>=1&&rN<=4);
		cscColPtr[j+1]=cscColPtr[j]+rN;
		for(int k=0; k<rN; ++k)
		{
			const int rInd=static_cast<int>(rand()/static_cast<float>(RAND_MAX)*n);
			assert(rInd>=0&&rInd<n);
			const int vInd=static_cast<int>(rand()/static_cast<float>(RAND_MAX)*5);
			assert(vInd>=0&&vInd<5);
			cscRowInd.push_back(rInd);
			cscValues.push_back(vals[vInd]);
		}
		thrust::sort(&cscRowInd[cscColPtr[j]],&cscRowInd[cscColPtr[j+1]]);
	}

	device_vector<int> d_cscRowInd(cscRowInd), d_cscColPtr(cscColPtr);
	device_vector<float> d_cscValues(cscValues);

	device_vector<int> d_csrRowPtr(n+1), d_csrColInd(cscRowInd.size());
	device_vector<float> d_csrValues(cscValues.size());

	if(g_handle==0)
		cusparseCreate(&g_handle);*/

	//return std::auto_ptr<CSCMatrix> (new CSCMatrix(m, n, d_cscValues, d_cscColPtr, d_csrRowInd));

	/*cusparseStatus_t status=cusparseScsr2csc(g_handle, m, n,
		raw_pointer_cast(&d_cscValues[0]), raw_pointer_cast(&d_cscRowInd[0]), raw_pointer_cast(&d_cscColPtr[0]),
		raw_pointer_cast(&d_csrValues[0]), raw_pointer_cast(&d_csrRowPtr[0]), raw_pointer_cast(&d_csrColInd[0]),
		1, CUSPARSE_INDEX_BASE_ZERO);

	if (status != CUSPARSE_STATUS_SUCCESS)
	{
		fprintf(stderr, "Matrix conversion failed\n");
		exit(EXIT_FAILURE);
	}*/


    /////////////////////////////
	//generating random indices//
	/////////////////////////////
	/*std::set<int> indices;
	std::vector<int> unfilledRows;
	for(int i=0; i<n; ++i)
	{
		unfilledRows.push_back(i);
	}

	//float const avgNnzInCol=values.size()/static_cast<float>(m_n);
	printf("Avg elems in column: %f\n", avgElsPerColumn);

	std::random_shuffle(unfilledRows.begin(), unfilledRows.end());
	for(int j=0; j<n; ++j)
	{
		const float f=rand()/static_cast<float>(RAND_MAX)*avgElsPerColumn;
		const int rN=static_cast<int>(f*2+0.5);
		//const int rN=round(2*f);

		//printf("%d %d %f\n", j, rN,2*f);
		std::set<int> rowIndeces;
		for(int ri=0; ri<rN; ++ri)
		{
			int rInd;
			if(!unfilledRows.empty())
			{
				rInd=unfilledRows.back();
				unfilledRows.pop_back();
			}else
			{
				for(;;)
				{
					rInd =static_cast<int>(rand()/static_cast<float>(RAND_MAX)*m);
					if(rowIndeces.find(rInd)==rowIndeces.end())
					{
						rowIndeces.insert(rInd);
						break;
					}
				}
			}
			assert(rInd*n+j>=0);
			indices.insert(rInd*n+j);
		//	if(indices.size()==values.size())
			//	break;
			//printf("%d ", rInd);
		}
		//if(indices.size()==values.size())
		//	break;
		//printf("\n");
	}
	printf("\nExact avg nnz per column is: %f\n", indices.size()/static_cast<float>(n));
	if(!unfilledRows.empty())
	{
		fprintf(stderr, "Too few reactions for reactants\n");
		exit(EXIT_FAILURE);
	}

	///////////////////////////////////
	//filling up the matix' structure//
	///////////////////////////////////
	int const nNz=indices.size();
	thrust::host_vector<float> values(nNz);
	thrust::host_vector<int> colInd(nNz);
	thrust::host_vector<int> rowPtr(mp1);


	int i=0;
	int ptrInd=0;
	int prevRow=0;
	rowPtr[i]=ptrInd;
	for(std::set<int>::const_iterator it=indices.begin(); it!=indices.end(); ++it, ++i)
	{
		//printf("%d\n", *it);
		values[i]=(float)(int)((rand()/static_cast<float>(RAND_MAX)-0.5)*4.5+0.5);
		int const ri=(*it)/n;
		if(ri!=prevRow)
		{
			rowPtr[(ptrInd++)+1]=i;
			prevRow=ri;
		}
		int const ci=(*it)%n;
		assert(0<=ci&&ci<n);
		//printf("%d %d %d\n", *it, ri, ci);
		colInd[i]=ci;
	}

	//printf("qqq %d %d\n", m_m, i);
	rowPtr[m]=i;*/

	/////////////////////////////////////////////////
	//create a new matrix based on host's structure//
	/////////////////////////////////////////////////
	//return std::auto_ptr<CSRMatrix> (new CSRMatrix(m, n, d_csrValues, d_csrRowPtr, d_csrColInd));
	/*return std::auto_ptr<CSRMatrix> (NULL);
}*/


//modify values by multiplying them by scalar
void CSRMatrix::mc(float x)
{
	thrust::constant_iterator<float> cit(x);
	thrust::transform(m_values.begin(), m_values.end(), cit, m_values.begin(),
		thrust::multiplies<float>());
}

class Square
{
public:
	inline __host__ __device__
	float operator()(float val)const
	{
		return val*val;
	}
};

void CSRMatrix::square()
{
	thrust::transform(m_values.begin(), m_values.end(), m_values.begin(),
		Square());
}

//multiply by vector
void CSRMatrix::mv(float *res, float const *x)const
{

	/*thrust::host_vector<float> h_prop(m_n);
	thrust::copy(thrust::device_ptr<float const>(x), thrust::device_ptr<float const>(x)+m_n, h_prop.begin());

	std::ostream_iterator<float, char,
             std::char_traits<char> > out(std::cout, " ");

	std::cout<<std::endl;
	std::copy(h_prop.begin(), h_prop.end(), out);*/

	//print();

	using thrust::raw_pointer_cast;
#ifdef DEVICE
	assert(g_handle);
	cusparseStatus_t status=cusparseScsrmv(g_handle, CUSPARSE_OPERATION_NON_TRANSPOSE,
		m_rows, m_cols, 1.0f, m_descr,
		raw_pointer_cast(m_values.data()),
		raw_pointer_cast(m_rowPtr.data()),
		raw_pointer_cast(m_colInd.data()),
		x, 0, res);
	if(status!=CUSPARSE_STATUS_SUCCESS)
	{
		fprintf(stderr, "status= %d\n", status);
		assert(false);
		exit(EXIT_FAILURE);
	}
#elif defined HOST
#error "host version is not implemented"
#endif


	//print();
}

class RowMult
{
	float const * const m_values;
	float const * const m_x;
	int const * const m_rowPtr;
	int const * const m_colInd;
public:
	RowMult(float const * values,
			int const * rowPtr,
			int const * colInd,
			float const * x):m_values(values), m_rowPtr(rowPtr), m_colInd(colInd), m_x(x){}
	inline __host__ __device__
	float operator()(int i)const
	{
		int s=m_rowPtr[i];
		int e=m_rowPtr[i+1];
		float res=0.f;
		for(int j=s; j<e; ++j)
		{
			res+=m_values[j]*m_x[m_colInd[j]];
		}
		return res;
	}
};

void CSRMatrix::mv_manual_row(float *res, float const *x)const
{
	thrust::counting_iterator<int> cit(0);
	//thrust::fill(device_ptr<float>(res), device_ptr<float>(res)+m_rows, 0.f);
	transform(cit, cit+m_rows, device_ptr<float>(res),
			RowMult(raw_pointer_cast(m_values.data()),
			raw_pointer_cast(m_rowPtr.data()),
			raw_pointer_cast(m_colInd.data()),
			x));
}

void CSRMatrix::print()const
{
	thrust::host_vector<float> h_v=m_values;

	std::ostream_iterator<float, char,
             std::char_traits<char> > fout(std::cout, " ");
	std::ostream_iterator<int, char,
             std::char_traits<char> > iout(std::cout, " ");



	std::cout<<std::endl<<"************************";
	std::cout<<std::endl<<"values:"<<std::endl;
	std::copy(h_v.begin(), h_v.end(), fout);

	thrust::host_vector<int> h_i=m_rowPtr;

	std::cout<<std::endl<<"rowPtr:"<<std::endl;
	std::copy(h_i.begin(), h_i.end(), iout);

	h_i=m_colInd;

	std::cout<<std::endl<<"colInd:"<<std::endl;
	std::copy(h_i.begin(), h_i.end(), iout);
	std::cout<<std::endl;
}

int CSRMatrix::size()const
{
	return m_values.size();
}

std::auto_ptr<CSCMatrix> CSRMatrix::ToCSC()const
{
	thrust::device_vector<float> values(size());
	thrust::device_vector<int> colPtr(m_cols+1);
	thrust::device_vector<int> rowInd(size());

	cusparseStatus_t status= cusparseScsr2csc(g_handle, m_rows, m_cols,
		raw_pointer_cast(&m_values[0]), raw_pointer_cast(&m_rowPtr[0]), raw_pointer_cast(&m_colInd[0]),
		raw_pointer_cast(&values[0]),  raw_pointer_cast(&rowInd[0]), raw_pointer_cast(&colPtr[0]),
		1, CUSPARSE_INDEX_BASE_ZERO);

	if (status != CUSPARSE_STATUS_SUCCESS)
	{
		fprintf(stderr, "Matrix conversion failed\n");
		exit(EXIT_FAILURE);
	}

	return std::auto_ptr<CSCMatrix>(new CSCMatrix(m_cols, values, colPtr, rowInd));

}

/*
std::auto_ptr<CSCMatrixPacked> CSRMatrix::ToCSCPacked()const
{
	thrust::device_vector<float> values(size());
	thrust::device_vector<int> colPtr(m_cols+1);
	thrust::device_vector<int> rowInd(size());

	cusparseStatus_t status= cusparseScsr2csc(g_handle, m_rows, m_cols,
		raw_pointer_cast(&m_values[0]), raw_pointer_cast(&m_rowPtr[0]), raw_pointer_cast(&m_colInd[0]),
		raw_pointer_cast(&values[0]), raw_pointer_cast(&rowInd[0]), raw_pointer_cast(&colPtr[0]),
		1, CUSPARSE_INDEX_BASE_ZERO);

	if (status != CUSPARSE_STATUS_SUCCESS)
	{
		fprintf(stderr, "Matrix conversion failed\n");
		exit(EXIT_FAILURE);
	}

	return std::auto_ptr<CSCMatrixPacked>(new CSCMatrixPacked(m_rows, m_cols, values, colPtr, rowInd));
}*/

std::auto_ptr<CSRMatrixPacked> CSRMatrix::ToCSRPacked()const
{
	return std::auto_ptr<CSRMatrixPacked>(new CSRMatrixPacked(m_rows, m_cols, m_values, m_rowPtr, m_colInd));
}

int CSRMatrix::Rows()const
{
	return m_rows;
}

class RowLengthFn
{
	int const *m_rowStartInd;
public:
	RowLengthFn(int const *rowStartInd):m_rowStartInd(rowStartInd)
	{}

	inline __host__ __device__
	int operator()(int i)const
	{
		return m_rowStartInd[i+1]-m_rowStartInd[i];
	}
};

class Comp
{
public:
	inline __host__ __device__
	bool operator()(int lhs, int rhs)const
	{
		return lhs>rhs;
	}
};



void CSRMatrix::sort()
{
	thrust::counting_iterator<int> it(0);
	device_vector<int> rIndices(m_rows);
	device_vector<int> rCapacity(m_rows+1,0);

	thrust::transform(it, it+m_rows, rCapacity.begin()+1,
		RowLengthFn(raw_pointer_cast(&m_rowPtr[0])));
	thrust::sequence(rIndices.begin(), rIndices.end());
	thrust::sort_by_key(rIndices.begin(), rIndices.end(), rCapacity.begin()+1, Comp());

	thrust::inclusive_scan(rCapacity.begin(),rCapacity.end(),rCapacity.begin());

//	thrust::host_vector r1(1), r2(1);
//	thrust::copy(rCapacity.begin()+m_rows,rCapacity.end(), r1.begin());
//	thrust::copy(rCapacity.begin()+m_rows,rCapacity.end(), r1.begin());
	assert(thrust::equal(rCapacity.begin()+m_rows,rCapacity.end(), m_rowPtr.end()-1));
}

std::auto_ptr<CSRMatrix> CSRMatrix::Clone()const
{
	vec_float_t values=m_values;
	vec_int_t rowPtr=m_rowPtr;
	vec_int_t colInd=m_colInd;
	return std::auto_ptr<CSRMatrix>(new CSRMatrix(m_rows, m_cols, values, rowPtr, colInd));
}
