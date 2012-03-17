#include "CSCMatrix.h"
#include "CSRMatrix.h"
#include "CSCHostMatrixData.cuh"
#include <thrust/for_each.h>
#include <thrust/sort.h>
#include <cassert>

using thrust::host_vector;
using thrust::device_vector;
using thrust::raw_pointer_cast;

using std::cout;
using std::endl;

class ColMult
{
	float const * const m_values;
	float const * const m_x;
	float * const m_res;
	int const * const m_colPtr;
	int const * const m_rowInd;
public:
	ColMult(float const * values,
			int const * colPtr,
			int const * rowInd,
			float * res,
			float const * x):m_values(values), m_colPtr(colPtr), m_rowInd(rowInd), m_x(x), m_res(res){}
	inline __device__
	void operator()(int j)const
	{
		int s=m_colPtr[j];
		int e=m_colPtr[j+1];
		for(int i=s; i<e; ++i)
		{
			int ri=m_rowInd[i];
			float toAdd=m_values[i]*m_x[j];
			atomicAdd(&m_res[ri], toAdd);
		}
	}
};

void CSCMatrix::mv(float *res, float const *x)const
{
	thrust::counting_iterator<int> cit(0);
	thrust::fill(thrust::device_ptr<float>(res), thrust::device_ptr<float>(res)+m_rows, 0.0f);
	thrust::for_each(cit, cit+m_cols,
		ColMult(raw_pointer_cast(m_values.data()),
			raw_pointer_cast(m_colPtr.data()),
			raw_pointer_cast(m_rowInd.data()),
			res, x));
}

void CSCMatrix::print()const
{
	thrust::host_vector<float> h_v=m_values;

	std::ostream_iterator<int, char,
             std::char_traits<char> > iout(std::cout, " ");
	std::ostream_iterator<float, char,
             std::char_traits<char> > fout(std::cout, " ");

	std::cout<<std::endl<<"************************";
	std::cout<<std::endl<<"values:"<<std::endl;
	std::copy(h_v.begin(), h_v.end(), fout);

	thrust::host_vector<int> h_i=m_colPtr;

	std::cout<<std::endl<<"colPtr:"<<std::endl;
	std::copy(h_i.begin(), h_i.end(), iout);

	h_i=m_rowInd;

	std::cout<<std::endl<<"rowInd:"<<std::endl;
	std::copy(h_i.begin(), h_i.end(), iout);
	std::cout<<std::endl;
}


std::auto_ptr<CSCMatrix> CSCMatrix::FillRand(int rows, int cols)
{

	float vals[]={-2.f, -1.f, 1.f, 2.f};

	//const int nNz=m_values.size();

	//assert(nNz==m_nNz);
	//assert(m_n+1==mp1);

	printf("Filling up the random matrix %dx%d...",rows, cols);

//	const int mp1=m+1;

	host_vector<int> usedInices(rows);
	thrust::sequence(usedInices.begin(), usedInices.end());
	std::random_shuffle(usedInices.begin(), usedInices.end());

	host_vector<int> cscRowInd, cscColPtr(cols+1);
	host_vector<float> cscValues;
//	int colPtr=0;
	cscColPtr[0]=0;
	for(int j=0; j<cols; ++j)
	{
		const int rN=static_cast<int>(rand()/static_cast<float>(RAND_MAX+1)*4+1);
		assert(rN>=1&&rN<=4);
		cscColPtr[j+1]=cscColPtr[j]+rN;
		for(int k=0; k<rN; ++k)
		{
			int rInd;
			if(usedInices.empty())
				rInd=static_cast<int>(rand()/static_cast<float>(RAND_MAX+1)*rows);
			else
			{
				rInd=usedInices.back();
				usedInices.pop_back();
			}
			assert(rInd>=0&&rInd<rows);
			const int vInd=static_cast<int>(rand()/static_cast<float>(RAND_MAX+1)*4);
			assert(vInd>=0&&vInd<4);
			cscRowInd.push_back(rInd);
			cscValues.push_back(vals[vInd]);
		}
		thrust::sort(&cscRowInd[cscColPtr[j]],&cscRowInd[cscColPtr[j+1]]);
	}

	if(!usedInices.empty())
	{
		std::cerr<<"There are not enough reactions; some reactants are unesed\n";
		exit(EXIT_FAILURE);
	}

	device_vector<int> d_cscRowInd(cscRowInd), d_cscColPtr(cscColPtr);
	device_vector<float> d_cscValues(cscValues);

	if(g_handle==0)
		cusparseCreate(&g_handle);

	printf(" done\n");

	return std::auto_ptr<CSCMatrix> (new CSCMatrix(rows, d_cscValues, d_cscColPtr, d_cscRowInd));
}

int CSCMatrix::size()const
{
	return m_values.size();
}

std::auto_ptr<CSRMatrix> CSCMatrix::ToCSR()const
{

	device_vector<int> d_csrRowPtr(m_rows+1), d_csrColInd(size());
	device_vector<float> d_csrValues(size());

	cusparseStatus_t status;
	if(g_handle==0)
	{
		status=cusparseCreate(&g_handle);
		if (status != CUSPARSE_STATUS_SUCCESS)
		{
			fprintf(stderr, "cusparseCreate initialization failed\n");
			fprintf(stderr, "Status code is: %d\n", status);
			exit(EXIT_FAILURE);
		}
	}

	//columns and rows numbers are swapped since we perform the back trancformation
	status=cusparseScsr2csc(g_handle, m_cols, m_rows,
		raw_pointer_cast(&m_values[0]), raw_pointer_cast(&m_colPtr[0]), raw_pointer_cast(&m_rowInd[0]),
		raw_pointer_cast(&d_csrValues[0]), raw_pointer_cast(&d_csrColInd[0]), raw_pointer_cast(&d_csrRowPtr[0]),
		1, CUSPARSE_INDEX_BASE_ZERO);

	if (status != CUSPARSE_STATUS_SUCCESS)
	{
		fprintf(stderr, "Matrix conversion failed\n");
		fprintf(stderr, "Status code is: %d\n", status);
		exit(EXIT_FAILURE);
	}
	return std::auto_ptr<CSRMatrix> (new CSRMatrix(m_rows, m_cols, d_csrValues, d_csrRowPtr, d_csrColInd));
}

int CSCMatrix::Rows()const
{
	return m_rows;
}

class OrderDetectorAndIndFiller
{
	int const * m_colPtr;
	int const * m_rowInd;
	float const *m_values;
public:
	OrderDetectorAndIndFiller(vec_float_t const &values, vec_int_t const &colPtr, vec_int_t const &rowInd):
#ifdef DEVICE
	  m_values(raw_pointer_cast(values.data())),
	  m_colPtr(raw_pointer_cast(colPtr.data())),
	  m_rowInd(raw_pointer_cast(rowInd.data())){}
#elif defined HOST
#error "host version is not implemented"
#endif

	__host__ __device__
	inline
	thrust::tuple<char, int, int> operator()(int i)const
	{
#ifdef DEVICE
		int startInd=m_colPtr[i];
		int endInd=m_colPtr[i+1];
		char firstOrder=0, secondOrder=0;
		int ind0=-1, ind1=-1;
		for(int i=startInd; i<endInd; ++i)
		{
			if(m_values[i]>=0)
				continue;

			if(firstOrder==0)
			{
				firstOrder=static_cast<char>(m_values[i]);
				ind0=m_rowInd[i];
			}
			else
			{
				secondOrder=static_cast<char>(m_values[i]);
				ind1=m_rowInd[i];
			}
		}
		if(firstOrder==-2)
			return thrust::make_tuple(3, ind0, ind1);
		if(firstOrder==-1&&secondOrder==-1)
			return thrust::make_tuple(2, ind0, ind1);
		if(firstOrder==-1)
			return thrust::make_tuple(1, ind0, ind1);
		return thrust::make_tuple(-1, -1, -1);

#elif defined HOST
#error "host version is not implemented"
#endif

	}
};

void CSCMatrix::DetectReactionOrdersAndFillReactantIndices(vec_char_t &reactionOrder,
														   vec_int_t &s0, vec_int_t &s1)const
{
	s0.resize(m_cols);
	s1.resize(m_cols);
	reactionOrder.resize(m_cols, -1);

	thrust::counting_iterator<int> cit(0);
	thrust::transform(cit, cit+m_cols,
		make_zip_iterator(thrust::make_tuple(reactionOrder.begin(), s0.begin(), s1.begin())),
		OrderDetectorAndIndFiller(m_values, m_colPtr, m_rowInd));
}

CSCMatrix::CSCMatrix(int rows,
		CSCHostMatrixData const &data):
		m_values(data.rates), m_rows(rows), m_cols(data.colPtr.size()-1), m_rowInd(data.rowInd), m_colPtr(data.colPtr)
{
    CheckCusparseForDeviceVersion();
}

CSCMatrix::CSCMatrix(int rows,
		thrust::device_vector<float> const &values,
		thrust::device_vector<int> const &colPtr,
		thrust::device_vector<int> const &rowInd):
m_rows(rows), m_cols(colPtr.size()), m_colPtr(colPtr), m_rowInd(rowInd), m_values(values)
{
    CheckCusparseForDeviceVersion();
}

void CSCMatrix::CheckCusparseForDeviceVersion()
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
