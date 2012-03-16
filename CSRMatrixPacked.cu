#include "CSRMatrixPacked.h"

#include <thrust/logical.h>
#include <thrust/functional.h>
#include <thrust/find.h>
#include <thrust/sequence.h>
#include <thrust/device_ptr.h>
#include <thrust/device_vector.h>
#include "functions.h"

using thrust::device_ptr;
using thrust::device_vector;
using thrust::minimum;
using thrust::maximum;
using thrust::raw_pointer_cast;
using std::cerr;
using std::cout;
using std::endl;



CSRMatrixPacked::CSRMatrixPacked(int rows, int cols,
		vec_float_t const &values,
		vec_int_t const &rowPtr,
		vec_int_t const &colInd):m_rows(rows), m_cols(cols), m_rowPtr(rowPtr)
{
	m_values_colInd.resize(values.size());

	thrust::transform(values.begin(), values.end(), colInd.begin(),
		m_values_colInd.begin(), ValueAndInd2Int());
}

CSRMatrixPacked::CSRMatrixPacked(CSRMatrixPacked const &other):
	m_cols(other.m_cols), m_rows(other.m_rows),
	m_values_colInd(other.m_values_colInd),
	m_rowPtr(other.m_rowPtr)
{
}


class RowMultPacked
{
	int const * const m_values_colInd;
	float const * const m_x;
	int const * const m_rowPtr;
	bool const *m_isCritical;
public:
	__host__ __device__
	RowMultPacked(int const * values_colInd,
			int const * rowPtr,
			float const * x,
			bool const *isCritical):m_values_colInd(values_colInd), m_rowPtr(rowPtr), m_x(x), m_isCritical(isCritical){}
	inline __host__ __device__
	float operator()(int i)const
	{
		int s=m_rowPtr[i];
		int e=m_rowPtr[i+1]; 
		float res=0.f;
		for(int j=s; j<e; ++j)
		{
			char v;
			int vci;
			split(m_values_colInd[j], v, vci);
			if(m_isCritical[vci])
				continue;
			if(v>0)
				continue;

			res+=(float)v*m_x[vci];
		}
		return res;
	}
};

class RowMultPackedSquared
{
	int const * const m_values_colInd;
	float const * const m_x;
	int const * const m_rowPtr;
	bool const *m_isCritical;
public:
	__host__ __device__
	RowMultPackedSquared(int const * values_colInd,
			int const * rowPtr,
			float const * x,
			bool const *isCritical):m_values_colInd(values_colInd), m_rowPtr(rowPtr), m_x(x), m_isCritical(isCritical){}
	inline __host__ __device__
	float operator()(int i)const
	{
		int s=m_rowPtr[i];
		int e=m_rowPtr[i+1];
		float res=0.f;
		for(int j=s; j<e; ++j)
		{
			char v;
			int vci;
			split(m_values_colInd[j], v, vci);
			if(m_isCritical[vci])
				continue;

			if(v>0)
				continue;

			res+=(float)v*(float)v*m_x[vci];
		}
		return res;
	}
};

__global__ void mv_kernel(int nRows, int const * values_colInd,
			int const * rowPtr,
			float const * x,
			float *resArr)
{
	int i=gridDim.x*blockDim.x*blockIdx.y+blockIdx.x*blockDim.x + threadIdx.x;

	if(i>nRows)
		return;

	int s=rowPtr[i];
	int e=rowPtr[i+1];
	float res=0.f;
	for(int j=s; j<e; ++j)
	{
		char v;
		int vci;
		split(values_colInd[j], v, vci);

		res+=(float)v*x[vci];
	}

	resArr[i]=res;

}

void CSRMatrixPacked::mv(float *res, float const *x, bool const *isCritical)const
{
	thrust::counting_iterator<int> cit(0);

    RowMultPacked fn(raw_pointer_cast(&m_values_colInd[0]),
    			raw_pointer_cast(&m_rowPtr[0]),
    			x,
				isCritical);

	transform(cit, cit+m_rows, device_ptr<float>(res), fn);

	/*const int thrCnt=256;
	dim3 bd(1024, (m_rows/(thrCnt*1024))+1,1);
	//cout<<"Kernel configuration is: "<<bd.x<<"x"<<bd.y<<"x"<<bd.z<<"x"<<thrCnt<<" for "<<m_rows<<" rows\n";

	//int bd=m_rows/thrCnt+1;
	mv_kernel<<<bd,thrCnt>>>(m_rows, raw_pointer_cast(&m_values_colInd[0]),
    			raw_pointer_cast(&m_rowPtr[0]),
    			x,
				res);

	cudaError_t err=cudaGetLastError();
	if(err!=cudaSuccess)
	{
		cerr<<"can not execute mv_kernel. The error is: "<<err<<endl;
		exit(EXIT_FAILURE);
	}*/

}


void CSRMatrixPacked::print()const
{
	/*thrust::host_vector<int> h_v=m_values_colInd;

	std::ostream_iterator<int, char,
             std::char_traits<char> > cout(std::cout, " ");
	std::ostream_iterator<int, char,
             std::char_traits<char> > iout(std::cout, " ");



	std::cout<<std::endl<<"************************";
	std::cout<<std::endl<<"values:"<<std::endl;
	std::copy(h_v.begin(), h_v.end(), cout);

	thrust::host_vector<int> h_i=m_rowPtr;

	std::cout<<std::endl<<"rowPtr:"<<std::endl;
	std::copy(h_i.begin(), h_i.end(), iout);

	h_i=m_colInd;

	std::cout<<std::endl<<"colInd:"<<std::endl;
	std::copy(h_i.begin(), h_i.end(), iout);
	std::cout<<std::endl;*/
}


class RowLenNZ:public thrust::unary_function<int, bool>
{
	int const * m_rowPtr;
public:
	RowLenNZ(int const * rowPtr):
		m_rowPtr(rowPtr)
	{}

	__host__ __device__
	bool operator()(int const &i)const
	{
		int s=m_rowPtr[i];
		int e=m_rowPtr[i+1];
		return s!=e;
	}

};

class RowAbsSumNZ:public thrust::unary_function<int, int>
{
	int const * m_values_colInd;
	int const * m_rowPtr;
public:
	RowAbsSumNZ(int const *values_colInd, int const * rowPtr):
		m_rowPtr(rowPtr), m_values_colInd(values_colInd)
	{}

	__host__ __device__
	bool operator()(int const &i)const
	{
		int s=m_rowPtr[i];
		int e=m_rowPtr[i+1];
		float res=0.f;
		for(int j=s; j<e; ++j)
		{
			char v;
			int vci;
			split(m_values_colInd[j], v, vci);

			res+=fabs((float)v);

		}
		if(res==0.0f)
				return false;
		return true;//res!=0.f;
	}

};

bool CSRMatrixPacked::check()const
{
	thrust::counting_iterator<int> cit(0);
	/*thrust::device_vector<int> aux(m_rows);
	thrust::sequence(aux.begin(), aux.end());*/

	RowLenNZ lenFN(raw_pointer_cast(&m_rowPtr[0]));
	bool res1=thrust::all_of(cit, cit+m_rows, lenFN);
	if(!res1)
	{
		cerr<<"There are empty rows in the matrix\n";
		return false;
	}

	RowAbsSumNZ fn(raw_pointer_cast(&m_values_colInd[0]),
    			raw_pointer_cast(&m_rowPtr[0]));
	return thrust::all_of(cit, cit+m_rows, fn);
}

int CSRMatrixPacked::Rows()const
{
	return m_rows;
}

class RawLenFn
{
	int const * m_rowPtr;
public:
	RawLenFn(int const * rowPtr):m_rowPtr(rowPtr)
	{}

	__device__ inline
	int operator()(int i)
	{
		return m_rowPtr[i+1]-m_rowPtr[i];
	}

};

/*int CSRMatrixPacked::MaxRawLen()const
{
	RawLenFn rlFn(raw_pointer_cast(&m_rowPtr[0]));

	thrust::counting_iterator<int> cit(0);

	return transform_reduce(cit, cit+m_rows,
		rlFn,
		0,  maximum<int>());
}*/

inline __device__ __host__
float GetG(int i, char type, int x)
{
	switch(type)
	{
	case 1:
		return 1.f;
	case 2:
		return 2.f;
	case 3:
		return 2.f+1.f/(x-1.f);
	default:
		//my_assert();
		return 0.f;
	}
}

class TauFn
{
	int const * m_values_colInd;
	int const * m_rowPtr;
	bool const *m_crit;
	int const *m_x;
	float const *m_a;
	char const *m_type;
	float m_eps;
	float const m_a0;
public:
	TauFn(int const * values_colInd,
		  int const * rowPtr,
		  bool const *crit, int const *x, char const *type, float const *a, float eps, float a0):
	  m_crit(crit), m_x(x), m_type(type), m_a(a), m_eps(eps),
	  m_values_colInd(values_colInd), m_rowPtr(rowPtr),
	m_a0(a0)
	{}
	inline __device__ __host__
	float operator()(int i)const
	{
		float mu=0.f;
		float sigma2=0.f;

		int s=m_rowPtr[i];
		int e=m_rowPtr[i+1];
		for(int j=s; j<e; ++j)
		{
			char v;
			int vci;
			split(m_values_colInd[j], v, vci);
			/*if(v>0)
				continue;*/

			float tmp=(float)v*m_a[vci]*(!m_crit[vci]);

			sigma2+=(float)v*tmp;
			mu+=tmp;
		}
		
		//if(m_crit[i])
		//	return 100000000.0f;
		float tmp=max(m_eps*m_x[i]/GetG(i, m_type[i], m_x[i]), 1.f);
		return min(tmp/fabs(mu), tmp*tmp/fabs(sigma2));
	}
private:
	
};

/*__global__ void TauKernel(int sz,
		  int const * values_colInd,
		  int const * rowPtr,
		  bool const *crit, int const *x, char const *type, float const *a, float eps, float a0, 
		  float* aux, TauFn fn)
{
	//const RowMultPackedSquared RMPS(values_colInd, rowPtr, a, crit);
	//const RowMultPacked RMP(values_colInd, rowPtr, a, crit);

	int tIdx=threadIdx.x;
	int i=gridDim.x*blockIdx.y + blockIdx.x;

	if(i>=sz)
		return; 

	
	float mu=0.f;;
	float sigma2=0.f;
	int s=rowPtr[i];
	int e=rowPtr[i+1]; 
	__shared__ float s_mu[32];
	__shared__ float s_sigma2[32];
	s_mu[tIdx]=0.f;
	s_sigma2[tIdx]=0.f;

	int j=tIdx+s;
	if(j<e)
	//for(int j=s; j<e; ++j)
	{
		
		char v;
		int vci;
		split(values_colInd[j], v, vci);
		if(!crit[vci]&&v<0)
		{
			float t=(float)v*a[vci];

			s_mu[tIdx]=t;
			s_sigma2[tIdx]=t*v;
		}
	}

		
	for (unsigned int s=blockDim.x/2; s>0; s>>=1) {
		if (tIdx < s) {
			s_mu[tIdx] += s_mu[tIdx + s];
			s_sigma2[tIdx] += s_sigma2[tIdx + s];
		}
	}

	if(tIdx!=0)
		return;

	mu=s_mu[0];
	sigma2=s_sigma2[0];
	
	float tmp=max(eps*x[i]/GetG(i, type[i], x[i]), 1.f);
	aux[i]=min(tmp/fabs(mu), tmp*tmp/fabs(sigma2));
	
}*/



float CSRMatrixPacked::CompTau(bool const *crit,
							   int const *x, char const *type, float const *a,
							   float eps, float a0)const
{
	thrust::counting_iterator<int> cit(0);

	TauFn fn(raw_pointer_cast(&m_values_colInd[0]), raw_pointer_cast(&m_rowPtr[0]),
		crit, x, type, a, eps, a0);

	/*const int THREADS_PER_BLOCK=32;//warp size
	//const int BLOCKS=m_rows%THREADS_PER_BLOCK?m_rows/THREADS_PER_BLOCK+1:m_rows/THREADS_PER_BLOCK;
	const dim3 BLOCKS(10000, m_rows/10000+1);
	TauKernel<<<BLOCKS,THREADS_PER_BLOCK>>>(m_rows, 
		raw_pointer_cast(&m_values_colInd[0]), raw_pointer_cast(&m_rowPtr[0]),
		crit, x, type, a, eps, a0,
		raw_pointer_cast(aux), fn);*/

	
	float res=transform_reduce(cit, cit+m_rows,
		fn,
		std::numeric_limits<float>::max(),  minimum<float>());
	//thrust::transform(cit, cit+m_rows, aux, fn);
//	float res=thrust::reduce(aux, aux+m_rows, std::numeric_limits<float>::max(),  minimum<float>());
	//float res=*thrust::min_element(aux, aux+m_rows);
	static int i=0;
	if(i==0)
		cout<<"Time Step for the first step is: "<<res<<endl;/*<<
			"configuration is "<<BLOCKS<<"x"<<THREADS_PER_BLOCK<<"="<<BLOCKS*THREADS_PER_BLOCK<<endl<<
			"max raw len is: "<<max_raw_len;*/
	++i;


	/*printf("mu: \n");
	print_vector<float, float>(mu);
	printf("\nsigma2: \n");
	print_vector<float, float>(sigma2);
	printf("\n");*/

	return res;
}

class FireReactionsFn
{
	int const * m_values_colInd;
	int const * m_rowPtr;
	int const *m_k;
	int *m_rc;
public:
	FireReactionsFn(vec_int_t const &values_colInd, vec_int_t const &rowPtr, vec_int_t const &k, vec_int_t &rc):
	  m_values_colInd(raw_pointer_cast(&values_colInd[0])), m_rowPtr(raw_pointer_cast(&rowPtr[0])),
	  m_k(raw_pointer_cast(&k[0])), m_rc(raw_pointer_cast(&rc[0])){}

	  inline __device__
	void operator()(int i)const
	{
		int s=m_rowPtr[i];
		int e=m_rowPtr[i+1];
		for(int j=s; j<e; ++j)
		{
			char v;
			int vci;
			split(m_values_colInd[j], v, vci);
			int k=m_k[vci];

			m_rc[i]+=k*v;

			//atomicAdd(&m_rc[i], k*v);//TODO: what's atomic doing here???
		}
	}
};

void CSRMatrixPacked::FireReactions(vec_int_t const &k, vec_int_t &rc)const
{

	/*printf("number of reactants before: \n");
		print_vector<int, int> (rc);
		printf("\n\n");*/

	thrust::counting_iterator<int> cit(0);
	thrust::for_each(cit, cit+k.size(), FireReactionsFn(m_values_colInd, m_rowPtr, k, rc));
}

//Raw device pointer to a combined values/column indices array
int const *CSRMatrixPacked::DeviceValColIndPtr()const{

	return thrust::raw_pointer_cast(m_values_colInd.data());
}

//Raw device pointer to a row array
int const *CSRMatrixPacked::DeviceRowPtr()const{
	return thrust::raw_pointer_cast(m_rowPtr.data());
}
