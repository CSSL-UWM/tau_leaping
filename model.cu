#include "model.h"
//#include "CSCMatrix.h"
//#include "CSRMatrix.h"
#include "CSRMatrixPacked.h"
#include <thrust/iterator/counting_iterator.h>
#include <thrust/random/normal_distribution.h>
#include <thrust/sort.h>
#include <thrust/transform_reduce.h>
#include <thrust/tuple.h>
#include <thrust/transform_scan.h>
#include <thrust/binary_search.h>
#include <thrust/distance.h>
#include <thrust/logical.h>
#include <thrust/random.h>
#ifdef EXTRA_CHECK
#include <thrust/count.h>
#endif
#include <cassert>
#include "functions.h"
#include "d_rng.cuh"
#include <fstream>
#include "MyTime.h"

#ifdef INTERNAL_TIMING
#include "CudaEventTimer.h"
#endif

using thrust::raw_pointer_cast;
using thrust::host_vector;
using std::vector;
using std::string;
using std::cout;
using std::endl;
using std::ofstream;
using std::cerr;

extern std::string g_logFN;
extern bool g_doLogging;


Model::Model(int channels_M, int reactants_N):
	m_StateData(channels_M, reactants_N),
	m_channelsCount(channels_M),
	m_reactantsCount(reactants_N),
	m_tau(0.f),
	m_ConstData(channels_M)
	/*,
	m_sortedIndices(channels_M)*/
{
	/*thrust::counting_iterator<int> cit(0);
	thrust::copy(cit, cit+m_reactionsCount, m_sortedIndices.begin());*/
	
}

Model::~Model()
{
	//delete m_reactionsDesc;
}

	/*
class ReactionSorter
{
	Model::vec_float_t const &m_r;
public:
	ReactionSorter(Model::vec_float_t const &r):m_r(r){}
	bool operator ()(int lhsI, int rhsI)const
	{
		return m_r[lhsI]<m_r[rhsI];
	}
};*/

class ReactionClassifier:public thrust::unary_function<float, bool>
{
	float const m_treshold;
public:
	ReactionClassifier(float treshold):m_treshold(treshold){}

	__host__ __device__ inline
	bool operator ()(float val)const
	{
		return val<m_treshold;
	}
};

/*
void Model::ClassifyReactions()
{
	//working variant based on sorting. Too complex though
	//thrust::sequence(m_classifyInd.begin(), m_classifyInd.end());
	//thrust::sort(m_classifyInd.begin(), m_classifyInd.end(), ReactionClassifier(m_reactionsPropensity));

	//copy noncritical reaction's indices to the begin of classifying array
	thrust::counting_iterator<int> cit(0);
	vec_int_t::iterator m_firstCritInd=
		thrust::copy_if(cit, cit+ReactionsCount(), m_reactionsPropensity.begin(), m_criticalInd.begin(), 
			thrust::not1(ReactionClassifier(10000)));

	
	//copy critical reaction's indices to the end of classifying array
	thrust::copy_if(cit, cit+ReactionsCount(), m_reactionsPropensity.begin(), m_firstCritInd, 
		ReactionClassifier(10000));


	{
		printf("Noncritical reactions:\n");
		for(vec_int_t::const_iterator it= m_criticalInd.begin(); it!=m_firstCritInd; ++it)
		{
			int ind=*it;
			printf("%f %d\n", m_reactionsPropensity[ind], ind);
		}

		printf("\nCritical reactions:\n");
		for(vec_int_t::const_iterator it= m_firstCritInd; it!=m_criticalInd.end(); ++it)
		{
			int ind=*it;
			printf("%f %d\n", m_reactionsPropensity[ind], ind);
		}
	}
}*/

/*void Model::FillRandom()
{
	thrust::host_vector<float> rand_floats(ReactionsCount());

	thrust::generate(rand_floats.begin(),rand_floats.end(), rand);
	m_reactionsPropensity=rand_floats;

	//if(m_reactionsDesc)
		//delete m_reactionsDesc;

	int nnz=ReactionsCount()*2;
	//m_reactionsDesc=new CSRMatrix(nnz, ReactantsCount(), ReactionsCount());

	//m_reactionsDesc->FillRand();
	m_reactionsDescCSC=CSCMatrix::FillRand(ReactantsCount(), ReactionsCount());
}*/


/*
class MuFunctor
{
	vec_int_t const &m_reactionsInd;//the reaction indices
	CSRMatrix const &m_reactionsDesc;//the reaction descriptons - how many particles participate in a reaction
	vec_float_t const &m_reactionsPropensity;
public:
	MuFunctor(vec_int_t const &reactionsInd,
			  CSRMatrix const &reactionsDesc,
			  vec_float_t const &reactionsPropensity):
	m_reactionsInd(reactionsInd), m_reactionsDesc(reactionsDesc), m_reactionsPropensity(reactionsPropensity){};
	__host__ __device__ inline 
	float operator()(int i)
	{
		return 1.0;
	}
};

void Model::ComputeMu()
{
	thrust::counting_iterator<int> cit(0);
	thrust::transform(cit, cit+ReactantsCount() ,m_Mu.begin(), 
		MuFunctor(m_reactionsInd, *m_reactionsDesc, m_reactionsPropensity));
}

void Model::ComputeSigma2()
{

}*/

typedef thrust::tuple<float, float> two_floats_t;
typedef thrust::tuple<int, int, int, char, float> input_t;

inline __host__ __device__
two_floats_t operator+(two_floats_t const &a1,
									  two_floats_t const &a2)
{
	return thrust::make_tuple(thrust::get<0>(a1)+thrust::get<0>(a2), 
		thrust::get<1>(a1)+thrust::get<1>(a2));
}

class PropFunctor
{
	int const *m_x;
	int const m_nc;
	float *m_a;
	bool *m_c;
public:
	PropFunctor(int nc, int const *x, float *a, bool *c):m_x(x),m_a(a), m_nc(nc), m_c(c){}

	
	inline __device__ __host__
	two_floats_t operator()(input_t const &arg)
	{
		int ind, s0, s1;
		char type;
		float rc;
		thrust::tie(ind, s0, s1, type, rc)=arg; 
		int const x0=m_x[s0], x1=m_x[s1];

		bool isCritical;//=IsCritical(type, m_nc, x0, x1);
		float prop;//=Propensity(rc, type, x0,x1);
		IsCriticalAndPropensity(type, x0, x1, m_nc, rc, &prop, &isCritical);

		m_a[ind]=prop;
		m_c[ind]=isCritical;

		return thrust::make_tuple(prop, isCritical?prop:0.f);
	}

};

float g_scanTime=0.f;

void Model::CalcPropensities()
{

#ifdef INTERNAL_TIMING
	CudaEventTimer timer;
	timer.Start();
#endif

	using thrust::raw_pointer_cast;

	//thrust::counting_iterator<int> cit(0);

	PropFunctor unary_op(5, 
		raw_pointer_cast(&m_StateData.m_reactantsPopulation[0]), 
		raw_pointer_cast(&m_StateData.m_reactionsPropensity[0]),
		raw_pointer_cast(&m_StateData.m_criticalClassifier[0]));
    thrust::plus<thrust::tuple<float, float> > assoc_op;

	thrust::counting_iterator<int> cit(0);

	thrust::transform_inclusive_scan(
		thrust::make_zip_iterator(thrust::make_tuple(
			cit, m_ConstData.m_s0.begin(), m_ConstData.m_s1.begin(), m_ConstData.m_ReactionOrder.begin(), 
			m_ConstData.m_reactionsConstants.begin())), 
		thrust::make_zip_iterator(thrust::make_tuple(
			cit+ChannelsCount(), m_ConstData.m_s0.end(), m_ConstData.m_s1.end(), m_ConstData.m_ReactionOrder.end(), 
			m_ConstData.m_reactionsConstants.end())),
		thrust::make_zip_iterator(thrust::make_tuple(m_StateData.m_rPropPartSum.begin(), m_StateData.m_rPropNonCritPartSum.begin())),
		unary_op,
		assoc_op);
	
#ifdef INTERNAL_TIMING
	g_scanTime+=timer.ElapsedTime();
#endif
}
/*
class MinTauSep
{
	float m_eps;
public:
	MinTauSep(float eps):m_eps(eps){}

	inline __host__ __device__
	float operator()(thrust::tuple<bool, float, float, float, float> const &arg)
	{
		bool crit;
		float mu, sigma2, x, g;
		thrust::tie(crit, mu, sigma2, x, g)=arg;

		if(crit)
			return 1.0f;
		float tmp=max(m_eps*x/g, 1.f);
		return min(tmp/fabs(mu), tmp/fabs(sigma2));
	}
};


float Model::CalcTau2()const
{
	MinTauSep fn(m_eps);
	return thrust::transform_reduce(
		make_zip_iterator(make_tuple(d_crit.begin(), d_mu.begin(), d_sigma2.begin(), d_x.begin(), d_g.begin())),
		make_zip_iterator(make_tuple(d_crit.end(), d_mu.end(), d_sigma2.end(), d_x.end(), d_g.end())),
		fn,
		std::numeric_limits<float>::max(),  thrust::minimum<float>());
}*/

float g_prop=0.f, g_tau=0.f, g_ft=0.f, g_fr=0.f, g_fr1=0.f,  g_fr2=0.f, g_frr=0.f, g_frt=0.f, g_sr=0.f;//, g_t2=0, g_t2=0;
float g_step=0.f;

struct IsNegative
{
	__device__ inline
	bool operator ()(int n)
	{
		return n<0;
	}
};

float Model::TauLeapingTimeStep(float eps, float a0)
{
	return m_ConstData.m_reactionsDesc->CompTau(raw_pointer_cast(&m_StateData.m_criticalClassifier[0]),
			raw_pointer_cast(&m_StateData.m_reactantsPopulation[0]),
			raw_pointer_cast(&m_ConstData.m_ReactionOrder[0]),
			raw_pointer_cast(&m_StateData.m_reactionsPropensity[0]), 
			eps, a0);
}

float Model::Step(float eps)
{
	//compute propensities at first...
	//MyTime::Time_t tb=MyTime::CTime();
//	cout<<"b\n";
#ifdef INTERNAL_TIMING
	CudaEventTimer timer1;
	timer1.Start(); 
#endif
	CalcPropensities();

	
#ifdef INTERNAL_TIMING
	g_prop+=timer1.ElapsedTime();
#endif

	//compute time step for tau-leaping algorithm
	float a0=m_StateData.m_rPropPartSum.back();//TODO: move to a function
	
#ifdef INTERNAL_TIMING
	CudaEventTimer timer2;
	timer2.Start();
#endif
	float const tau1=TauLeapingTimeStep(eps, a0);
	
#ifdef INTERNAL_TIMING
	g_tau+=timer2.ElapsedTime();
#endif
	
	//choose direct or tau-leaping branch and invoke 
	//appropriate algorithm
	int i=0;
	float tau;
	if(tau1<m_N/a0)
	{
		//direct branch
		printf("tau1=%f\n", tau1);
		printf("firing 100 direct steps, partial sums: \n");
		
		//print_vector<float, float> (m_StateData.m_rPropPartSum);
	   // printf("\n");

		for(int j=0; j<100; ++j)//TODO:: why 100???
			tau=DirectMethodStep();
	}else
	{//tau-leaping branch

		//find a reaction to be fired
		float tau2;
		int k2Fire;
		
#ifdef INTERNAL_TIMING
		CudaEventTimer timer3;
		timer3.Start();
#endif
		
		GetReaction2Fire(m_StateData.m_rPropNonCritPartSum, tau2, k2Fire);
		
#ifdef INTERNAL_TIMING
		g_fr+=timer3.ElapsedTime();
#endif
		tau=min(tau1, tau2);
		++i;
		//	mylog<<i<<" "<<m_tau<<" "<<tau<<endl;
		
		//set times to be fired for a reaction
		
#ifdef INTERNAL_TIMING
		timer3.Start();
#endif
		//thrust::fill(m_k.begin(), m_k.end(), 0);
		SetNumberOfFiring(k2Fire!=ChannelsCount()&&tau2<tau1, k2Fire, tau);
		
#ifdef INTERNAL_TIMING
		g_sr+=timer3.ElapsedTime();
#endif

		//fire reactions
		
#ifdef INTERNAL_TIMING
		CudaEventTimer timer4;
		timer4.Start();
#endif
		FireReactions();
		
#ifdef INTERNAL_TIMING
		g_ft+=timer4.ElapsedTime();
#endif
		m_tau+=tau;

#ifdef EXTRA_CHECK
		int res=thrust::count_if(m_reactantsPopulation.begin(), m_reactantsPopulation.end(), IsNegative());
		if(res>0)
		{
			cerr<<"Negative population occured"<<endl;
			exit(EXIT_FAILURE);
		}
#endif
		
	}
	/*printf("number of firings: \n");
		print_vector<int, int> (m_k); 
		printf("\n");

	printf("number of reactants: \n");
		print_vector<int, int> (m_reactantsPopulation);
		printf("\n\n");*/

//	cout<<"e\n";

	return tau;
	
}

void Model::TransformMatrix(CSCMatrix const *reactionsDescCSC)
{
	m_ConstData.Init(reactionsDescCSC);
}

void Model::InitSeeds()
{

	CreateRandomSeedsUint4();
	
}

float Model::DirectMethodStep()
{
	float tau2;
	int k2Fire;
	
#ifdef INTERNAL_TIMING
	CudaEventTimer timer;
	timer.Start();
#endif
	GetReaction2Fire(m_StateData.m_rPropPartSum, tau2, k2Fire);
	
#ifdef INTERNAL_TIMING
	g_fr+=timer.ElapsedTime();
#endif

	thrust::fill(m_StateData.m_k.begin(), m_StateData.m_k.end(), 0);
	//thrust::host_vector<char> auxK(1, k2Fire);
	//thrust::copy(auxK.begin(), auxK.end(), m_k.begin()+k2Fire);

	FireReactions();
	m_tau+=tau2;

	return tau2;

}

void Model::FireReactions()
{
	m_ConstData.m_reactionsDesc->FireReactions(m_StateData.m_k, m_StateData.m_reactantsPopulation);
}

class NumberOfFiringFn
{
	bool m_fireCritical;
	int m_k2Fire;
	float const *m_prop;
	float const m_tau;
	uint4 *m_seeds;
public:
	NumberOfFiringFn(bool fireCritical, int k2Fire, vec_float_t const &prop, float const tau, thrust::device_vector<uint4> &seeds):
	  m_fireCritical(fireCritical), m_k2Fire(k2Fire), m_prop(thrust::raw_pointer_cast(&prop[0])), m_tau(tau), m_seeds(thrust::raw_pointer_cast(&seeds[0])){}

	inline __device__
	int operator()(bool isCritical, int reactionId)
	{
		return Poisson(reactionId)*(!isCritical)+
			(reactionId==m_k2Fire)*m_fireCritical*isCritical;
	}
private:
	inline __device__
	int Poisson(int reactionId)
	{
		//uint4 &globSeed=m_seeds[reactionId];
		uint4 locSeed=m_seeds[reactionId];
		const float prop=m_prop[reactionId];

		const float mean=prop*m_tau;

		const float lam=1.f/exp(mean);
		const int UL=300;

//		const int BLOCK_SIZE=blockDim.x*blockDim.y*blockDim.z;

//		__shared__ uint4 s_Seed[1024];

		//const int inBlockId=blockDim.x*threadIdx.y+threadIdx.x;
		//s_Seed[inBlockId]=m_seeds[reactionId];
		
		

		if(mean<80)
		{//return a Poisson-distribute RN
		//	*((int *)0)=5;
		
			float val=1.f;
			for(int n=1;n<UL;++n)
			{
				//return lam;
				val*=HybridTausRng(&locSeed);
				if(val<lam)
				{
					m_seeds[reactionId]=locSeed;
					return n-1;
				}
			}
			m_seeds[reactionId]=locSeed;
			return UL;
		}
		else
		{//return a normal-distribute approximation
			return abs(boxMuller(mean, mean, m_seeds[reactionId], m_seeds[2*reactionId]));
		}
		
	}
};

class PropSorter
{
public:
	__device__ inline
	bool operator()(int lhs, int rhs)const
	{
		return lhs<rhs;
	}
};

#ifdef _DEBUG
thrust::pair<vec_int_t::const_iterator, vec_int_t::const_iterator> 
#else
void
#endif
Model::SetNumberOfFiring(bool fireCritical, int k2Fire, float tau)
{
/*	thrust::sort(m_sortedIndices.begin(), m_sortedIndices.end(), 
		PropSorter());*/
	thrust::counting_iterator<int> cit(0);
	thrust::transform(m_StateData.m_criticalClassifier.begin(),m_StateData.m_criticalClassifier.end(),cit, m_StateData.m_k.begin(),
		NumberOfFiringFn(fireCritical, k2Fire, m_StateData.m_reactionsPropensity, tau, m_StateData.m_seeds));

#ifdef _DEBUG
	return thrust::minmax_element(m_StateData.m_k.begin(), m_StateData.m_k.end());
#endif
}


#ifdef INTERNAL_TIMING
CudaEventTimer g_timer1;
#endif
void Model::GetReaction2Fire(vec_float_t const &rPropPartSum, float &tauCrit, int &k2Fire)const
{
	//g_timer1.Start();
	tauCrit=std::numeric_limits<float>::max();
	//MyTime::Time_t tb=MyTime::CTime();
	
//	g_timer1.Start();
	float a0c=rPropPartSum.back();//ValueFromDevVector(rPropPartSum.end()-1);
	//g_frt+=g_timer1.ElapsedTime();
	
	//g_timer1.Start();
	float const rnd1=rand()/static_cast<float>(RAND_MAX);
	//g_frr+=g_timer1.ElapsedTime();
	
	//g_timer1.Start();
	thrust::device_vector<float>::const_iterator ub=
		thrust::upper_bound(rPropPartSum.begin(), rPropPartSum.end(), a0c*rnd1);
	//g_fr1+=g_timer1.ElapsedTime();
	
	k2Fire=thrust::distance(rPropPartSum.begin(), ub);
	if(k2Fire==rPropPartSum.size())
		return;

	//assert(k2Fire)

	/*float valPrev=*ub;
	assert(a0c*rnd1<valPrev);
	float valNext=*(ub+1);
	assert(a0c*rnd1>=valNext);*/

	//TODO: enable!!
	//if(&rPropPartSum==&m_StateData.m_rPropNonCritPartSum)
		//assert(m_StateData.m_criticalClassifier[k2Fire]==true);//check if we really found a critical reaction...

	//timer.Start();
	float const rnd2=rand()/static_cast<float>(RAND_MAX);
	//g_frr+=timer.ElapsedTime();
	assert(rnd2>=0.f&&rnd2<=1.f);

	//g_timer1.Start();
	//tauCrit=1.f/ValueFromDevVector(m_reactionsPropensity.begin()+k2Fire)*log(1.f/rnd2);
	tauCrit=1.f/m_StateData.m_reactionsPropensity[k2Fire]*log(1.f/rnd2);
	//g_fr2+=g_timer1.ElapsedTime();	
}

void Model::LogStatistics(std::ofstream &logF, vector<char> const &r2Track, int nIter, stat_t &statistics)const
{
//	assert(logF.is_open());
	thrust::host_vector<int> population=m_StateData.m_reactantsPopulation;

	//logF<<"#"<<m_tau<<endl;
	//logF<<m_tau<<"\t";
	int rN=0;
	for(size_t i=0; i<population.size(); ++i)
	{
		if(!r2Track[i])
			continue;
		//logF<<population[i]<<"\t";
		//statistics[rN][nIter]=population[i];
		statistics[m_tau].push_back(population[i]);
		++rN;
	}
	//logF<<endl;
}

void StatHeader(ofstream &logF)
{
	//open a log file
	unlink(g_logFN.c_str());
	logF.open(g_logFN.c_str());
	if(!logF.is_open())
	{
		cerr<<"Can not open the file: "<<g_logFN.c_str();
		exit(EXIT_FAILURE);
	}
	//logF<<"time\tS1\tS2\tS3"<<endl;
	
}

double g_totalProp=0,
		g_totalTau=0,
		g_totalFr=0,
		g_totalSr=0,
		g_totalFt=0;

float Model::PerformCalc(float endTau, float eps, vector<char> const &r2Track, stat_t &statistics, int maxIterations)
{
	ofstream logF;
	//if(g_doLogging)
		//StatHeader(logF);
	
	int i=0;
	double prevTau=-1.;
	bool normalExit=true;
	
#ifdef INTERNAL_TIMING
	CudaEventTimer timer;
#endif

	float totalTime=0.f;
	int iterationNumber=0;
	

	while(m_tau<endTau)
	{//main cycle...
		iterationNumber++;
		if(iterationNumber>maxIterations)
			break;
		if(g_doLogging)
			LogStatistics(logF, r2Track, iterationNumber-1, statistics);

		//might be uncommented for progress output
		if(m_tau-prevTau>endTau/100)
		{
	//		printf("step #%e, current time: %f\n", i, m_tau);
			prevTau=m_tau;
		}
		
		//preforming one step of the algorithm

		
#ifdef INTERNAL_TIMING
		timer.Start();
#endif
		float tau=Step(eps);

		
#ifdef INTERNAL_TIMING
		g_step+=timer.ElapsedTime();
#endif

		totalTime+=tau;

#ifdef EXTRA_CHECK
		if(thrust::none_of(m_k.begin(), m_k.end(), thrust::identity<int>()))
		{
			normalExit=false;
			break;
		}
#endif
		++i;
	}

	//benchmarking statistics
	//if(g_doLogging)
	{

#ifdef INTERNAL_TIMING
		g_totalProp+=g_prop;
		g_totalTau+=g_tau;
		g_totalFr+=g_fr;
		g_totalSr+=g_sr;
		g_totalFt+=g_ft;
#endif

		//cout<<endl;
	//	cout<<"Done, "<<i<<" iterations made"<<endl;
		//cout<<"Total real time: "<<totalTime<<endl;
		if(!normalExit)
			cout<<"Bad Exit"<<endl;
#ifdef INTERNAL_TIMING
	/*	cout<<"For propensity calculation: "<<g_prop/1000<<" s"<<endl;
		//cout<<"\tfor scanning: "<<g_scanTime/1000.f<<" s"<<endl;
		cout<<"For time step calculation: "<<g_tau/1000<<" s"<<endl;
		cout<<"For finding critical reaction: "<<g_fr/1000<<" s"<<endl;
		//cout<<"\taux operation "<<g_frt/1000<<" s"<<endl;
		//cout<<"\tStage1 "<<g_fr1/1000<<" s"<<endl;
		//cout<<"\tStage2 "<<g_fr2/1000<<" s"<<endl;
		//cout<<"\trandoms "<<g_frr/1000<<" s"<<endl;
		cout<<"For setting number of reactions to fire: "<<g_sr/1000<<" s"<<endl;
		cout<<"For firing reactions: "<<g_ft/1000<<" s"<<endl;
		cout<<"-------****************---------"<<endl;
		cout<<"Summary: "<<
			(g_prop+g_tau+g_fr+g_ft+g_sr)/1000<<" s"<<endl;
		cout<<endl;
		cout<<"For overal step calculation: "<<g_step/1000<<" s"<<endl;*/
#endif

		//cout<<"test2\n";
		
	}
	return totalTime;
}

__host__ void Model::CreateRandomSeedsUint4()
{
	thrust::host_vector<uint4> h_seeds;
	try
	{
		//printf("\tCreateRandomSeedsUint4: host memory allocating...\n");
		h_seeds.resize(ChannelsCount()*2);
	}
	catch(std::exception const &ec)
	{
		fprintf(stderr, "%s while trying to allocate %d Mbytes of memory\n", ec.what(), ChannelsCount()*sizeof(uint4)/(1024*1024));
		exit(EXIT_FAILURE);
	}
	
	//printf("\tCreateRandomSeedsUint4: host random numbers generating...\n");
	if(!FillRandoms(static_cast<char *>(static_cast<void *>(&h_seeds[0])), 
		sizeof(uint4)/sizeof(char)*ChannelsCount()))
	{
		fprintf(stderr, "Can not allocate random seeds");
		exit(EXIT_FAILURE);
	}

	//printf("\tCreateRandomSeedsUint4: memory copying...\n");
	thrust::copy(h_seeds.begin(), h_seeds.end(), m_StateData.m_seeds.begin());
	//printf("\tCreateRandomSeedsUint4: done\n");
}

Model *Model::Clone()const
{
	Model *res=new Model(ChannelsCount(), ReactantsCount());
	assert(m_ConstData.m_reactionsDesc.get());
	m_ConstData.CopyTo(res->m_ConstData);
	res->m_StateData=m_StateData;
	res->m_tau=m_tau;

	return res;
}

