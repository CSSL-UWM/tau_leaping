#include "model.h"
#include <time.h>
#include <cassert>
#include "cl_parser.h"
#include <unordered_set>
#include <stdexcept>

#include "CSCMatrix.h"

#include "MyTime.h"

#ifdef _WIN32
namespace std{
	using namespace tr1;
}
#endif

cusparseHandle_t g_handle;

void test();
using std::endl;
using std::cerr;
using std::cout;
using std::string;
using std::unordered_set;
using std::vector;


std::string g_logFN="stat.txt";
bool g_doLogging=false;


extern double g_totalProp,
		g_totalTau,
		g_totalFr,
		g_totalSr,
		g_totalFt; 


int main(int argc, char* argv[])
{
	/*bool fromFile=false;
	const char *fn=NULL;
	float endTime=1.0f;
	float eps=0.03f;
	unordered_set<string> reactants2Track;*/

	cudaSetDevice(1);

	CLArguments args;

	//parsing command line arguments
	ParseCL(argc, argv, args);
	g_logFN=args.logFN;

	//filling out the set of reactants to be tracked
	if(!args.reactants2Track.empty())
	{
		cout<<"Tracking reactants: "<<endl;
		for(unordered_set<string>::const_iterator it=args.reactants2Track.begin(); it!=args.reactants2Track.end(); ++it)
		{
			cout<<*it<<endl;
		}	
	}

	cout<<"Processing file "<<args.fileName<<endl;
	cout<<"Finish time: "<<args.endTime<<endl;
	cout<<"Running "<<args.numParallel<<" system(s) in parallel"<<endl;
	cout<<"Logging intervals: "<<args.loggingIntervals<<endl;
	cout<<"Iterations limit: "<<args.maxIterations<<endl;

#ifdef EXTRA_CHECK
	cout<<"Extra checking on"<<endl;
#else
	cout<<"Extra checking off"<<endl;
#endif

	cout<<endl<<"logging statistics into "<<g_logFN<<" file"<<endl;

	time_t tm;
	srand(static_cast<unsigned int>(time(&tm)));
	//std::auto_ptr<Model> model;
	vector<Model *> models(1);
	vector<string> reactantNames;
	//std::auto_ptr<CSCMatrix> reactionsDescCSC;
	//if(fromFile)
	{//reading model from a StochKit xml file
		//assert(!args.fileName.empty()&&"there must be a valid file name");
		try
		{
			if(args.fileName.empty())
				throw std::runtime_error("Random systems are no longer supported");
			CSCMatrix *ptr;
			models[0]=Model::ReadXML(args.fileName.c_str(), reactantNames);
			//reactionsDescCSC.reset(ptr);
		}catch(std::exception const &ec)
		{
			cerr<<ec.what()<<endl;
			exit(EXIT_FAILURE);
		}
	}

	//storing reactants to be tracked into an appropriate structure
	vector<char> r2Track(models[0]->ReactantsCount(), 0);
	{
		bool trackAll=args.reactants2Track.size()==1&&
			(args.reactants2Track.find("-1")!=args.reactants2Track.end()||
			 args.reactants2Track.find("all")!=args.reactants2Track.end())?
			 true:false;
		for(int i=0; i<models[0]->ReactantsCount(); ++i)
		{
			string const& rName=reactantNames[i];
			if(args.reactants2Track.find(rName)!=args.reactants2Track.end()||trackAll)
			{
				r2Track[i]=1;
			}
		}
	}

	//no longer need reactants' names
	vector<string>().swap(reactantNames);

	//assert(models[0]->m_ConstData.m_reactionsDesc.get());

	//if we have multiple systems running simultaneously, 
	//copy other systems from the first one
	for(int i=1; i<args.numParallel; ++i)
	{
		models.push_back(models[0]->Clone());
	}

	//preparing for calculation...
	for(int i=0; i<args.numParallel; ++i)
	{
		models[i]->InitSeeds();
	}

	std::vector<stat_t> statistics(args.numParallel);

	//do the calculation
	MyTime::Time_t tb=MyTime::CTime();

	float compTime=0.f;
	for(int i=0; i<args.numParallel; ++i)
	{
		if(i%(args.numParallel/10)==0)
			cout<<"realization: "<<i<<endl;
		
		compTime+=models[i]->PerformCalc(args.endTime, args.eps, r2Track, statistics[i], args.maxIterations);
	}

	{
		cout<<"For propensity calculation: "<<g_totalProp/1000<<" s"<<endl;
		//cout<<"\tfor scanning: "<<g_scanTime/1000.f<<" s"<<endl;
		cout<<"For time step calculation: "<<g_totalTau/1000<<" s"<<endl;
		cout<<"For finding critical reaction: "<<g_totalFr/1000<<" s"<<endl;
		//cout<<"\taux operation "<<g_frt/1000<<" s"<<endl;
		//cout<<"\tStage1 "<<g_fr1/1000<<" s"<<endl;
		//cout<<"\tStage2 "<<g_fr2/1000<<" s"<<endl;
		//cout<<"\trandoms "<<g_frr/1000<<" s"<<endl;
		cout<<"For setting number of reactions to fire: "<<g_totalSr/1000<<" s"<<endl;
		cout<<"For firing reactions: "<<g_totalFt/1000<<" s"<<endl;
		cout<<"-------****************---------"<<endl;
		cout<<"Summary: "<<
			(g_totalProp+g_totalTau+g_totalFr+g_totalFt+g_totalSr)/1000<<" s"<<endl;
		cout<<endl;
		//cout<<"For overal step calculation: "<<g_step/1000<<" s"<<endl;
	}
	
	cout<<"It took "<<MyTime::ElapsedTime(tb, MyTime::CTime())/1000.f<<" s to complete calculations"<<endl;
	cout<<"Computational time, s: "<<compTime<<endl;

	StepStat_t stepStat;

	cout<<"Computing statistics at time steps..."<<endl;
	AllStat2StepStat(statistics, stepStat, args.endTime/args.loggingIntervals);
	cout<<" done"<<endl;

	StdDevStat_t meanStdDev;
	StepStat2StdDev(stepStat, meanStdDev);

	StatToFile(args.logFN, meanStdDev);

	for(size_t i=0; i<models.size(); ++i)
	{
		delete models[i];
	}
		
	return 0;
}

