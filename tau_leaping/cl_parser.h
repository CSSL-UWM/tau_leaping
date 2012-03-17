#pragma once
#include <unordered_set>
#include <string>

//parsing similar command-line arguments for 
//tau-leaping and "improved" direct methods

#ifdef _WIN32
namespace std{
	using namespace tr1;
}
#endif

struct CLArguments
{
	std::string fileName;
	float endTime;
	float eps;
	std::unordered_set<std::string> reactants2Track;
	std::string logFN;
	int numParallel;
	int cudaDeviceNumber;
	int ai;
	int loggingIntervals;
	int maxIterations;

};

void ParseCL(int argc, char *argv[],
			 CLArguments &args);
