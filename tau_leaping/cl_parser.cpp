//parsing similar command-line arguments for 
//tau-leaping and "improved" direct methods

#include "cl_parser.h"
#include <iostream>
#include <vector>
#include <limits>

#ifndef _WIN32
#include <string.h>
#endif

using std::unordered_set;
using std::string;
using std::cout;
using std::endl;
using std::vector;

void ParseCL(int argc, char *argv[],
			 CLArguments &args)
{
    args.eps=0.03f;
	args.numParallel=1;
	args.cudaDeviceNumber=0;
	args.loggingIntervals=1;
	args.maxIterations=std::numeric_limits<int>::max();
	for(int a=1; a<argc; ++a)
	{
		if(strcmp(argv[a],"-f")==0)//input file name
		{
			args.fileName=argv[++a];
			//fromFile=true;
		}

		if(strcmp(argv[a],"-n")==0)//how many systems to run in parallel
		{
			const char *ch=argv[++a];
			args.numParallel=atoi(ch);
			continue;
		}

		if(strcmp(argv[a],"-a")==0)//stop after this number of iterations
		{
			const char *ch=argv[++a];
			args.maxIterations=atoi(ch);
			continue;
		}

		if(strcmp(argv[a],"-d")==0)//which cuda device to use
		{
			const char *ch=argv[++a];
			args.cudaDeviceNumber=atoi(ch);
			continue;
		}

		if(strcmp(argv[a],"-b")==0)
		{
			args.logFN=argv[++a];
			/*g_logFN+="_";
			time_t ctm;
			time(&ctm);
			// Convert time to struct tm form
			tm *newTime = localtime( &ctm );
			g_logFN+=asctime(newTime);
			g_logFN=g_logFN.substr(0, g_logFN.length()-1);
			for(int i=0; i<g_logFN.size();++i)
			{
			if(g_logFN[i]==':')
			g_logFN[i]='_';
			}*/

			args.logFN+=".txt";
		}
		if(strcmp(argv[a],"-t")==0)//end time
		{
			const char *val=argv[++a];
			args.endTime=(float)atof(val);
		}
		if(strcmp(argv[a],"-e")==0)//epsilon
		{
			const char *val=argv[++a];
			args.eps=(float)atof(val);
			continue;
		}

		if(strcmp(argv[a],"-i")==0)//logging intervals
		{
			const char *val=argv[++a];
			args.loggingIntervals=atoi(val);
			continue;
		}
		if(strcmp(argv[a], "-r")==0)//trackable reactants
		{
			string str;

			for(int j=a+1; j<argc; ++j)
			{
				const char *val=argv[j];
				if(val[0]=='-')
					break;
				str.append(val);
				str.append(" ");
			}

			//			cmatch res;
			//str = "<h2>Egg prices</h2>";
			//regex rx("<h(.)>([^<]+)");
			/*std::regex pattern("\\w+");
			std::sregex_token_iterator end;
			for (std::sregex_token_iterator i(str.begin(),
				str.end(), pattern);
				i != end;
			++i)
			{
				//std::cout << i->str() << std::endl;
				reactants2Track.insert(i->str());
			}*/
			const char *splitter=" ,";
			//vector<char> cstr(str.begin(), str.end());
			//char *cstr=str.c_str();
			char *pch = strtok (&str[0], splitter);
			while (pch != NULL)
			{
			    args.reactants2Track.insert(pch);
				++a;
			    pch = strtok (NULL, splitter);
			}
		}
	}
}
