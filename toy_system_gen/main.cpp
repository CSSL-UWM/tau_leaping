#include <cstdlib>
#include <rapidxml/rapidxml.hpp>
#pragma warning (push)
#pragma warning (disable: 4100)
#include <rapidxml/rapidxml_print.hpp>
#pragma warning (pop)
#include <iostream>
#include <string>
#include <cassert>
#include <ctime>
#include <sstream>
#include <unordered_set>
#include <vector>
#include <cstring>

#ifndef _WIN32
#include <tr1/random>
#else
#include <random>
#endif

using namespace rapidxml;
using std::cout;
using std::string;
using std::endl;
using std::cerr;
using std::stringstream;
using std::vector;
using std::swap;

namespace std{
	using namespace tr1;
}


using std::uniform_real;
using std::uniform_int;
using std::unordered_set;

std::mt19937 g_rEng;//Mersenne twister std::tr1 RNG


const string stochKitXMLComment="<!-- Stochkit 2.0 text input format\n"
"\n"
"reference:\n"
"M. Rathinam, L.R. Petzold, Y. Cao, and D.T. Gillespie. J Chem Phys 119, 12784 (2003).\n"
"-->\n";

const string myComment="<!-- Input data for testing StochKit vs. GPU implementation. Ivan Komarov 2012-->\n";

enum SystemType{NormalModel, CyclicChainModel, ColloidalAggregationModel};

struct Params
{
	int channels;
	int reactants;
	int foReactions;
	int soReactions;
	int soReactionsSpec;
	int minSpPpl;
	int maxSpPpl;
	SystemType type;
};

int pi(int k1, int k2)
{
	return	(k1 + k2)/2*(k1 + k2 + 1) + k2;
}

struct Reaction
{
	int r1Ind;
	int r1Stoich;
	int r2Ind;
	int r2Stoich;
	
	int p1Ind;
	int p1Stoich;
	int p2Ind;
	int p2Stoich;

	int type;

	int id;

	/*operator size_t()const
	{
		return pi(r1Ind, r2Ind);
	}*/


	bool operator==(Reaction const&r)const
	{
		if(type!=r.type)
			return false;
		if(type==1||type==3)
		{
			if(r1Ind!=r.r1Ind||p1Ind!=r.p1Ind)
				return false;
		}else
		{
			assert(type==2);
			if(!((r1Ind==r.r1Ind&&r2Ind==r.r2Ind)||(r1Ind==r.r2Ind&&r2Ind==r.r1Ind)))
				return false;
		}
		return true;
	}
	std::string ToDesc()const
	{
		assert(Check());
		stringstream res;
		//res<<'\'';
		//res.putback('>');
		const char *arrow="to";
		if(r1Stoich==2)
			res<<"2 * ";
		res<<"S"<<r1Ind;

		if(r2Stoich!=0)
		{
			res<<" + ";
			if(r2Stoich==2)
				res<<" 2 * ";
			else
				res<<" ";
			res<<"S"<<r2Ind;
		}

		res<<" "<<arrow;

		if(p1Stoich==0&&p2Stoich==0) 
			res<<" null";
		else
		{
			if(p1Stoich==2)
				res<<" 2 * ";
			else
				res<<" ";
			res<<"S"<<p1Ind;

			if(p2Stoich!=0)
			{
				res<<" + ";
				if(p2Stoich==2)
					res<<" 2 * ";
				res<<"S"<<p2Ind;
			}
		}
		//res<<'\'';
		return res.str();
	}

	bool Check()const
	{
		if(r1Ind==r2Ind||r1Ind==p1Ind||r1Ind==p2Ind)
			return false;

		if(r2Ind!=0&&(r2Ind==p1Ind||r2Ind==p2Ind))
			return false;

		if(p1Ind!=0&&p1Ind==p2Ind)
			return false;

		if(r1Ind==0)
			return false;

		if(r1Stoich!=1&&r1Stoich!=2)
			return false;

		if(r2Ind!=0)
			if(r2Stoich!=1&&r2Stoich!=2)
				return false;

		if(p1Ind!=0)
			if(p1Stoich!=1&&p1Stoich!=2)
				return false;

		if(p2Ind!=0)
			if(p2Stoich!=1&&p2Stoich!=2)
				return false;

		if(!(type==1||type==2||type==3))
			return false;

		if(id<=0)
			return false;

		return true;
	}
};

namespace std {

   template <>
   struct hash<Reaction> : public unary_function<Reaction, size_t>
   {
       size_t operator()(const Reaction& r) const
       {
           return pi(r.r1Ind, r.r2Ind);
       }
   };
}
/*
namespace std
{
	inline size_t hash<Gateway<T> >::operator()(const Gateway<T> & gT) const 
    {
          return gT.getCore()->hash();          
    }
}*/

void PrintParams(Params const &params)
{
	cerr<<"Generating input data for "<<params.channels<<" channels and "<<
		params.reactants<<" reactants"<<endl;
	cerr<<"System type: ";
	switch (params.type){
	case NormalModel:
		cerr<<"Normal"<<endl;break;
	case CyclicChainModel:
		cerr<<"Weak"<<endl;break;
	case ColloidalAggregationModel:
		cerr<<"Strong"<<endl; break;
	default:
		cerr<<"Unknown reaction type";
		exit(1);
	}
	cerr<<"Making "<<params.foReactions<<" of first type of reactions,"<<endl
		<<"\t "<<params.soReactions<<" of second type of reactions,"<<endl
		<<"\t "<<params.soReactionsSpec<<" of third type of reactions"<<endl;
	cerr<<"Reactant's population varies from "<<params.minSpPpl<<" to "<<params.maxSpPpl<<endl;
}

void ExitOnWrongOptions()
{
	cerr<<"Use -c <channels number> -r <reactants number> options at least\n";
	exit(EXIT_FAILURE);
}

void ParseArguments(int argc, char**argv, Params &params)
{
	if(argc<4)
	{
		ExitOnWrongOptions();
	}

	Params res={0,0,-1,-1,-1, 0, 10000, NormalModel};
	for(int i=1; i<argc; ++i)
	{
		const char *key=argv[i];
		if(strcmp(key,"-c")==0)
		{
			if(i+1>=argc)
				ExitOnWrongOptions();
			res.channels=atoi(argv[++i]);
			if(res.channels<=0)
			{
				cerr<<"Incorrect number of channels"<<endl;
				exit(EXIT_FAILURE);
			}
			continue;
		}
		if(strcmp(key,"-r")==0)
		{
			if(i+1>=argc)
				ExitOnWrongOptions();
			res.reactants=atoi(argv[++i]);
			if(res.reactants<=3)
			{
				cerr<<"Incorrect number of reactants"<<endl;
				exit(EXIT_FAILURE);
			}
			continue;
		}

		if(strcmp(key,"-t1")==0)
		{
			if(i+1>=argc)
				ExitOnWrongOptions();
			res.foReactions=atoi(argv[++i]);
			if(res.foReactions<=0||res.foReactions>res.channels)
			{
				cerr<<"Incorrect number of reactants of the first type"<<endl;
				exit(EXIT_FAILURE);
			}
			continue;
		}

		if(strcmp(key,"-t2")==0)
		{
			if(i+1>=argc)
				ExitOnWrongOptions();
			res.soReactions=atoi(argv[++i]);
			if(res.soReactions<=0||res.soReactions>res.channels)
			{
				cerr<<"Incorrect number of reactants of the second type"<<endl;
				exit(EXIT_FAILURE);
			}
			continue;
		}

		if(strcmp(key,"-t3")==0)
		{
			if(i+1>=argc)
				ExitOnWrongOptions();
			res.soReactionsSpec=atoi(argv[++i]);
			if(res.soReactionsSpec<=0||res.soReactionsSpec>res.channels)
			{
				cerr<<"Incorrect number of reactants of the third type"<<endl;
				exit(EXIT_FAILURE);
			}
			continue;
		}

		if(strcmp(key,"-mnp")==0)
		{
			if(i+1>=argc)
				ExitOnWrongOptions(); 
			res.minSpPpl=atoi(argv[++i]);
			
			continue;
		}

		if(strcmp(key,"-mxp")==0)
		{
			if(i+1>=argc)
				ExitOnWrongOptions();
			res.maxSpPpl=atoi(argv[++i]);
			
			continue;
		}

		if(strcmp(key,"-w")==0)
		{
			res.type=CyclicChainModel;
			
			continue;
		}

		if(strcmp(key,"-s")==0)
		{
			res.type=ColloidalAggregationModel;
			
			continue;
		}

		cerr<<"Unknown parameter: "<<key<<endl;
		exit(EXIT_FAILURE);
	}
	
	if((res.foReactions>=0?res.foReactions:0+
		res.soReactions>=0?res.soReactions:0+
		res.soReactionsSpec>=0?res.soReactionsSpec:0)>res.channels)
	{
		cerr<<"Inconcistent number of reactions of differnt types\n";
		exit(EXIT_FAILURE);
	}

	
	if(res.foReactions==-1)
	{
		int diff=res.channels-(res.soReactions>=0?res.soReactions:0+
			res.soReactionsSpec>=0?res.soReactionsSpec:0);
		if(res.soReactions>=0&&res.soReactionsSpec>=0)
			res.foReactions=diff;
		else
		{
			uniform_int<int> rnd(0, diff);
			res.foReactions=rnd(g_rEng);
		}
	}

	if(res.soReactions==-1)
	{
		int diff=res.channels-(res.foReactions>=0?res.foReactions:0+
			res.soReactionsSpec>=0?res.soReactionsSpec:0);
		if(res.foReactions>=0&&res.soReactionsSpec>=0)
			res.soReactions=diff;
		else
		{
			uniform_int<int> rnd1(0, diff);
			res.soReactions=rnd1(g_rEng);
		}
	}

	if(res.soReactionsSpec==-1)
	{
		res.soReactionsSpec=res.channels-(res.foReactions+res.soReactions);
	}

	if(res.foReactions+res.soReactions+res.soReactionsSpec!=res.channels)
	{
		cerr<<"Inconcistent number of reactions of differnt types\n";
		exit(EXIT_FAILURE);
	}

	if(res.type==CyclicChainModel){
		res.channels=res.reactants;
	}

	if(res.type==ColloidalAggregationModel){
		res.channels=res.reactants*res.reactants/2;
	}


	params=res;
}

void AddRandRnConstants(xml_document<>& doc, xml_node<> *paramList, int rn)
{
	assert(paramList);

	char buff[250];

	xml_node<> *param = doc.allocate_node(node_element, "Parameter");

	sprintf(buff,"c%d",rn+1);
	char *nodeVal = doc.allocate_string(buff);        // Allocate string and copy name into it
	xml_node<> *node = doc.allocate_node(node_element, "Id", nodeVal);
	param->append_node(node);

	uniform_real<float> unif(0.f, 1.f);//std::tr1 uniform distribution

	float rConst=unif(g_rEng);
#ifdef __GNUG__
	rConst/=static_cast<float>(std::numeric_limits<unsigned>::max());
#endif 
	sprintf(buff,"%f",rConst);
	//cerr<<rConst<<" "<<buff<<endl;
	if(rConst<0.f||rConst>1.f)
	{
		cerr<<"Check (standard) uniform_real implementation. There is a bug in some GCC libraries"<<endl;
		exit(EXIT_FAILURE);
	}
	nodeVal = doc.allocate_string(buff);        // Allocate string and copy name into it
	node = doc.allocate_node(node_element, "Expression", nodeVal);
	param->append_node(node);

	paramList->append_node(param);
}


xml_node<>* MakeSrNode(xml_document<> &doc, int rInd, int rStoich)
{
	char buff[250];

	xml_node<> *srNode = doc.allocate_node(node_element, "SpeciesReference");
	sprintf(buff,"S%d",rInd);
	char *nodeVal = doc.allocate_string(buff);        // Allocate string and copy name into it
	xml_attribute<> *attr = doc.allocate_attribute("id", nodeVal);
	srNode->append_attribute(attr);

	sprintf(buff,"%d",rStoich);
	nodeVal = doc.allocate_string(buff);        // Allocate string and copy name into it
	attr = doc.allocate_attribute("stoichiometry", nodeVal);
	srNode->append_attribute(attr);
	return srNode;
}
void AddRandRns(xml_document<>& doc, xml_node<> *paramList, Reaction const &reaction)
{
	assert(paramList);
	assert(reaction.Check());

	char buff[250];

	xml_node<> *param = doc.allocate_node(node_element, "Reaction");

	sprintf(buff,"R%d",reaction.id);
	char *nodeVal = doc.allocate_string(buff);        // Allocate string and copy name into it
	xml_node<> *node = doc.allocate_node(node_element, "Id", nodeVal);
	param->append_node(node);

	sprintf(buff,"%s",reaction.ToDesc().c_str());
	nodeVal = doc.allocate_string(buff);        // Allocate string and copy name into it
	node = doc.allocate_node(node_element, "Description", nodeVal);
	param->append_node(node);

	node = doc.allocate_node(node_element, "Type", "mass-action");
	param->append_node(node);

	sprintf(buff,"c%d",reaction.id);
	nodeVal = doc.allocate_string(buff);        // Allocate string and copy name into it
	node = doc.allocate_node(node_element, "Rate", nodeVal);
	param->append_node(node);

	{
		xml_node<> *rNode = doc.allocate_node(node_element, "Reactants");
		param->append_node(rNode);

		rNode->append_node(MakeSrNode(doc, reaction.r1Ind, reaction.r1Stoich));

		if(reaction.r2Ind!=0)
		{
			assert(reaction.r2Stoich!=0);
			rNode->append_node(MakeSrNode(doc, reaction.r2Ind, reaction.r2Stoich));
		}

		xml_node<> *pNode = doc.allocate_node(node_element, "Products");
		param->append_node(pNode);

		if(reaction.p1Ind!=0)
		{
			assert(reaction.p1Stoich!=0);
			pNode->append_node(MakeSrNode(doc, reaction.p1Ind, reaction.p1Stoich));
		}

		if(reaction.p2Ind!=0)
		{
			assert(reaction.p2Stoich!=0);
			pNode->append_node(MakeSrNode(doc, reaction.p2Ind, reaction.p2Stoich));
		}
		
	}


	paramList->append_node(param);
}

void AddRandRtPopulations(xml_document<>& doc, xml_node<> *paramList, int rt, int minP, int maxP)
{
	assert(paramList);

	char buff[250];

	xml_node<> *param = doc.allocate_node(node_element, "Species");

	sprintf(buff,"S%d",rt+1);
	char *nodeVal = doc.allocate_string(buff);        // Allocate string and copy name into it
	xml_node<> *node = doc.allocate_node(node_element, "Id", nodeVal);
	param->append_node(node);

	sprintf(buff,"Species #%d",rt+1);
	nodeVal = doc.allocate_string(buff);        // Allocate string and copy name into it
	node = doc.allocate_node(node_element, "Description", nodeVal);
	param->append_node(node);

	uniform_int<int> unif(minP, maxP);//std::tr1 uniform distribution
	sprintf(buff,"%d",unif(g_rEng));
	nodeVal = doc.allocate_string(buff);        // Allocate string and copy name into it
	node = doc.allocate_node(node_element, "InitialPopulation", nodeVal);
	param->append_node(node);

	paramList->append_node(param);
}

int GetNextInd(std::vector<int> &unusedReactants, int fullCount)
{
	
	if(unusedReactants.empty())
	{
		return uniform_int<int>(1, fullCount)(g_rEng);
	}else
	{
		int i=uniform_int<int>(0, unusedReactants.size()-1)(g_rEng);
		swap(unusedReactants[i], unusedReactants[unusedReactants.size()-1]);
		int ind=unusedReactants.back();
		unusedReactants.pop_back();
		return ind;
	}
}

Reaction MakeReaction(int type, vector<int> &unusedReactants, int fullRC)
{
	assert(type==1||type==2||type==3);
	//uniform_int<int> rtDistr(0, unusedReactants.size()-1);
	Reaction res;
	if(type==1)
	{
		Reaction rn={GetNextInd(unusedReactants, fullRC),1, 0, 0, 0, uniform_int<int>(0,2)(g_rEng)};
		if(rn.p1Stoich)
		{
			do{
				rn.p1Ind=GetNextInd(unusedReactants, fullRC);
			}while(rn.p1Ind==rn.r1Ind);

			if(rn.p1Stoich==1)
			{
				rn.p2Stoich=uniform_int<int>(0,1)(g_rEng);
				if(rn.p2Stoich)
				{
					do{
						rn.p2Ind=GetNextInd(unusedReactants, fullRC);
					}while(rn.p2Ind==rn.r1Ind||rn.p2Ind==rn.p1Ind);
				}
			}
		}
		res=rn;
		res.type=1;
		
	}else
	if(type==2)
	{
		Reaction rn={GetNextInd(unusedReactants, fullRC), 1, 0, 1, 
			0, 1, 0, uniform_int<int>(0,1)(g_rEng)};

		do{
			rn.r2Ind=GetNextInd(unusedReactants, fullRC);
		}while(rn.r2Ind==rn.r1Ind);

		do{
			rn.p1Ind=GetNextInd(unusedReactants, fullRC);
		}while(rn.p1Ind==rn.r1Ind||rn.p1Ind==rn.r2Ind);

		if(rn.p2Stoich)
		{
			do{
				rn.p2Ind=GetNextInd(unusedReactants, fullRC);
			}while(rn.p2Ind==rn.r1Ind||rn.p2Ind==rn.r2Ind||rn.p2Ind==rn.p1Ind);
		}
		res=rn;
		res.type=2;
	}else
	if(type==3)
	{
		Reaction rn={GetNextInd(unusedReactants, fullRC), 2, 0, 0, 0, 1};

		do{
			rn.p1Ind=GetNextInd(unusedReactants, fullRC);
		}while(rn.p1Ind==rn.r1Ind);
		
		res=rn;
		res.type=3;
	}else
	{
		cerr<<"unknown type of reatcion!"<<endl;
		exit(EXIT_FAILURE);
	}
//	assert(res.Check());
	return res;
}
void AddReaction(int maxCycles, vector<int> &unusedReactants, int rc, 
				 unordered_set<Reaction> &reactionList, int &rId, int type)
{
	const char *err="There are not enough rectants to produce unique reactions";

	for(int i=0; i<maxCycles; ++i)
	{
		Reaction r=MakeReaction(type, unusedReactants, rc);
		unordered_set<Reaction>::const_iterator cit=reactionList.find(r);
		if(cit==reactionList.end())
		{
			r.id=rId++;
			assert(r.Check());
			reactionList.insert(r);
			return;
		}
	}

	cerr<<err<<endl;
	exit(EXIT_FAILURE);
		
}

void GenReactionsList(Params const &params, unordered_set<Reaction> &reasctionList)
{
	int const CYCLES=100;
	//uniform_real<float> distr(0.f, 1.f);
	vector<int> unusedReactants(params.reactants);
	reasctionList.clear();
	
	/*Reaction r={1,1};
	r.type=1;
	reasctionList.insert(r);
	reasctionList.insert(r);*/
	for(int i=0; i<params.reactants; ++i)
	{
		unusedReactants[i]=i+1;
	}
	int rId=1;
	
	for(int rn=0; rn<params.foReactions; ++rn)
	{
		AddReaction(CYCLES, unusedReactants, params.channels, reasctionList, rId, 1);
	}

	for(int rn=0; rn<params.soReactions; ++rn)
	{
		AddReaction(CYCLES, unusedReactants, params.channels, reasctionList, rId, 2);
	}

	for(int rn=0; rn<params.soReactionsSpec; ++rn)
	{
		AddReaction(CYCLES, unusedReactants, params.channels, reasctionList, rId, 3);
	}
	assert(rId==params.channels+1&&reasctionList.size()==static_cast<size_t>(params.channels));
	if(!unusedReactants.empty())
	{
		cerr<<"Not all reactants were used, "<<
			unusedReactants.size()<<" left. Too few channels(reactions)"<<endl;
		exit(EXIT_FAILURE);
	}
}

void GenCyclicChainReactionsList(Params const &params, unordered_set<Reaction> &reactionList){
	int i=1;
	for(; i<params.channels; ++i){
		Reaction rn={i, 1, 0, 0, //reactants
			i+1, 1, 0, 0, //products
			1,
			i
		};
		reactionList.insert(rn);
	}
	Reaction rn={i, 1, 0, 0, //reactants
			1, 1, 0, 0, //products
			1,
			i
	};
	reactionList.insert(rn);
}

void GenColloidalAggregationReactionsList(Params const &params, unordered_set<Reaction> &reactionList){
	int i=1;
	for(int n=1; n<=params.reactants/2; ++n){
		Reaction rn={n, 2, 0, 0, //reactants
			2*n, 1, 0, 0, //products
			3,
			i++
		};
		auto res=reactionList.insert(rn);
		assert(res.second);

		for(int m=n+1; m<=params.reactants-n; ++m){
			Reaction rn={n, 1, m, 1, //reactants
				m+n, 1, 0, 0, //products
				2,
				i++
			};
			auto res=reactionList.insert(rn);
			assert(res.second);
		}
	}

	for(int p=1; p<=params.reactants; ++p){
		for(int q=1; q<=p/2; ++q){
			if(q==p-q){
				Reaction rn={p, 1, 0, 0, //reactants
					q, 2, 0, 0, //products
					1,
					i++
				};
				auto res=reactionList.insert(rn);
				assert(res.second);
			}else{
				Reaction rn={p, 1, 0, 0, //reactants
					q, 1, p-q, 1, //products
					1,
					i++
				};
				auto res=reactionList.insert(rn);
				assert(res.second);
			}
		}
	}

}

void MakeXMLDoc(xml_document<>& doc, Params const &params)
{

	unordered_set<Reaction> reasctions;
	switch(params.type){
	case NormalModel:
		GenReactionsList(params, reasctions);
		break;
	case CyclicChainModel:
		GenCyclicChainReactionsList(params, reasctions);
		break;
	case ColloidalAggregationModel:
		GenColloidalAggregationReactionsList(params, reasctions);
		break;
	}

	

	assert(reasctions.size()==static_cast<size_t>(params.channels));

	char buff[255];
	xml_node<> *modelNode = doc.allocate_node(node_element, "Model");

	xml_node<> *node = doc.allocate_node(node_element, "Description", "Test input data");
	modelNode->append_node(node);

	sprintf(buff, "%d", params.channels);
	char *node_name = doc.allocate_string(buff);        // Allocate string and copy name into it
	node = doc.allocate_node(node_element, "NumberOfReactions", node_name);
	modelNode->append_node(node);

	sprintf(buff, "%d", params.reactants);
    node_name = doc.allocate_string(buff);        // Allocate string and copy name into it
	node = doc.allocate_node(node_element, "NumberOfSpecies", node_name);
	modelNode->append_node(node);

	cerr<<"generating reactions' constants..."<<endl;
	{
		xml_node<> *paramList = doc.allocate_node(node_element, "ParametersList");
		modelNode->append_node(paramList);
		for(int rn=0; rn<params.channels; ++rn)
		{
			AddRandRnConstants(doc, paramList, rn);
		}
	}

	cerr<<"generating reactions' list..."<<endl;
	{
		xml_node<> *reactionsList = doc.allocate_node(node_element, "ReactionsList");
		modelNode->append_node(reactionsList);


		//for(int rn=0; rn<params.channels; ++rn)
		for(unordered_set<Reaction>::const_iterator cit=reasctions.begin(); cit!=reasctions.end(); ++cit)
		{
			//Reaction r={1,2,4,1,3,1,5,1};
			AddRandRns(doc, reactionsList, *cit);
		}
	}

	cerr<<"generating spicies' populations..."<<endl;
	{
		xml_node<> *speciesList = doc.allocate_node(node_element, "SpeciesList");
		modelNode->append_node(speciesList);
		for(int rt=0; rt<params.reactants; ++rt)
		{
			AddRandRtPopulations(doc, speciesList, rt, params.minSpPpl, params.maxSpPpl);
		}
	}

	doc.append_node(modelNode);
}


int main(int argc, char **argv)
{
	try{
	//seed the RNG engine with current time for the release version
#ifndef _DEBUG
	g_rEng.seed(static_cast<unsigned long>(time(NULL)));
#else
	g_rEng.seed(2);//just to get more "interesting" input data
#endif


	//parsing command line
	Params params;
	ParseArguments(argc, argv, params);

	//parmeters printing
	PrintParams(params);

	//creating an XML document with random reactions 
	xml_document<> doc;
	MakeXMLDoc(doc, params);

	cerr<<"Writing XML...."<<endl;
	cout<<stochKitXMLComment<<endl<<myComment<<endl;
	

	if(params.type==CyclicChainModel)
		cout<<"<!--Weakly coupled (Cyclic Chain) system -->\n";
	if(params.type==ColloidalAggregationModel)
		cout<<"<!--Strongly coupled (Colloidal Aggregation) system -->\n";

	//printing document
	cout<<doc;   

	cerr<<"Finished successfully"<<endl;

	}
	catch(std::exception const &ec)
	{
		cerr<<"Error occured: "<<ec.what()<<endl;
		return EXIT_FAILURE;
	}

	return EXIT_SUCCESS;
}
