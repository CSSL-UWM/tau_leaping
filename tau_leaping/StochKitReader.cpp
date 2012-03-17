#include "StochKitReader.h"
#include "functions.h"

using std::vector;
using rapidxml::xml_document;
using rapidxml::xml_node;
using std::unordered_map;
using std::string;
using std::cerr;
using std::endl;
using std::pair;
using rapidxml::xml_attribute;

StochKitReader::StochKitReader(void)
{
}

StochKitReader::~StochKitReader(void)
{
}

rapidxml::xml_node<> *StochKitReader::FirstNode(rapidxml::xml_node<> *parentNode, const char *nodeName, bool throwOnError)
{
	xml_node<> *node=parentNode->first_node(nodeName);
	if(node==NULL&&throwOnError)
	{
		char err[250];
		sprintf(err, "Model::readNode: can not find %s node", nodeName);
		throw std::runtime_error(err);
	}
	return node;
}


void StochKitReader::ReadReactionsList(rapidxml::xml_node<> *modelNode,
			std::unordered_map<std::string, float> const &rIdConstMap,
			std::vector<ReactionProp> &reactionProps) const
{
	
	xml_node<> *listParentNode=FirstNode(modelNode, "ReactionsList", true);

	xml_node<> *node=FirstNode(listParentNode, "Reaction", true);
	int i=0;
	do
	{
		++i;
		ReactionProp rp;
		const char *rName=FirstNode(node, "Rate", true)->value();
		unordered_map<string, float>::const_iterator cit=rIdConstMap.find(rName);
		if(cit==rIdConstMap.end())
		{
			cerr<<"Can not find the rate for the reaction "<<rName;
			exit(EXIT_FAILURE);
		}
		rp.rate=cit->second;

		//reading reactants
		xml_node<> *spRef=FirstNode(FirstNode(node, "Reactants", true), "SpeciesReference", false);
		if(spRef)
		{
			do
			{
				const char *sName=spRef->first_attribute("id")->value();
				char const stoich=atoi(spRef->first_attribute("stoichiometry")->value());
				rp.reactants[sName]=stoich;
			}while(spRef=spRef->next_sibling("SpeciesReference"));
		}

		//reading props
		{
			auto productNode=FirstNode(node, "Products", false);
			if(productNode)//there may be no products at all...
			{
				spRef=FirstNode(productNode, "SpeciesReference", false);
				if(spRef)
				{
					do
					{
						const char *sName=spRef->first_attribute("id")->value();
						char const stoich=atoi(spRef->first_attribute("stoichiometry")->value());
						rp.products[sName]=stoich;
					}while(spRef=spRef->next_sibling("SpeciesReference"));
				}
			}
		}


		reactionProps.push_back(rp);
	}while(node=node->next_sibling("Reaction"));
	
}
void StochKitReader::ReadReactionsRates(xml_node<> *modelNode, unordered_map<std::string, float>& rIdConstMap) const
{
	xml_node<> *listParentNode=FirstNode(modelNode, "ParametersList", true);

	xml_node<> *node=FirstNode(listParentNode, "Parameter", true);

	do
	{
		const char *name=FirstNode(node, "Id", true)->value();
		float rConst=static_cast<float>(atof(FirstNode(node, "Expression", true)->value()));
		rIdConstMap[name]=rConst;
	}while((node=node->next_sibling("Parameter")));
}

void StochKitReader::ReadReactantsPopulations(xml_node<> *modelNode, 
											  unordered_map<std::string, pair<int, int> > &reactantsPopulations) const
{

	xml_node<> *node=FirstNode(FirstNode(modelNode, "SpeciesList", true), "Species", true);

	int i=0;
	do
	{
		const char *name=FirstNode(node, "Id", true)->value();
		int count =atoi(FirstNode(node, "InitialPopulation", true)->value());
		reactantsPopulations[name]=std::make_pair(i++, count);

	}while(node=node->next_sibling("Species"));
	
}
void StochKitReader::ReadXML(const char *fn, vector<ReactionProp> &reactionProps, 
							 std::unordered_map<string, std::pair<int, int> > &reactantsPopulations)const 
{

	char err[250];
	if(!FileExists(fn))
	{
		sprintf(err, "Model::ReadXML: can not open file %s", fn);
		throw std::runtime_error(err);
	}

	vector<char> buff;
	FileToMem(fn, buff);
	
	xml_document<> doc;
	doc.parse<0>(&buff[0]);

	xml_node<> *modelNode = doc.first_node("Model");
	if(modelNode==NULL)
	{
		sprintf(err, "Model::ReadXML: can not find %s node", "Model");
		throw std::runtime_error(err);
	}

	xml_node<> *node=FirstNode(modelNode, "NumberOfReactions", true);
	const int numOfReactions=atoi(node->value());

	node=FirstNode(modelNode, "NumberOfSpecies", true);
	const int numOfSpecies=atoi(node->value());

	unordered_map<string, float> rIdConstMap;
	ReadReactionsRates(modelNode, rIdConstMap);

	/*if(rIdConstMap.size()!=numOfReactions)
	{
		cerr<<"Mismatch reaction number: rIdConstMap size is "<<rIdConstMap.size()<<"; "<<numOfReactions<<" expected\n";
		exit(EXIT_FAILURE);
	}*/
	
	
	ReadReactionsList(modelNode, rIdConstMap, reactionProps);

	rIdConstMap.clear();

	if(reactionProps.size()!=numOfReactions)
	{
		cerr<<"Mismatch reaction number 2"<<endl;
		exit(EXIT_FAILURE);
	}

	ReadReactantsPopulations(modelNode, reactantsPopulations);
	
	if(reactantsPopulations.size()!=numOfSpecies)
	{
		cerr<<"Mismatch reaction number 3"<<endl;
		exit(EXIT_FAILURE);
	}
}
