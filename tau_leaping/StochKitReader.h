#pragma once

//class for reading reaction network out of "StochKit" XML file

#include <rapidxml/rapidxml.hpp>
#include <string>
#include <vector>

#include <unordered_map>

#ifdef _WIN32
namespace std{
	using namespace tr1;
}
#endif

//auxilary structure
struct ReactionProp
{
	float rate;
	std::unordered_map<std::string, char> reactants;
	std::unordered_map<std::string, char> products;
};


class StochKitReader
{
public:
	StochKitReader(void);
	~StochKitReader(void);
	void ReadXML(const char *fn, std::vector<ReactionProp> &reactionProps,
		std::unordered_map<std::string, std::pair<int, int> > &reactantsPopulations)const;

private:
		static rapidxml::xml_node<> *FirstNode(rapidxml::xml_node<> *parentNode, const char *nodeName, bool throwOnError);
		void ReadReactionsList(rapidxml::xml_node<> *modelNode,
			std::unordered_map<std::string, float> const &rIdConstMap,
			std::vector<ReactionProp> &reactionProps) const;
		void ReadReactantsPopulations(rapidxml::xml_node<> *modelNode,
			std::unordered_map<std::string, std::pair<int, int> > &reactantsPopulations) const;

		void ReadReactionsRates(rapidxml::xml_node<> *modelNode, std::unordered_map<std::string, float>& rIdConstMap)const;
};
