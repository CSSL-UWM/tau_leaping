#pragma once

#include <vector>
#include <map>
#include <string>

typedef std::map<double, std::vector<int> > stat_t;

typedef std::map<double, std::vector<std::vector<double> > > StepStat_t;

typedef std::map<double, std::vector<std::pair<double, double> > > StdDevStat_t;

void AllStat2StepStat(std::vector<stat_t> const &allStat, StepStat_t &stepStat, double timeStep);
void StepStat2StdDev(StepStat_t const &stepStat, StdDevStat_t &meanStdDev);
void StatToFile(std::string fileName, StdDevStat_t const &meanStdDev);