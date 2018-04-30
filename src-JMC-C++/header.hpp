/*
 * header.h
 *
 *  Created on: 9 déc. 2011
 *      Author: ppudlo
 */

#pragma once


#include <string>
#include <vector>
#include <iostream>

#include "history.hpp"

using namespace std;

class MutParameterC {
public:
	string name;
	int groupe;
	int category; //0 pour mutmoy, 1 pour Pmoy, 2 pour snimoy, 3 pour musmoy, 4 pour k1moy et 5 pour k2moy
	double value;
	PriorC prior;
	//	MutParameterC(MutParameterC const & source);
	MutParameterC& operator=(MutParameterC const& source);

	void ecris();
};


class HeaderC {
public:
	string message, datafilename, entete, entetehist, entetemut, entetemut0, entetestat;
	string pathbase;
	int nparamtot, nstat, nstatsnp, nscenarios, nconditions, ngroupes, nparamut, nsimfile;
	vector<string> statname;
	ScenarioC scen;
	vector<HistParameterC> histparam;
	vector<ConditionC> condition;
	bool drawuntil, reference;
	vector<MutParameterC> mutparam;
	vector<float> stat_obs;
	float reffreqmin;
	float threshold;

	void libere();
	void assignloc(bool sim, int gr);

	int readHeaderDebut(ifstream& file);
	int readHeaderScenarios(ifstream& file);
	int readHeaderHistParam(ifstream& file);
	int readHeaderLoci(ifstream& file);
	int readHeaderGroupPrior(ifstream& file);
	int readHeaderAllStat(ifstream& file, string headerfilename);
	int readHeaderGroupStat(ifstream& file);
	int readHeaderLine(string& ss, bool statname);
	int readHeaderEntete(ifstream& file);
	int buildSuperScen();
	int buildMutParam();
	int readHeader(string headerfilename);

	int readHeadersimDebut(ifstream& file);
	int readHeadersimScenario(ifstream& file);
	int readHeadersimHistParam(ifstream& file);
	int readHeadersimLoci(ifstream& file);
	int readHeadersimGroupPrior(ifstream& file);
	int readHeadersimGroupSNP();
	int readHeadersimFin();
	int readHeadersim(string headersimfilename);

	string calstatobs(string statobsfilename);
private:
	//HeaderC(const HeaderC & source) {};
	//HeaderC & operator=(const HeaderC & source) { return *this;} ;
};
