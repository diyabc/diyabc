/*
 * header.h
 *
 *  Created on: 9 déc. 2011
 *      Author: ppudlo
 */

#ifndef HEADER_H_
#define HEADER_H_
#include <string>
#include <fstream>
#include <vector>
#include "data\.hpp"
#include "history\.hpp"
#include "particuleC\.hpp"
#include "header\.hpp"

class MutParameterC
{
public:
    std::string name;
    int groupe;
    int category;   //0 pour mutmoy, 1 pour Pmoy, 2 pour snimoy, 3 pour musmoy, 4 pour k1moy et 5 pour k2moy
    double value;
    PriorC prior;
//	MutParameterC(MutParameterC const & source);
	MutParameterC & operator=(MutParameterC const & source);
	
    void ecris();
};


class HeaderC
{
public:
    std::string message,datafilename,entete,entetehist,entetemut,entetemut0,entetestat;
    std::string pathbase;
    int nparamtot,nstat,nstatsnp,nscenarios,nconditions,ngroupes,nparamut,nsimfile;
    std::vector<std::string> statname;
    ScenarioC scen;
    std::vector<HistParameterC> histparam;
    std::vector<ConditionC> condition;
    bool drawuntil,reference;
    std::vector<MutParameterC> mutparam;
    std::vector<float> stat_obs;
	float reffreqmin;
	float threshold;
  
	void libere();
	void assignloc(int gr);

	int readHeaderDebut(std::ifstream & file);
	int readHeaderScenarios(std::ifstream & file);
	int readHeaderHistParam(std::ifstream & file);
	int readHeaderLoci(std::ifstream & file);
	int readHeaderGroupPrior(std::ifstream & file);
	int readHeaderAllStat(std::ifstream & file);
	int readHeaderGroupStat(std::ifstream & file);
	int readHeaderEntete(std::ifstream & file);
	int buildSuperScen();
	int buildMutParam();
	int readHeader(std::string headerfilename);
	
	int readHeadersimDebut(std::ifstream & file);
	int readHeadersimScenario(std::ifstream & file);
	int readHeadersimHistParam(std::ifstream & file);
	int readHeadersimLoci(std::ifstream & file);
	int readHeadersimGroupPrior(std::ifstream & file);
	int readHeadersimGroupSNP();
	int readHeadersimFin();
	int readHeadersim(std::string headersimfilename);

	string calstatobs(std::string statobsfilename);
 private:
	//HeaderC(const HeaderC & source) {};
	//HeaderC & operator=(const HeaderC & source) { return *this;} ;
};


#endif /* HEADER_H_ */
