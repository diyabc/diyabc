/*
 * reftable.h
 *
 *  Created on: 9 déc. 2011
 *      Author: ppudlo
 */

#pragma once

#include <string>
#include <vector>
#include <fstream>

#include "history.hpp"
#include "header.hpp"

using namespace std;

class enregC {
public:
	int numscen;
	vector<float> param;
	vector<float> paramvv;
	vector<float> stat;
	long double dist;
	std::string message;
	friend bool operator<(const enregC& lhs, const enregC& rhs);
};


class ReftableC {
public:
	int nrec, nscen, nreclus, nrec0;
	vector<int> nrecscen;
	long posnrec;
	std::string datapath, filename, filelog, filename0, filerefscen;
	int nstat, po, nparamax, nscenchoisi, *scenchoisi, scenteste, nparamut, *nhistparam;
	vector<int> nparam;


	float *param, *sumstat;
	int histparamlength;
	std::vector<std::vector<HistParameterC>> histparam;
	MutParameterC* mutparam;
	std::fstream fifo;
	int nstatOK, nsel, nenr;
	enregC* enrsel;
	float* stat_obs;
	long double* var_stat;

	ReftableC(): scenchoisi(NULL),
	             nhistparam(NULL), param(NULL), sumstat(NULL), histparam(0),
	             mutparam(NULL), enrsel(NULL), stat_obs(NULL), var_stat(NULL) {
		histparamlength = 0;
	};

	~ReftableC() {
		if (not nrecscen.empty()) nrecscen.clear();
		if (not nparam.empty()) nparam.clear();
		if (scenchoisi != NULL) {
			delete [] scenchoisi;
			scenchoisi = NULL;
		}
		if (scenchoisi != NULL) {
			delete [] scenchoisi;
			scenchoisi = NULL;
		}
		if (nhistparam != NULL) {
			delete [] nhistparam;
			nhistparam = NULL;
		}
		if (param != NULL) {
			delete [] param;
			param = NULL;
		}
		if (sumstat != NULL) {
			delete [] sumstat;
			sumstat = NULL;
		}
		if (mutparam != NULL) {
			delete [] mutparam;
			mutparam = NULL;
		}
		if (enrsel != NULL) {
			delete [] enrsel;
			enrsel = NULL;
		}
		if (stat_obs != NULL) {
			delete [] stat_obs;
			stat_obs = NULL;
		}
		if (var_stat != NULL) {
			delete [] var_stat;
			var_stat = NULL;
		}
	};

	void sethistparamname(HeaderC const& header);
	int readheader(std::string fname, std::string flogname, std::string reftabscen);
	int writeheader();
	int readrecord(enregC* enr);
	int writerecords(int nenr, enregC* enr);
	int openfile();
	int openfile2();
	int testfile(std::string reftablefilename, int npart);
	int closefile();
	void bintotxt();
	void bintotxt2();
	void bintocsv(HeaderC const& header);
	void concat();
	void concat2();
	// calcule les variances des statistiques résumées
	// sur les 100000 premiers enregistrements de la table de référence
	int cal_varstat();
	// alloue / desalloue la mémoire pour enrsel
	void alloue_enrsel(int nsel);
	void desalloue_enrsel(int nsel);
	// calcule la distance de chaque jeu de données simulé au jeu observé
	// et sélectionne les nsel enregistrements les plus proches (copiés dans enregC *enrsel)
	void cal_dist(int nrec, int nsel, float* stat_obs, bool scenarioteste, bool allscenarios);
	int readparam(vector<float>& param);
};
