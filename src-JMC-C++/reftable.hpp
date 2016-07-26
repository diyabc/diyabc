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
	vector<float> stat;
	long double dist;
	std::string message;
	friend bool operator<(const enregC& lhs, const enregC& rhs) {
		return lhs.dist < rhs.dist;
	}
};


class ReftableC {
public:
	int nrec, nscen, nreclus, nrec0;
	vector<int> nrecscen;
	std::string datapath, filename, filelog, filename0, filerefscen;
	int nstat, nparamax, nscenchoisi, scenteste, nparamut;
	vector<int> nhistparam;
	vector<int> scenchoisi;
	vector<int> nparam;

	int histparamlength = 0;
	std::vector<std::vector<HistParameterC>> histparam;
	vector<MutParameterC> mutparam;
	std::fstream fifo;
	vector<enregC> enrsel;
	vector<long double> var_stat;

	void sethistparamname(HeaderC const& header);
	int readheader(std::string fname, std::string flogname, std::string reftabscen);
	int writeheader();
	int readrecord(enregC& enr);
	int writerecords(int nenr, vector<enregC>& enr);
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
	void desalloue_enrsel();
	// calcule la distance de chaque jeu de données simulé au jeu observé
	// et sélectionne les nsel enregistrements les plus proches (copiés dans enregC *enrsel)
	void cal_dist(int nrec, int nsel, float* stat_obs, bool scenarioteste, bool allscenarios);
	int readparam(vector<float>& param);
};
