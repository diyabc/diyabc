/*
 * bias.cpp
 *
 *  Created on: 23 march 2011
 *      Author: cornuet
 */

#include <string>
#include <vector>
#include <cstdlib>
#include <cmath>
#include <time.h>
#include <iomanip>


#include "bias.hpp"
#include "estimparam.hpp"
#include "reftable.hpp"
#include "header.hpp"
#include "mesutils.hpp"
#include "particleset.hpp"
#include "modchec.hpp"


/*
#ifndef HEADER
#include "header.cpp"
#define HEADER
#endif
#ifndef MATRICES
#include "matrices.cpp"
#define MATRICES
#endif
#ifndef MESUTILS
#include "mesutils.cpp"
#define MESUTILS
#endif
#ifndef PARTICLESET
#include "particleset.cpp"
#define PARTICLESET
#endif
*/
using namespace std;

extern ReftableC rt;
extern ParticleSetC ps;
extern enregC* enreg;
extern bool multithread;
extern string progressfilename, path;
extern string scurfile;
extern HeaderC header;
extern int nparamcom, nparcompo, nparscaled, **numpar, numtransf, npar, nstatOKsel, nenr;
extern string ident, headerfilename;
extern bool original, composite, scaled;
extern string *nomparamO, *nomparamC, *nomparamS;
//extern long double **phistar,**phistarcompo,**phistarscaled;
extern int debuglevel;
extern ofstream fprog, fpar;
extern long double **simpar, **simparcompo, **simparscaled;
extern int nsimpar;
extern long double** phistarOK;
extern vector<ScenarioC> scenario;
extern vector<LocusGroupC> groupe;

extern bool deltanul;

parstatC **paramest, **paramestcompo, **paramestscaled;
parstatC *paramestS, *paramestcompoS, *paramestscaledS;

long double **simparsel, **simparcomposel, **simparscaledsel, *paramestmoyS, *paramestmedS, *paramestcompomoyS, *paramestcompomedS, *paramestscaledmoyS, *paramestscaledmedS;
/*
 * paramestS[numero du parametre] contient les statistiques évaluées à partir des  10⁵ premiers enregistrements de la table de reference :
 * moyenne, médiane, mode, quantiles de la distribution a priori de chaque paramètre tiré selon son prior 
 * avec tirage conforme aux éventuelles conditions
*/
enreC* enreg2;
long double **br_O, *rrmise_O, *rmad_O, **rmse_O, *cov50_O, *cov95_O, **fac2_O, **rmb_O, *rmedad_O, **rmae_O, ***bmed_O, **bmedr_O, ***bmeda_O;
long double **br_C, *rrmise_C, *rmad_C, **rmse_C, *cov50_C, *cov95_C, **fac2_C, **rmb_C, *rmedad_C, **rmae_C, ***bmed_C, **bmedr_C, ***bmeda_C;
long double **br_S, *rrmise_S, *rmad_S, **rmse_S, *cov50_S, *cov95_S, **fac2_S, **rmb_S, *rmedad_S, **rmae_S, ***bmed_S, **bmedr_S, ***bmeda_S;
long double **br_OS, *rrmise_OS, *rmad_OS, **rmse_OS, *cov50_OS, *cov95_OS, **fac2_OS, **rmb_OS, *rmedad_OS, **rmae_OS, ***bmed_OS, **bmedr_OS, ***bmeda_OS;
long double **br_CS, *rrmise_CS, *rmad_CS, **rmse_CS, *cov50_CS, *cov95_CS, **fac2_CS, **rmb_CS, *rmedad_CS, **rmae_CS, ***bmed_CS, **bmedr_CS, ***bmeda_CS;
long double **br_SS, *rrmise_SS, *rmad_SS, **rmse_SS, *cov50_SS, *cov95_SS, **fac2_SS, **rmb_SS, *rmedad_SS, **rmae_SS, ***bmed_SS, **bmedr_SS, ***bmeda_SS;
long double **paretoil, **paretoilcompo, **paretoilscaled;
long double **paretoilS, **paretoilcompoS, **paretoilscaledS;
ofstream ftrace;
MwcGen mw;
int* numero;
double zeroplus = 1E-100;
bool valinfinie = false;


/**
 * tire un échantillon avec remise dans les simpar
 * 
 */

long double** echantillon(int nsel, int nparam, long double** simpar) {
	long double** sp;
	sp = new long double*[nsel];
	for (int i = 0; i < nsel; i++) {
		sp[i] = new long double[nparam];
		for (int j = 0; j < nparam; j++) sp[i][j] = simpar[numero[i]][j];
	}
	return sp;
}

/*
    string pseudoprior(string s) {
        string spr;
        long double x=atof(s.c_str());
        long double mini=0.99999*x,maxi=1.00001*x;
        spr="UN["+DoubleToString(mini)+","+DoubleToString(maxi)+",0.0,0.0";
        return spr;
    }
*/

/**
 *  modifie les paramètres historiques en accord avec les priors/valeurs des pseudo-observés
 */
void resethistparam(string s) {
	string name, sprior, smini, smaxi;
	vector<string> ss;
	splitwords(s, "=", ss);
	//int n = ss.size();
	name = ss[0];
	int i = 0;
	while ((i < scenario[rt.scenteste - 1].nparam)and (name != scenario[rt.scenteste - 1].histparam[i].name)) i++;
	//cout<<"resethistparam   parametre "<<name<<"  ("<<i<<")   "<<ss[1]<<"\n";
	if (ss[1].find("[") != string::npos) {
		//cout<<"resethistparam   avant readprior\n ";
		scenario[rt.scenteste - 1].histparam[i].prior.readprior(ss[1]);
		//cout<<"resethistparam   apres readprior\n ";
		scenario[rt.scenteste - 1].histparam[i].prior.fixed = false;
		//cout<<"resethistparam   avant ecris\n ";
		//scenario[rt.scenteste-1].histparam[i].prior.ecris();
	}
	else {
		scenario[rt.scenteste - 1].histparam[i].value = atof(ss[1].c_str());
		scenario[rt.scenteste - 1].histparam[i].prior.fixed = true;
		//cout<<scenario[rt.scenteste-1].histparam[i].value<<"\n";
	}
}

/**
 *  modifie les conditions sur les paramètres historiques en accord avec les priors des pseudo-observés
 */
void resetcondition(int j, string s) {
	scenario[rt.scenteste - 1].condition[j].readcondition(s);

}

/**
 *  modifie les paramètres mutationnels en accord avec les priors des pseudo-observés
 */
void resetmutparam(string s) {
	cout << "debut de resetmutparam\n";
	string numgr, s1, sg;
	vector<string> ss;
	int gr, i0, i1;
	numgr = s.substr(1, s.find("(") - 1);
	gr = atoi(numgr.c_str());
	i0 = s.find("(");
	i1 = s.find(")");
	cout << "i0=" << i0 << "  i1=" << i1 << "\n";
	s1 = s.substr(i0 + 1, i1 - i0 - 1);
	cout << "groupe " << gr << "  " << s1 << "\n";
	splitwords(s1, " ", ss);
	//int n =(int)ss.size();
	if (groupe[gr].type == 0) {
		//cout<<"mutmoy : \n";
		//groupe[gr].priormutmoy.ecris();
		if (groupe[gr].priormutmoy.constant) groupe[gr].mutmoy = groupe[gr].priormutmoy.mini;
		else {
			if (ss[0].find("[") == string::npos) {
				groupe[gr].mutmoy = atof(ss[0].c_str());
				groupe[gr].priormutmoy.fixed = true;
				cout << "mutmoy=" << groupe[gr].mutmoy << "\n";
			}
			else {
				groupe[gr].priormutmoy.readprior(ss[0]);
				groupe[gr].priormutmoy.fixed = false;

			}
		}
		//cout<<"mutmoy : \n";
		//groupe[gr].priormutmoy.ecris();
		if (ss[1].find("[") == string::npos) groupe[gr].priormutloc.sdshape = atof(ss[1].c_str());
		else groupe[gr].priormutloc.readprior(ss[1]);

		//cout<<"Pmoy : \n";
		//groupe[gr].priorPmoy.ecris();
		if (groupe[gr].priorPmoy.constant) groupe[gr].Pmoy = groupe[gr].priorPmoy.mini;
		else {
			if (ss[2].find("[") == string::npos) {
				groupe[gr].Pmoy = atof(ss[2].c_str());
				groupe[gr].priorPmoy.fixed = true;
			}
			else {
				groupe[gr].priorPmoy.readprior(ss[2]);
				groupe[gr].priorPmoy.fixed = false;
			}
		}
		//cout<<"Pmoy : \n";
		//groupe[gr].priorPmoy.ecris();
		if (ss[3].find("[") == string::npos) groupe[gr].priorPloc.sdshape = atof(ss[3].c_str());
		else groupe[gr].priorPloc.readprior(ss[3]);

		//cout<<"snimoy : \n";
		//groupe[gr].priorsnimoy.ecris();
		if (groupe[gr].priorsnimoy.constant) groupe[gr].snimoy = groupe[gr].priorsnimoy.mini;
		else {
			if (ss[4].find("[") == string::npos) {
				groupe[gr].snimoy = atof(ss[4].c_str());
				groupe[gr].priorsnimoy.fixed = true;
			}
			else {
				groupe[gr].priorsnimoy.readprior(ss[4]);
				groupe[gr].priorsnimoy.fixed = false;
			}
		}
		//cout<<"snimoy : \n";
		//groupe[gr].priorsnimoy.ecris();
		if (ss[5].find("[") == string::npos) groupe[gr].priorsniloc.sdshape = atof(ss[5].c_str());
		else groupe[gr].priorsniloc.readprior(ss[5]);
	}
	else if (groupe[gr].type == 1) {
		//cout<<"resetmutparam type sequence\n";
		//cout<<"modele "<< groupe[gr].mutmod<<"\n";
		if (groupe[gr].priormusmoy.constant) groupe[gr].musmoy = groupe[gr].priormusmoy.mini;
		else {
			if (ss[0].find("[") == string::npos) {
				groupe[gr].musmoy = atof(ss[0].c_str());
				groupe[gr].priormusmoy.fixed = true;
			}
			else {
				groupe[gr].priormusmoy.readprior(ss[0]);
				groupe[gr].priormusmoy.fixed = false;
			}
		}
		//cout<<"ss[0] : >"<<ss[0]<<"<      musmoy="<<groupe[gr].musmoy<<"\n";
		if (ss[1].find("[") == string::npos) groupe[gr].priormusloc.sdshape = atof(ss[1].c_str());
		else groupe[gr].priormusloc.readprior(ss[1]);
		//cout<<"musloc"<<"\n";
		if (groupe[gr].mutmod > 0) {
			if (groupe[gr].priork1moy.constant) groupe[gr].k1moy = groupe[gr].priork1moy.mini;
			else {
				if (ss[2].find("[") == string::npos) {
					groupe[gr].k1moy = atof(ss[2].c_str());
					groupe[gr].priork1moy.fixed = true;
				}
				else {
					groupe[gr].priork1moy.readprior(ss[2]);
					groupe[gr].priork1moy.fixed = false;
				}
			}
			//cout<<"k1moy="<<groupe[gr].k1moy<<"\n";
			if (ss[3].find("[") == string::npos) groupe[gr].priork1loc.sdshape = atof(ss[3].c_str());
			else groupe[gr].priork1loc.readprior(ss[3]);
			//cout<<"k1loc\n";
			if (groupe[gr].mutmod > 2) {
				if (groupe[gr].priork1moy.constant) groupe[gr].k1moy = groupe[gr].priork1moy.mini;
				else {
					if (ss[4].find("[") == string::npos) {
						groupe[gr].k2moy = atof(ss[4].c_str());
						groupe[gr].priork2moy.fixed = true;
					}
					else {
						groupe[gr].priork2moy.readprior(ss[4]);
						groupe[gr].priork2moy.fixed = false;
					}
				}
				//cout<<"k2moy="<<groupe[gr].k2moy<<"\n";
				if (ss[5].find("[") == string::npos) groupe[gr].priork2loc.sdshape = atof(ss[5].c_str());
				else groupe[gr].priork2loc.readprior(ss[5]);
			}
			//cout<<"fin\n";
		}
	}
}


/**
* initialise les tableaux et les différentes statistiques de biais, rmse... pour les paramètres originaux
*/
void initbiasO(int ntest, int nsel) {
	paramest = new parstatC*[ntest];
	paretoil = new long double*[nsel];
	for (int i = 0; i < nsel; i++)paretoil[i] = new long double[nparamcom];
	//////////////// mean relative bias
	br_O = new long double* [3];
	for (int k = 0; k < 3; k++) {
		br_O[k] = new long double[nparamcom];
		for (int j = 0; j < nparamcom; j++) br_O[k][j] = 0.0;
	}
	////////////  RRMISE
	rrmise_O = new long double [nparamcom];
	for (int j = 0; j < nparamcom; j++) rrmise_O[j] = 0.0;
	////////////  RMAD
	rmad_O = new long double [nparamcom];
	for (int j = 0; j < nparamcom; j++) rmad_O[j] = 0.0;
	////////////  RMSE
	rmse_O = new long double*[3];
	for (int k = 0; k < 3; k++) {
		rmse_O[k] = new long double[nparamcom];
		for (int j = 0; j < nparamcom; j++) rmse_O[k][j] = 0.0;
	}
	/////////////// coverages
	cov50_O = new long double[nparamcom];
	cov95_O = new long double[nparamcom];
	for (int j = 0; j < nparamcom; j++) {
		cov50_O[j] = 0.0;
		cov95_O[j] = 0.0;
	}
	///////////////// factors 2
	fac2_O = new long double*[3];
	for (int k = 0; k < 3; k++) {
		fac2_O[k] = new long double[nparamcom];
		for (int j = 0; j < nparamcom; j++) fac2_O[k][j] = 0.0;
	}
	/////////////////// medianes du biais relatif
	bmed_O = new long double**[3];
	for (int k = 0; k < 3; k++) {
		bmed_O[k] = new long double*[nparamcom];
		for (int j = 0; j < nparamcom; j++) {
			bmed_O[k][j] = new long double[ntest];
		}
	}
	rmb_O = new long double*[3];
	for (int k = 0; k < 3; k++) rmb_O[k] = new long double[nparamcom];
	///////////////// Relative Median Absolute Deviation
	bmedr_O = new long double*[nparamcom];
	for (int j = 0; j < nparamcom; j++) bmedr_O[j] = new long double[ntest];
	rmedad_O = new long double [nparamcom];
	////////////  RMAE
	rmae_O = new long double*[3];
	for (int k = 0; k < 3; k++) rmae_O[k] = new long double[nparamcom];
	bmeda_O = new long double**[3];
	for (int k = 0; k < 3; k++) {
		bmeda_O[k] = new long double*[nparamcom];
		for (int j = 0; j < nparamcom; j++) bmeda_O[k][j] = new long double[ntest];
	}
}

void initbiasOS(int ntest, int nsel) {
	paretoilS = new long double*[nsel];
	for (int i = 0; i < nsel; i++)paretoilS[i] = new long double[nparamcom];
	paramestmoyS = new long double[nparamcom];
	for (int j = 0; j < nparamcom; j++) paramestmoyS[j] = 0.0;
	paramestmedS = new long double[nparamcom];
	for (int j = 0; j < nparamcom; j++) paramestmedS[j] = 0.0;
	//////////////// mean relative bias
	br_OS = new long double* [3];
	for (int k = 0; k < 3; k++) {
		br_OS[k] = new long double[nparamcom];
		for (int j = 0; j < nparamcom; j++) br_OS[k][j] = 0.0;
	}
	////////////  RRMISE
	rrmise_OS = new long double [nparamcom];
	for (int j = 0; j < nparamcom; j++) rrmise_OS[j] = 0.0;
	////////////  RMAD
	rmad_OS = new long double [nparamcom];
	for (int j = 0; j < nparamcom; j++) rmad_OS[j] = 0.0;
	////////////  RMSE
	rmse_OS = new long double*[3];
	for (int k = 0; k < 3; k++) {
		rmse_OS[k] = new long double[nparamcom];
		for (int j = 0; j < nparamcom; j++) rmse_OS[k][j] = 0.0;
	}
	/////////////// coverages
	cov50_OS = new long double[nparamcom];
	cov95_OS = new long double[nparamcom];
	for (int j = 0; j < nparamcom; j++) {
		cov50_OS[j] = 0.0;
		cov95_OS[j] = 0.0;
	}
	///////////////// factors 2
	fac2_OS = new long double*[3];
	for (int k = 0; k < 3; k++) {
		fac2_OS[k] = new long double[nparamcom];
		for (int j = 0; j < nparamcom; j++) fac2_OS[k][j] = 0.0;
	}
	/////////////////// medianes du biais relatif
	bmed_OS = new long double**[3];
	for (int k = 0; k < 3; k++) {
		bmed_OS[k] = new long double*[nparamcom];
		for (int j = 0; j < nparamcom; j++) {
			bmed_OS[k][j] = new long double[ntest];
		}
	}
	rmb_OS = new long double*[3];
	for (int k = 0; k < 3; k++) rmb_OS[k] = new long double[nparamcom];
	///////////////// Relative Median Absolute Deviation
	bmedr_OS = new long double*[nparamcom];
	for (int j = 0; j < nparamcom; j++) bmedr_OS[j] = new long double[ntest];
	rmedad_OS = new long double [nparamcom];
	////////////  RMAE
	rmae_OS = new long double*[3];
	for (int k = 0; k < 3; k++) rmae_OS[k] = new long double[nparamcom];
	bmeda_OS = new long double**[3];
	for (int k = 0; k < 3; k++) {
		bmeda_OS[k] = new long double*[nparamcom];
		for (int j = 0; j < nparamcom; j++) bmeda_OS[k][j] = new long double[ntest];
	}
}

/**
* initialise les tableaux et les différentes statistiques de biais, rmse... pour les paramètres composites
*/
void initbiasC(int ntest, int nsel) {
	paramestcompo = new parstatC*[ntest];
	paretoilcompo = new long double*[nsel];
	for (int i = 0; i < nsel; i++)paretoilcompo[i] = new long double[nparcompo];
	//////////////// mean relative bias
	br_C = new long double* [3];
	for (int k = 0; k < 3; k++) {
		br_C[k] = new long double[nparcompo];
		for (int j = 0; j < nparcompo; j++) br_C[k][j] = 0.0;
	}
	////////////  RRMISE
	rrmise_C = new long double [nparcompo];
	for (int j = 0; j < nparcompo; j++) rrmise_C[j] = 0.0;
	////////////  RMAD
	rmad_C = new long double [nparcompo];
	for (int j = 0; j < nparcompo; j++) rmad_C[j] = 0.0;
	////////////  RMSE
	rmse_C = new long double*[3];
	for (int k = 0; k < 3; k++) {
		rmse_C[k] = new long double[nparcompo];
		for (int j = 0; j < nparcompo; j++) rmse_C[k][j] = 0.0;
	}
	/////////////// coverages
	cov50_C = new long double[nparcompo];
	cov95_C = new long double[nparcompo];
	for (int j = 0; j < nparcompo; j++) {
		cov50_C[j] = 0.0;
		cov95_C[j] = 0.0;
	}
	///////////////// factors 2
	fac2_C = new long double*[3];
	for (int k = 0; k < 3; k++) {
		fac2_C[k] = new long double[nparcompo];
		for (int j = 0; j < nparcompo; j++) fac2_C[k][j] = 0.0;
	}
	/////////////////// medianes du biais relatif
	bmed_C = new long double**[3];
	for (int k = 0; k < 3; k++) {
		bmed_C[k] = new long double*[nparcompo];
		for (int j = 0; j < nparcompo; j++) {
			bmed_C[k][j] = new long double[ntest];
		}
	}
	rmb_C = new long double*[3];
	for (int k = 0; k < 3; k++) rmb_C[k] = new long double[nparcompo];
	///////////////// Relative Median Absolute Deviation
	bmedr_C = new long double*[nparcompo];
	for (int j = 0; j < nparcompo; j++) bmedr_C[j] = new long double[ntest];
	rmedad_C = new long double [nparcompo];
	////////////  RMAE
	rmae_C = new long double*[3];
	for (int k = 0; k < 3; k++) rmae_C[k] = new long double[nparcompo];
	bmeda_C = new long double**[3];
	for (int k = 0; k < 3; k++) {
		bmeda_C[k] = new long double*[nparcompo];
		for (int j = 0; j < nparcompo; j++) bmeda_C[k][j] = new long double[ntest];
	}
}

void initbiasCS(int ntest, int nsel) {
	paretoilcompoS = new long double*[nsel];
	for (int i = 0; i < nsel; i++)paretoilcompoS[i] = new long double[nparcompo];
	paramestcompomoyS = new long double[nparamcom];
	for (int j = 0; j < nparamcom; j++) paramestcompomoyS[j] = 0.0;
	paramestcompomedS = new long double[nparamcom];
	for (int j = 0; j < nparamcom; j++) paramestcompomedS[j] = 0.0;
	//////////////// mean relative bias
	br_CS = new long double* [3];
	for (int k = 0; k < 3; k++) {
		br_CS[k] = new long double[nparcompo];
		for (int j = 0; j < nparcompo; j++) br_CS[k][j] = 0.0;
	}
	////////////  RRMISE
	rrmise_CS = new long double [nparcompo];
	for (int j = 0; j < nparcompo; j++) rrmise_CS[j] = 0.0;
	////////////  RMAD
	rmad_CS = new long double [nparcompo];
	for (int j = 0; j < nparcompo; j++) rmad_CS[j] = 0.0;
	////////////  RMSE
	rmse_CS = new long double*[3];
	for (int k = 0; k < 3; k++) {
		rmse_CS[k] = new long double[nparcompo];
		for (int j = 0; j < nparcompo; j++) rmse_CS[k][j] = 0.0;
	}
	/////////////// coverages
	cov50_CS = new long double[nparcompo];
	cov95_CS = new long double[nparcompo];
	for (int j = 0; j < nparcompo; j++) {
		cov50_CS[j] = 0.0;
		cov95_CS[j] = 0.0;
	}
	///////////////// factors 2
	fac2_CS = new long double*[3];
	for (int k = 0; k < 3; k++) {
		fac2_CS[k] = new long double[nparcompo];
		for (int j = 0; j < nparcompo; j++) fac2_CS[k][j] = 0.0;
	}
	/////////////////// medianes du biais relatif
	bmed_CS = new long double**[3];
	for (int k = 0; k < 3; k++) {
		bmed_CS[k] = new long double*[nparcompo];
		for (int j = 0; j < nparcompo; j++) {
			bmed_CS[k][j] = new long double[ntest];
		}
	}
	rmb_CS = new long double*[3];
	for (int k = 0; k < 3; k++) rmb_CS[k] = new long double[nparcompo];
	///////////////// Relative Median Absolute Deviation
	bmedr_CS = new long double*[nparcompo];
	for (int j = 0; j < nparcompo; j++) bmedr_CS[j] = new long double[ntest];
	rmedad_CS = new long double [nparcompo];
	////////////  RMAE
	rmae_CS = new long double*[3];
	for (int k = 0; k < 3; k++) rmae_CS[k] = new long double[nparcompo];
	bmeda_CS = new long double**[3];
	for (int k = 0; k < 3; k++) {
		bmeda_CS[k] = new long double*[nparcompo];
		for (int j = 0; j < nparcompo; j++) bmeda_CS[k][j] = new long double[ntest];
	}
}


void initbiasS(int ntest, int nsel) {
	paramestscaled = new parstatC*[ntest];
	paretoilscaled = new long double*[nsel];
	for (int i = 0; i < nsel; i++)paretoilscaled[i] = new long double[nparscaled];
	//////////////// mean relative bias
	br_S = new long double* [3];
	for (int k = 0; k < 3; k++) {
		br_S[k] = new long double[nparscaled];
		for (int j = 0; j < nparscaled; j++) br_S[k][j] = 0.0;
	}
	////////////  RRMISE
	rrmise_S = new long double [nparscaled];
	for (int j = 0; j < nparscaled; j++) rrmise_S[j] = 0.0;
	////////////  RMAD
	rmad_S = new long double [nparscaled];
	for (int j = 0; j < nparscaled; j++) rmad_S[j] = 0.0;
	////////////  RMSE
	rmse_S = new long double*[3];
	for (int k = 0; k < 3; k++) {
		rmse_S[k] = new long double[nparscaled];
		for (int j = 0; j < nparscaled; j++) rmse_S[k][j] = 0.0;
	}
	/////////////// coverages
	cov50_S = new long double[nparscaled];
	cov95_S = new long double[nparscaled];
	for (int j = 0; j < nparscaled; j++) {
		cov50_S[j] = 0.0;
		cov95_S[j] = 0.0;
	}
	///////////////// factors 2
	fac2_S = new long double*[3];
	for (int k = 0; k < 3; k++) {
		fac2_S[k] = new long double[nparscaled];
		for (int j = 0; j < nparscaled; j++) fac2_S[k][j] = 0.0;
	}
	/////////////////// medianes du biais relatif
	bmed_S = new long double**[3];
	for (int k = 0; k < 3; k++) {
		bmed_S[k] = new long double*[nparscaled];
		for (int j = 0; j < nparscaled; j++) {
			bmed_S[k][j] = new long double[ntest];
		}
	}
	rmb_S = new long double*[3];
	for (int k = 0; k < 3; k++) rmb_S[k] = new long double[nparscaled];
	///////////////// Relative Median Absolute Deviation
	bmedr_S = new long double*[nparscaled];
	for (int j = 0; j < nparscaled; j++) bmedr_S[j] = new long double[ntest];
	rmedad_S = new long double [nparscaled];
	////////////  RMAE
	rmae_S = new long double*[3];
	for (int k = 0; k < 3; k++) rmae_S[k] = new long double[nparscaled];
	bmeda_S = new long double**[3];
	for (int k = 0; k < 3; k++) {
		bmeda_S[k] = new long double*[nparscaled];
		for (int j = 0; j < nparscaled; j++) bmeda_S[k][j] = new long double[ntest];
	}
}

void initbiasSS(int ntest, int nsel) {
	paretoilscaledS = new long double*[nsel];
	for (int i = 0; i < nsel; i++)paretoilscaledS[i] = new long double[nparscaled];
	paramestscaledmoyS = new long double[nparamcom];
	for (int j = 0; j < nparamcom; j++) paramestscaledmoyS[j] = 0.0;
	paramestscaledmedS = new long double[nparamcom];
	for (int j = 0; j < nparamcom; j++) paramestscaledmedS[j] = 0.0;
	//////////////// mean relative bias
	br_SS = new long double* [3];
	for (int k = 0; k < 3; k++) {
		br_SS[k] = new long double[nparscaled];
		for (int j = 0; j < nparscaled; j++) br_SS[k][j] = 0.0;
	}
	////////////  RRMISE
	rrmise_SS = new long double [nparscaled];
	for (int j = 0; j < nparscaled; j++) rrmise_SS[j] = 0.0;
	////////////  RMAD
	rmad_SS = new long double [nparscaled];
	for (int j = 0; j < nparscaled; j++) rmad_SS[j] = 0.0;
	////////////  RMSE
	rmse_SS = new long double*[3];
	for (int k = 0; k < 3; k++) {
		rmse_SS[k] = new long double[nparscaled];
		for (int j = 0; j < nparscaled; j++) rmse_SS[k][j] = 0.0;
	}
	/////////////// coverages
	cov50_SS = new long double[nparscaled];
	cov95_SS = new long double[nparscaled];
	for (int j = 0; j < nparscaled; j++) {
		cov50_SS[j] = 0.0;
		cov95_SS[j] = 0.0;
	}
	///////////////// factors 2
	fac2_SS = new long double*[3];
	for (int k = 0; k < 3; k++) {
		fac2_SS[k] = new long double[nparscaled];
		for (int j = 0; j < nparscaled; j++) fac2_SS[k][j] = 0.0;
	}
	/////////////////// medianes du biais relatif
	bmed_SS = new long double**[3];
	for (int k = 0; k < 3; k++) {
		bmed_SS[k] = new long double*[nparscaled];
		for (int j = 0; j < nparscaled; j++) {
			bmed_SS[k][j] = new long double[ntest];
		}
	}
	rmb_SS = new long double*[3];
	for (int k = 0; k < 3; k++) rmb_SS[k] = new long double[nparscaled];
	///////////////// Relative Median Absolute Deviation
	bmedr_SS = new long double*[nparscaled];
	for (int j = 0; j < nparscaled; j++) bmedr_SS[j] = new long double[ntest];
	rmedad_SS = new long double [nparscaled];
	////////////  RMAE
	rmae_SS = new long double*[3];
	for (int k = 0; k < 3; k++) rmae_SS[k] = new long double[nparscaled];
	bmeda_SS = new long double**[3];
	for (int k = 0; k < 3; k++) {
		bmeda_SS[k] = new long double*[nparscaled];
		for (int j = 0; j < nparscaled; j++) bmeda_SS[k][j] = new long double[ntest];
	}
}

/**
* ecrit les paramvv,les paramest et les sumstatdans un fichier
*/
void tracebiais(int p) {
	ftrace.precision(5);
	for (int j = 0; j < nparamcom; j++) ftrace << enreg2[p].paramvv[j] << '\t';
	for (int j = 0; j < nparamcom; j++) ftrace << paramest[p][j].moy << '\t' << paramest[p][j].med << '\t' << paramest[p][j].mod << '\t';
	for (int j = 0; j < header.nstat; j++) ftrace << enreg2[p].stat[j] << '\t';
	ftrace << '\n';
}


/**
* calcule les différentes statistiques de biais, rmse...
*/
void biaisrelO(int ntest, int nsel, int p) {
	//cout<<"debut biaisrelO  p="<<p<<"\n";
	long double s, d, ds;
	simparsel = echantillon(nsel, nparamcom, simpar);
	if (p > 0) delete []paramestS;
	paramestS = calparstat(nsel, nparamcom, simparsel);
	for (int j = 0; j < nparamcom; j++) {
		paramestmoyS[j] += paramestS[j].moy;
		paramestmedS[j] += paramestS[j].med;
	}
	//////////////// mean relative bias
	for (int k = 0; k < 3; k++) {
		for (int j = 0; j < nparamcom; j++) {
			switch (k) {
				case 0: d = paramest[p][j].moy - enreg2[p].paramvv[j];
					ds = paramestS[j].moy - enreg2[p].paramvv[j];
					break;
				case 1: d = paramest[p][j].med - enreg2[p].paramvv[j];
					ds = paramestS[j].med - enreg2[p].paramvv[j];
					break;
				case 2: d = paramest[p][j].mod - enreg2[p].paramvv[j];
					ds = paramestS[j].mod - enreg2[p].paramvv[j];
					break;
			}
			br_O[k][j] += d / enreg2[p].paramvv[j]; // il restera à diviser par ntest
			br_OS[k][j] += ds / enreg2[p].paramvv[j]; // il restera à diviser par ntest
		}
		//cout <<br_OS[k][0]<<"    ";
	}
	//cout<<"\n";
	//cout<<("1\n");
	////////////  RRMISE
	for (int j = 0; j < nparamcom; j++) {
		s = 0.0;
		for (int i = 0; i < nsel; i++) {
			d = paretoil[i][j] - enreg2[p].paramvv[j];
			s += d * d;
		}
		rrmise_O[j] += s / enreg2[p].paramvv[j] / enreg2[p].paramvv[j] / double(nsel);// il restera à diviser par ntest
		s = 0.0;
		for (int i = 0; i < nsel; i++) {
			d = simparsel[i][j] - enreg2[p].paramvv[j];
			s += d * d;
		}
		rrmise_OS[j] += s / enreg2[p].paramvv[j] / enreg2[p].paramvv[j] / double(nsel);// il restera à diviser par ntest
	}
	//cout<<("2\n");
	////////////  RMAD
	for (int j = 0; j < nparamcom; j++) {
		s = 0.0;
		for (int i = 0; i < nsel; i++) s += fabs(paretoil[i][j] - enreg2[p].paramvv[j]);
		rmad_O[j] += s / fabs(enreg2[p].paramvv[j]) / (long double)nsel;// il restera à diviser par ntest
		s = 0.0;
		for (int i = 0; i < nsel; i++) s += fabs(simparsel[i][j] - enreg2[p].paramvv[j]);
		rmad_OS[j] += s / fabs(enreg2[p].paramvv[j]) / (long double)nsel;// il restera à diviser par ntest
	}
	//cout<<("3\n");
	////////////  RMSE
	for (int k = 0; k < 3; k++) {
		for (int j = 0; j < nparamcom; j++) {
			switch (k) {
				case 0: d = paramest[p][j].moy - enreg2[p].paramvv[j];
					ds = paramestS[j].moy - enreg2[p].paramvv[j];
					break;
				case 1: d = paramest[p][j].med - enreg2[p].paramvv[j];
					ds = paramestS[j].med - enreg2[p].paramvv[j];
					break;
				case 2: d = paramest[p][j].mod - enreg2[p].paramvv[j];
					ds = paramestS[j].mod - enreg2[p].paramvv[j];
					break;
			}
			rmse_O[k][j] += d * d / enreg2[p].paramvv[j] / enreg2[p].paramvv[j];
			// il restera à prendre la racine carrée et à diviser par ntest
			rmse_OS[k][j] += ds * ds / enreg2[p].paramvv[j] / enreg2[p].paramvv[j];
		}
	}
	//cout<<("4\n");
	/////////////// coverages
	long double atest = 1.0 / (long double)ntest;
	for (int j = 0; j < nparamcom; j++) {
		if ((paramest[p][j].q025 <= enreg2[p].paramvv[j])and (paramest[p][j].q975 >= enreg2[p].paramvv[j])) cov95_O[j] += atest;
		if ((paramest[p][j].q250 <= enreg2[p].paramvv[j])and (paramest[p][j].q750 >= enreg2[p].paramvv[j])) cov50_O[j] += atest;
		if ((paramestS[j].q025 <= enreg2[p].paramvv[j])and (paramestS[j].q975 >= enreg2[p].paramvv[j])) cov95_OS[j] += atest;
		if ((paramestS[j].q250 <= enreg2[p].paramvv[j])and (paramestS[j].q750 >= enreg2[p].paramvv[j])) cov50_OS[j] += atest;
	}
	//cout<<("5\n");
	///////////////// factors 2
	for (int k = 0; k < 3; k++) {
		for (int j = 0; j < nparamcom; j++) {
			switch (k) {
				case 0: if ((paramest[p][j].moy >= 0.5 * enreg2[p].paramvv[j])and (paramest[p][j].moy <= 2.0 * enreg2[p].paramvv[j])) fac2_O[k][j] += atest;
					if ((paramestS[j].moy >= 0.5 * enreg2[p].paramvv[j])and (paramestS[j].moy <= 2.0 * enreg2[p].paramvv[j])) fac2_OS[k][j] += atest;
					break;
				case 1: if ((paramest[p][j].med >= 0.5 * enreg2[p].paramvv[j])and (paramest[p][j].med <= 2.0 * enreg2[p].paramvv[j])) fac2_O[k][j] += atest;
					if ((paramestS[j].med >= 0.5 * enreg2[p].paramvv[j])and (paramestS[j].med <= 2.0 * enreg2[p].paramvv[j])) fac2_OS[k][j] += atest;
					break;
				case 2: if ((paramest[p][j].mod >= 0.5 * enreg2[p].paramvv[j])and (paramest[p][j].mod <= 2.0 * enreg2[p].paramvv[j])) fac2_O[k][j] += atest;
					if ((paramestS[j].mod >= 0.5 * enreg2[p].paramvv[j])and (paramestS[j].mod <= 2.0 * enreg2[p].paramvv[j])) fac2_OS[k][j] += atest;
					break;
			}
		}
	}
	//cout<<("6\n");
	/////////////////// medianes du biais relatif
	for (int k = 0; k < 3; k++) {
		for (int j = 0; j < nparamcom; j++) {
			switch (k) {
				case 0: bmed_O[k][j][p] = (paramest[p][j].moy - enreg2[p].paramvv[j]) / enreg2[p].paramvv[j];
					bmed_OS[k][j][p] = (paramestS[j].moy - enreg2[p].paramvv[j]) / enreg2[p].paramvv[j];
					break;
				case 1: bmed_O[k][j][p] = (paramest[p][j].med - enreg2[p].paramvv[j]) / enreg2[p].paramvv[j];
					bmed_OS[k][j][p] = (paramestS[j].med - enreg2[p].paramvv[j]) / enreg2[p].paramvv[j];
					break;
				case 2: bmed_O[k][j][p] = (paramest[p][j].mod - enreg2[p].paramvv[j]) / enreg2[p].paramvv[j];
					bmed_OS[k][j][p] = (paramestS[j].mod - enreg2[p].paramvv[j]) / enreg2[p].paramvv[j];
					break;
			}
		}
		//Il restera à calculer la médiane des bmed[k][j]
	}
	//cout<<("7\n");
	///////////////// Relative Median Absolute Deviation
	long double* cc;
	cc = new long double[nsel];

	for (int j = 0; j < nparamcom; j++) {
		for (int i = 0; i < nsel; i++) cc[i] = (fabs(paretoil[i][j] - enreg2[p].paramvv[j])) / enreg2[p].paramvv[j];
		bmedr_O[j][p] = cal_medL(nsel, cc);
		//Il restera à calculer la médiane des bmedr_O[j]
		for (int i = 0; i < nsel; i++) cc[i] = (fabs(simparsel[i][j] - enreg2[p].paramvv[j])) / enreg2[p].paramvv[j];
		bmedr_OS[j][p] = cal_medL(nsel, cc);
	}
	delete []cc;
	for (int i = 0; i < nsel; i++) delete [] simparsel[i];
	delete []simparsel;
	//cout<<("8\n");
	////////////  RMAE
	for (int k = 0; k < 3; k++) {
		for (int j = 0; j < nparamcom; j++) {
			switch (k) {
				case 0: bmeda_O[k][j][p] = fabs(paramest[p][j].moy - enreg2[p].paramvv[j]) / enreg2[p].paramvv[j];
					bmeda_OS[k][j][p] = fabs(paramestS[j].moy - enreg2[p].paramvv[j]) / enreg2[p].paramvv[j];
					break;
				case 1: bmeda_O[k][j][p] = fabs(paramest[p][j].med - enreg2[p].paramvv[j]) / enreg2[p].paramvv[j];
					bmeda_OS[k][j][p] = fabs(paramestS[j].med - enreg2[p].paramvv[j]) / enreg2[p].paramvv[j];
					break;
				case 2: bmeda_O[k][j][p] = fabs(paramest[p][j].mod - enreg2[p].paramvv[j]) / enreg2[p].paramvv[j];
					bmeda_OS[k][j][p] = fabs(paramestS[j].mod - enreg2[p].paramvv[j]) / enreg2[p].paramvv[j];
					break;
			}
		}
	}
	//Il restera à calculer la médiane des bmeda_O[k][j]
	//cout<<"fin de biaisrelO\n";
}

/**
* calcule les différentes statistiques de biais, rmse... des paramètres compo
*/
void biaisrelC(int ntest, int nsel, int p) {
	long double s, d, ds;
	simparcomposel = echantillon(nsel, nparcompo, simparcompo);
	if (p > 0) delete []paramestcompoS;
	paramestcompoS = calparstat(nsel, nparcompo, simparcomposel);
	for (int j = 0; j < nparcompo; j++) {
		paramestcompomoyS[j] += paramestcompoS[j].moy;
		paramestcompomedS[j] += paramestcompoS[j].med;
	}
	//if(debuglevel==11) cout<<"biaisrelC 1\n";
	//////////////// mean relative bias
	for (int k = 0; k < 3; k++) {
		for (int j = 0; j < nparcompo; j++) {
			switch (k) {
				case 0: d = paramestcompo[p][j].moy - enreg2[p].paramvvC[j];
					ds = paramestcompoS[j].moy - enreg2[p].paramvvC[j];
					break;
				case 1: d = paramestcompo[p][j].med - enreg2[p].paramvvC[j];
					ds = paramestcompoS[j].med - enreg2[p].paramvvC[j];
					break;
				case 2: d = paramestcompo[p][j].mod - enreg2[p].paramvvC[j];
					ds = paramestcompoS[j].mod - enreg2[p].paramvvC[j];
					break;
			}
			br_C[k][j] += d / enreg2[p].paramvvC[j]; // il restera à diviser par ntest
			br_CS[k][j] += ds / enreg2[p].paramvvC[j]; // il restera à diviser par ntest
		}
	}
	//if(debuglevel==11) cout<<"biaisrelC 2\n";
	////////////  RRMISE
	for (int j = 0; j < nparcompo; j++) {
		s = 0.0;
		for (int i = 0; i < nsel; i++) {
			d = paretoilcompo[i][j] - enreg2[p].paramvvC[j];
			s += d * d;
		}
		rrmise_C[j] += s / enreg2[p].paramvvC[j] / enreg2[p].paramvvC[j] / double(nsel);// il restera à diviser par ntest
		s = 0.0;
		for (int i = 0; i < nsel; i++) {
			d = simparcomposel[i][j] - enreg2[p].paramvvC[j];
			s += d * d;
		}
		rrmise_CS[j] += s / enreg2[p].paramvvC[j] / enreg2[p].paramvvC[j] / double(nsel);// il restera à diviser par ntest
	}
	//if(debuglevel==11) cout<<"biaisrelC 3\n";
	////////////  RMAD
	for (int j = 0; j < nparcompo; j++) {
		s = 0.0;
		for (int i = 0; i < nsel; i++) s += fabs(paretoilcompo[i][j] - enreg2[p].paramvvC[j]);
		rmad_C[j] += s / fabs(enreg2[p].paramvvC[j]) / (long double)nsel;// il restera à diviser par ntest
		s = 0.0;
		for (int i = 0; i < nsel; i++) s += fabs(simparcomposel[i][j] - enreg2[p].paramvvC[j]);
		rmad_CS[j] += s / fabs(enreg2[p].paramvvC[j]) / (long double)nsel;// il restera à diviser par ntest
	}
	//if(debuglevel==11) cout<<"biaisrelC 4\n";
	////////////  RMSE
	for (int k = 0; k < 3; k++) {
		for (int j = 0; j < nparcompo; j++) {
			switch (k) {
				case 0: d = paramestcompo[p][j].moy - enreg2[p].paramvvC[j];
					ds = paramestcompoS[j].moy - enreg2[p].paramvvC[j];
					break;
				case 1: d = paramestcompo[p][j].med - enreg2[p].paramvvC[j];
					ds = paramestcompoS[j].med - enreg2[p].paramvvC[j];
					break;
				case 2: d = paramestcompo[p][j].mod - enreg2[p].paramvvC[j];
					ds = paramestcompoS[j].mod - enreg2[p].paramvvC[j];
					break;
			}
			rmse_C[k][j] += d * d / enreg2[p].paramvvC[j] / enreg2[p].paramvvC[j];
			// il restera à prendre la racine carrée et à diviser par ntest
			rmse_CS[k][j] += ds * ds / enreg2[p].paramvvC[j] / enreg2[p].paramvvC[j];
		}
	}
	//if(debuglevel==11) cout<<"biaisrelC 5\n";
	/////////////// coverages
	long double atest = 1.0 / (long double)ntest;
	for (int j = 0; j < nparcompo; j++) {
		if ((paramestcompo[p][j].q025 <= enreg2[p].paramvvC[j])and (paramestcompo[p][j].q975 >= enreg2[p].paramvvC[j])) cov95_C[j] += atest;
		if ((paramestcompo[p][j].q250 <= enreg2[p].paramvvC[j])and (paramestcompo[p][j].q750 >= enreg2[p].paramvvC[j])) cov50_C[j] += atest;
		if ((paramestcompoS[j].q025 <= enreg2[p].paramvvC[j])and (paramestcompoS[j].q975 >= enreg2[p].paramvvC[j])) cov95_CS[j] += atest;
		if ((paramestcompoS[j].q250 <= enreg2[p].paramvvC[j])and (paramestcompoS[j].q750 >= enreg2[p].paramvvC[j])) cov50_CS[j] += atest;
	}
	//if(debuglevel==11) cout<<"biaisrelC 6\n";
	///////////////// factors 2
	for (int k = 0; k < 3; k++) {
		for (int j = 0; j < nparcompo; j++) {
			switch (k) {
				case 0: if ((paramestcompo[p][j].moy >= 0.5 * enreg2[p].paramvvC[j])and (paramestcompo[p][j].moy <= 2.0 * enreg2[p].paramvvC[j])) fac2_C[k][j] += atest;
					if ((paramestcompoS[j].moy >= 0.5 * enreg2[p].paramvvC[j])and (paramestcompoS[j].moy <= 2.0 * enreg2[p].paramvvC[j])) fac2_CS[k][j] += atest;
					break;
				case 1: if ((paramestcompo[p][j].med >= 0.5 * enreg2[p].paramvvC[j])and (paramestcompo[p][j].med <= 2.0 * enreg2[p].paramvvC[j])) fac2_C[k][j] += atest;
					if ((paramestcompoS[j].med >= 0.5 * enreg2[p].paramvvC[j])and (paramestcompoS[j].med <= 2.0 * enreg2[p].paramvvC[j])) fac2_CS[k][j] += atest;
					break;
				case 2: if ((paramestcompo[p][j].mod >= 0.5 * enreg2[p].paramvvC[j])and (paramestcompo[p][j].mod <= 2.0 * enreg2[p].paramvvC[j])) fac2_C[k][j] += atest;
					if ((paramestcompoS[j].mod >= 0.5 * enreg2[p].paramvvC[j])and (paramestcompoS[j].mod <= 2.0 * enreg2[p].paramvvC[j])) fac2_CS[k][j] += atest;
					break;
			}
		}
	}
	//if(debuglevel==11) cout<<"biaisrelC 7\n";
	/////////////////// medianes du biais relatif
	for (int k = 0; k < 3; k++) {
		for (int j = 0; j < nparcompo; j++) {
			switch (k) {
				case 0: bmed_C[k][j][p] = (paramestcompo[p][j].moy - enreg2[p].paramvvC[j]) / enreg2[p].paramvvC[j];
					bmed_CS[k][j][p] = (paramestcompoS[j].moy - enreg2[p].paramvvC[j]) / enreg2[p].paramvvC[j];
					break;
				case 1: bmed_C[k][j][p] = (paramestcompo[p][j].med - enreg2[p].paramvvC[j]) / enreg2[p].paramvvC[j];
					bmed_CS[k][j][p] = (paramestcompoS[j].med - enreg2[p].paramvvC[j]) / enreg2[p].paramvvC[j];
					break;
				case 2: bmed_C[k][j][p] = (paramestcompo[p][j].mod - enreg2[p].paramvvC[j]) / enreg2[p].paramvvC[j];
					bmed_CS[k][j][p] = (paramestcompoS[j].mod - enreg2[p].paramvvC[j]) / enreg2[p].paramvvC[j];
					break;
			}
		}
		//Il restera à calculer la médiane des bmed[k][j]
	}
	//if(debuglevel==11) cout<<"biaisrelC 8\n";
	///////////////// Relative Median Absolute Deviation
	long double* cc;
	cc = new long double[nsel];

	for (int j = 0; j < nparcompo; j++) {
		for (int i = 0; i < nsel; i++) cc[i] = (fabs(paretoilcompo[i][j] - enreg2[p].paramvvC[j])) / enreg2[p].paramvvC[j];
		bmedr_C[j][p] = cal_medL(nsel, cc);
		//Il restera à calculer la médiane des bmedr_C[j]
		for (int i = 0; i < nsel; i++) cc[i] = (fabs(simparcomposel[i][j] - enreg2[p].paramvvC[j])) / enreg2[p].paramvvC[j];
		bmedr_CS[j][p] = cal_medL(nsel, cc);
	}
	delete []cc;
	for (int i = 0; i < nsel; i++) delete [] simparcomposel[i];
	delete []simparcomposel;
	//if(debuglevel==11) cout<<"biaisrelC 9\n";
	////////////  RMAE
	for (int k = 0; k < 3; k++) {
		for (int j = 0; j < nparcompo; j++) {
			switch (k) {
				case 0: bmeda_C[k][j][p] = fabs(paramestcompo[p][j].moy - enreg2[p].paramvvC[j]) / enreg2[p].paramvvC[j];
					bmeda_CS[k][j][p] = fabs(paramestcompoS[j].moy - enreg2[p].paramvvC[j]) / enreg2[p].paramvvC[j];
					break;
				case 1: bmeda_C[k][j][p] = fabs(paramestcompo[p][j].med - enreg2[p].paramvvC[j]) / enreg2[p].paramvvC[j];
					bmeda_CS[k][j][p] = fabs(paramestcompoS[j].med - enreg2[p].paramvvC[j]) / enreg2[p].paramvvC[j];
					break;
				case 2: bmeda_C[k][j][p] = fabs(paramestcompo[p][j].mod - enreg2[p].paramvvC[j]) / enreg2[p].paramvvC[j];
					bmeda_CS[k][j][p] = fabs(paramestcompoS[j].mod - enreg2[p].paramvvC[j]) / enreg2[p].paramvvC[j];
					break;
			}
		}
	}
	//Il restera à calculer la médiane des bmeda_O[k][j]
	//if(debuglevel==11) cout<<"biaisrelC 10\n";
}


/**
* calcule les différentes statistiques de biais, rmse... des paramètres scaled
*/
void biaisrelS(int ntest, int nsel, int p) {
	//cout<<"debut de biaisrelS\n";
	long double s, d, ds;
	simparscaledsel = echantillon(nsel, nparscaled, simparscaled);
	if (p > 0) delete []paramestscaledS;
	paramestscaledS = calparstat(nsel, nparscaled, simparscaledsel);
	for (int j = 0; j < nparscaled; j++) {
		paramestscaledmoyS[j] += paramestscaledS[j].moy;
		paramestscaledmedS[j] += paramestscaledS[j].med;
	}
	//////////////// mean relative bias
	for (int k = 0; k < 3; k++) {
		for (int j = 0; j < nparscaled; j++) {
			switch (k) {
				case 0: d = paramestscaled[p][j].moy - enreg2[p].paramvvS[j];
					ds = paramestscaledS[j].moy - enreg2[p].paramvvS[j];
					break;
				case 1: d = paramestscaled[p][j].med - enreg2[p].paramvvS[j];
					ds = paramestscaledS[j].med - enreg2[p].paramvvS[j];
					break;
				case 2: d = paramestscaled[p][j].mod - enreg2[p].paramvvS[j];
					ds = paramestscaledS[j].mod - enreg2[p].paramvvS[j];
					break;
			}
			br_S[k][j] += d / enreg2[p].paramvvS[j]; // il restera à diviser par ntest
			br_SS[k][j] += ds / enreg2[p].paramvvS[j]; // il restera à diviser par ntest
		}
	}

	////////////  RRMISE
	for (int j = 0; j < nparscaled; j++) {
		s = 0.0;
		for (int i = 0; i < nsel; i++) {
			d = paretoilscaled[i][j] - enreg2[p].paramvvS[j];
			s += d * d;
		}
		rrmise_S[j] += s / enreg2[p].paramvvS[j] / enreg2[p].paramvvS[j] / double(nsel);// il restera à diviser par ntest
		s = 0.0;
		for (int i = 0; i < nsel; i++) {
			d = simparscaledsel[i][j] - enreg2[p].paramvvS[j];
			s += d * d;
		}
		rrmise_SS[j] += s / enreg2[p].paramvvS[j] / enreg2[p].paramvvS[j] / double(nsel);// il restera à diviser par ntest
	}

	////////////  RMAD
	for (int j = 0; j < nparscaled; j++) {
		s = 0.0;
		for (int i = 0; i < nsel; i++) s += fabs(paretoilscaled[i][j] - enreg2[p].paramvvS[j]);
		rmad_S[j] += s / fabs(enreg2[p].paramvvS[j]) / (long double)nsel;// il restera à diviser par ntest
		s = 0.0;
		for (int i = 0; i < nsel; i++) s += fabs(simparscaledsel[i][j] - enreg2[p].paramvvS[j]);
		rmad_SS[j] += s / fabs(enreg2[p].paramvvS[j]) / (long double)nsel;// il restera à diviser par ntest
	}
	////////////  RMSE
	for (int k = 0; k < 3; k++) {
		for (int j = 0; j < nparscaled; j++) {
			switch (k) {
				case 0: d = paramestscaled[p][j].moy - enreg2[p].paramvvS[j];
					ds = paramestscaledS[j].moy - enreg2[p].paramvvS[j];
					break;
				case 1: d = paramestscaled[p][j].med - enreg2[p].paramvvS[j];
					ds = paramestscaledS[j].med - enreg2[p].paramvvS[j];
					break;
				case 2: d = paramestscaled[p][j].mod - enreg2[p].paramvvS[j];
					ds = paramestscaledS[j].mod - enreg2[p].paramvvS[j];
					break;
			}
			rmse_S[k][j] += d * d / enreg2[p].paramvvS[j] / enreg2[p].paramvvS[j];
			// il restera à prendre la racine carrée et à diviser par ntest
			rmse_SS[k][j] += ds * ds / enreg2[p].paramvvS[j] / enreg2[p].paramvvS[j];
		}
	}
	/////////////// coverages
	long double atest = 1.0 / (long double)ntest;
	for (int j = 0; j < nparscaled; j++) {
		if ((paramestscaled[p][j].q025 <= enreg2[p].paramvvS[j])and (paramestscaled[p][j].q975 >= enreg2[p].paramvvS[j])) cov95_S[j] += atest;
		if ((paramestscaled[p][j].q250 <= enreg2[p].paramvvS[j])and (paramestscaled[p][j].q750 >= enreg2[p].paramvvS[j])) cov50_S[j] += atest;
		if ((paramestscaledS[j].q025 <= enreg2[p].paramvvS[j])and (paramestscaledS[j].q975 >= enreg2[p].paramvvS[j])) cov95_SS[j] += atest;
		if ((paramestscaledS[j].q250 <= enreg2[p].paramvvS[j])and (paramestscaledS[j].q750 >= enreg2[p].paramvvS[j])) cov50_SS[j] += atest;
	}
	///////////////// factors 2
	for (int k = 0; k < 3; k++) {
		for (int j = 0; j < nparscaled; j++) {
			switch (k) {
				case 0: if ((paramestscaled[p][j].moy >= 0.5 * enreg2[p].paramvvS[j])and (paramestscaled[p][j].moy <= 2.0 * enreg2[p].paramvvS[j])) fac2_S[k][j] += atest;
					if ((paramestscaledS[j].moy >= 0.5 * enreg2[p].paramvvS[j])and (paramestscaledS[j].moy <= 2.0 * enreg2[p].paramvvS[j])) fac2_SS[k][j] += atest;
					break;
				case 1: if ((paramestscaled[p][j].med >= 0.5 * enreg2[p].paramvvS[j])and (paramestscaled[p][j].med <= 2.0 * enreg2[p].paramvvS[j])) fac2_S[k][j] += atest;
					if ((paramestscaledS[j].med >= 0.5 * enreg2[p].paramvvS[j])and (paramestscaledS[j].med <= 2.0 * enreg2[p].paramvvS[j])) fac2_SS[k][j] += atest;
					break;
				case 2: if ((paramestscaled[p][j].mod >= 0.5 * enreg2[p].paramvvS[j])and (paramestscaled[p][j].mod <= 2.0 * enreg2[p].paramvvS[j])) fac2_S[k][j] += atest;
					if ((paramestscaledS[j].mod >= 0.5 * enreg2[p].paramvvS[j])and (paramestscaledS[j].mod <= 2.0 * enreg2[p].paramvvS[j])) fac2_SS[k][j] += atest;
					break;
			}
		}
	}
	/////////////////// medianes du biais relatif
	for (int k = 0; k < 3; k++) {
		for (int j = 0; j < nparscaled; j++) {
			switch (k) {
				case 0: bmed_S[k][j][p] = (paramestscaled[p][j].moy - enreg2[p].paramvvS[j]) / enreg2[p].paramvvS[j];
					bmed_SS[k][j][p] = (paramestscaledS[j].moy - enreg2[p].paramvvS[j]) / enreg2[p].paramvvS[j];
					break;
				case 1: bmed_S[k][j][p] = (paramestscaled[p][j].med - enreg2[p].paramvvS[j]) / enreg2[p].paramvvS[j];
					bmed_SS[k][j][p] = (paramestscaledS[j].med - enreg2[p].paramvvS[j]) / enreg2[p].paramvvS[j];
					break;
				case 2: bmed_S[k][j][p] = (paramestscaled[p][j].mod - enreg2[p].paramvvS[j]) / enreg2[p].paramvvS[j];
					bmed_SS[k][j][p] = (paramestscaledS[j].mod - enreg2[p].paramvvS[j]) / enreg2[p].paramvvS[j];
					break;
			}
		}
		//Il restera à calculer la médiane des bmed[k][j]
	}
	///////////////// Relative Median Absolute Deviation
	long double* cc;
	cc = new long double[nsel];

	for (int j = 0; j < nparscaled; j++) {
		for (int i = 0; i < nsel; i++) cc[i] = (fabs(paretoilscaled[i][j] - enreg2[p].paramvvS[j])) / enreg2[p].paramvvS[j];
		bmedr_S[j][p] = cal_medL(nsel, cc);
		//Il restera à calculer la médiane des bmedr_S[j]
		for (int i = 0; i < nsel; i++) cc[i] = (fabs(simparscaledsel[i][j] - enreg2[p].paramvvS[j])) / enreg2[p].paramvvS[j];
		bmedr_SS[j][p] = cal_medL(nsel, cc);
	}
	delete []cc;
	for (int i = 0; i < nsel; i++) delete []simparscaledsel[i];
	delete simparscaledsel;
	////////////  RMAE
	for (int k = 0; k < 3; k++) {
		for (int j = 0; j < nparscaled; j++) {
			switch (k) {
				case 0: bmeda_S[k][j][p] = fabs(paramestscaled[p][j].moy - enreg2[p].paramvvS[j]) / enreg2[p].paramvvS[j];
					bmeda_SS[k][j][p] = fabs(paramestscaledS[j].moy - enreg2[p].paramvvS[j]) / enreg2[p].paramvvS[j];
					break;
				case 1: bmeda_S[k][j][p] = fabs(paramestscaled[p][j].med - enreg2[p].paramvvS[j]) / enreg2[p].paramvvS[j];
					bmeda_SS[k][j][p] = fabs(paramestscaledS[j].med - enreg2[p].paramvvS[j]) / enreg2[p].paramvvS[j];
					break;
				case 2: bmeda_S[k][j][p] = fabs(paramestscaled[p][j].mod - enreg2[p].paramvvS[j]) / enreg2[p].paramvvS[j];;
					bmeda_SS[k][j][p] = fabs(paramestscaledS[j].mod - enreg2[p].paramvvS[j]) / enreg2[p].paramvvS[j];
					break;
			}
		}
	}
	//Il restera à calculer la médiane des bmeda_O[k][j]
	//cout<<"fin de biaisrelS\n";
}

void finbiaisrelO(int ntest) {
	for (int j = 0; j < nparamcom; j++) {
		paramestmoyS[j] /= (long double)ntest;
		paramestmedS[j] /= (long double)ntest;
	}
	for (int k = 0; k < 3; k++) {
		for (int j = 0; j < nparamcom; j++) {
			br_O[k][j] /= (long double)ntest;
			br_OS[k][j] /= (long double)ntest;
		}
	}
	for (int j = 0; j < nparamcom; j++) {
		rrmise_O[j] = sqrt(rrmise_O[j] / (long double)ntest);
		rrmise_OS[j] = sqrt(rrmise_OS[j] / (long double)ntest);
	}
	for (int j = 0; j < nparamcom; j++) {
		rmad_O[j] /= (long double)ntest;
		rmad_OS[j] /= (long double)ntest;
	}
	for (int k = 0; k < 3; k++) {
		for (int j = 0; j < nparamcom; j++) {
			rmse_O[k][j] = sqrt(rmse_O[k][j] / (long double)ntest);
			rmse_OS[k][j] = sqrt(rmse_OS[k][j] / (long double)ntest);
		}
	}
	for (int k = 0; k < 3; k++) {
		for (int j = 0; j < nparamcom; j++) {
			rmb_O[k][j] = cal_medL(ntest, bmed_O[k][j]);
			rmb_OS[k][j] = cal_medL(ntest, bmed_OS[k][j]);
		}
	}
	for (int j = 0; j < nparamcom; j++) {
		rmedad_O[j] = cal_medL(ntest, bmedr_O[j]);
		rmedad_OS[j] = cal_medL(ntest, bmedr_OS[j]);
	}
	for (int k = 0; k < 3; k++) {
		for (int j = 0; j < nparamcom; j++) {
			rmae_O[k][j] = cal_medL(ntest, bmeda_O[k][j]);
			rmae_OS[k][j] = cal_medL(ntest, bmeda_OS[k][j]);
		}
	}
}

void finbiaisrelC(int ntest) {
	if (debuglevel == 11) cerr << "debut de finbiaisrelC\n";
	for (int j = 0; j < nparamcom; j++) {
		paramestcompomoyS[j] /= (long double)ntest;
		paramestcompomedS[j] /= (long double)ntest;
	}
	if (debuglevel == 11) cerr << "finbiaisrelC 1\n";
	for (int k = 0; k < 3; k++) {
		for (int j = 0; j < nparcompo; j++) {
			br_C[k][j] /= (long double)ntest;
			br_CS[k][j] /= (long double)ntest;
		}
	}
	if (debuglevel == 11) cerr << "finbiaisrelC 2\n";
	for (int j = 0; j < nparcompo; j++) {
		rrmise_C[j] = sqrt(rrmise_C[j] / (long double)ntest);
		rrmise_CS[j] = sqrt(rrmise_CS[j] / (long double)ntest);
	}
	if (debuglevel == 11) cerr << "finbiaisrelC 3\n";
	for (int j = 0; j < nparcompo; j++) {
		rmad_C[j] /= (long double)ntest;
		rmad_CS[j] /= (long double)ntest;
	}
	if (debuglevel == 11) cerr << "finbiaisrelC 4\n";
	for (int k = 0; k < 3; k++) {
		for (int j = 0; j < nparcompo; j++) {
			rmse_C[k][j] = sqrt(rmse_C[k][j] / (long double)ntest);
			rmse_CS[k][j] = sqrt(rmse_CS[k][j] / (long double)ntest);
		}
	}
	if (debuglevel == 11) cerr << "finbiaisrelC 5\n";
	for (int k = 0; k < 3; k++) {
		for (int j = 0; j < nparcompo; j++) {
			rmb_C[k][j] = cal_medL(ntest, bmed_C[k][j]);
			rmb_CS[k][j] = cal_medL(ntest, bmed_CS[k][j]);
		}
	}
	if (debuglevel == 11) cerr << "finbiaisrelC 6\n";
	for (int j = 0; j < nparcompo; j++) {
		rmedad_C[j] = cal_medL(ntest, bmedr_C[j]);
		rmedad_CS[j] = cal_medL(ntest, bmedr_CS[j]);
	}
	if (debuglevel == 11) cerr << "finbiaisrelC 7\n";
	for (int k = 0; k < 3; k++) {
		for (int j = 0; j < nparcompo; j++) {
			rmae_C[k][j] = cal_medL(ntest, bmeda_C[k][j]);
			rmae_CS[k][j] = cal_medL(ntest, bmeda_CS[k][j]);
		}
	}
	if (debuglevel == 11) cerr << "fin de finbiaisrelC\n";
}

void finbiaisrelS(int ntest) {
	for (int j = 0; j < nparamcom; j++) {
		paramestscaledmoyS[j] /= (long double)ntest;
		paramestscaledmedS[j] /= (long double)ntest;
	}
	for (int k = 0; k < 3; k++) {
		for (int j = 0; j < nparscaled; j++) {
			br_S[k][j] /= (long double)ntest;
			br_SS[k][j] /= (long double)ntest;
		}
	}
	for (int j = 0; j < nparscaled; j++) {
		rrmise_S[j] = sqrt(rrmise_S[j] / (long double)ntest);
		rrmise_SS[j] = sqrt(rrmise_SS[j] / (long double)ntest);
	}
	for (int j = 0; j < nparscaled; j++) {
		rmad_S[j] /= (long double)ntest;
		rmad_SS[j] /= (long double)ntest;
	}
	for (int k = 0; k < 3; k++) {
		for (int j = 0; j < nparscaled; j++) {
			rmse_S[k][j] = sqrt(rmse_S[k][j] / (long double)ntest);
			rmse_SS[k][j] = sqrt(rmse_SS[k][j] / (long double)ntest);
		}
	}
	for (int k = 0; k < 3; k++) {
		for (int j = 0; j < nparscaled; j++) {
			rmb_S[k][j] = cal_medL(ntest, bmed_S[k][j]);
			rmb_SS[k][j] = cal_medL(ntest, bmed_SS[k][j]);
		}
	}
	for (int j = 0; j < nparscaled; j++) {
		rmedad_S[j] = cal_medL(ntest, bmedr_S[j]);
		rmedad_SS[j] = cal_medL(ntest, bmedr_SS[j]);
	}
	for (int k = 0; k < 3; k++) {
		for (int j = 0; j < nparscaled; j++) {
			rmae_S[k][j] = cal_medL(ntest, bmeda_S[k][j]);
			rmae_SS[k][j] = cal_medL(ntest, bmeda_SS[k][j]);
		}
	}
}

string fmLD(long double x, int largeur, int precision) {
	string s;
	size_t pos;
	int lon;
	if (x < 1E80) {
		s = LongDoubleToString(x);
		if (s.find(".") == string::npos) s = s + ".0";
		s = s + "000000000000000";
		//cout<<"x="<<x<<"   s="<<s;
		pos = s.find(".");
		s = s.substr(0, pos + precision + 1);
		//cout<<"   s'="<<s;
		lon = s.length();
		while (lon < largeur) {
			s = " " + s;
			lon++;
		}
		//cout<<"   s final="<<s<<"\n";
	}
	else {
		s = "+inf";
		lon = s.length();
		while (lon < largeur) {
			s = " " + s;
			lon++;
		}
		valinfinie = true;
		cout << "valinfinie=" << valinfinie << "\n";
	}
	return s;
}

string fmLDP(long double x, int largeur, int precision) {
	string s;
	size_t pos;
	int lon;
	if (x < 1E80) {
		s = LongDoubleToString(x);
		if (s.find(".") == string::npos) s = s + ".0";
		s = s + "000000000000000";
		//cout<<"x="<<x<<"   s="<<s;
		pos = s.find(".");
		s = "(" + s.substr(0, pos + precision + 1) + ")";
		//cout<<"   s'="<<s;
		lon = s.length();
		while (lon < largeur) {
			s = " " + s;
			lon++;
		}
		//cout<<"   s final="<<s<<"\n";
	}
	else {
		s = "+inf";
		lon = s.length();
		while (lon < largeur) {
			s = " " + s;
			lon++;
		}
		valinfinie = true;
		cout << "valinfinie=" << valinfinie << "\n";
	}
	return s;
}


/**
 *   Mise en forme des résultats pour les paramètres originaux
 */
void ecriresO(int ntest, int nsel) {
	time_t rawtime;
	struct tm* timeinfo;
	time(&rawtime);
	timeinfo = localtime(&rawtime);
	string nomfiresult;
	long double *x, *mo;
	x = new long double[ntest];
	mo = new long double[4];
	nomfiresult = path + ident + "_bias_original.txt";
	cout << "Les résultats sont dans " << nomfiresult << "\n";
	ofstream f1;
	f1.open(nomfiresult.c_str());
	f1 << "DIYABC :                 Bias and Mean Square Error Analysis                         " << asctime(timeinfo) << "\n";
	f1 << "Data file       : " << header.datafilename << "\n";
	f1 << "Reference table : " << rt.filename << "\n";
	switch (numtransf) {
		case 1: f1 << "No transformation of parameters\n";
			break;
		case 2: f1 << "Transformation LOG of parameters\n";
			break;
		case 3: f1 << "Transformation LOGIT of parameters\n";
			break;
		case 4: f1 << "Transformation LOG(TG) of parameters\n";
			break;
	}
	f1 << "Chosen scenario : " << rt.scenteste << "\n";
	f1 << "Number of simulated data sets : " << rt.nreclus << "\n";
	f1 << "Number of selected data sets  : " << nsel << "\n";
	f1 << "Results based on " << ntest << " test data sets\n\n";
	if (valinfinie) f1 << "Some statisics are noted '+inf' due to zero values of drawn parameters. Check prior minimum values.\n";
	if (valinfinie) cout << "VALEUR INFINIE\n";
	else cout << "PAS DE VALEUR INFINIE\n";
	f1 << "\n";
	f1 << "                                                                        Averages\n";
	f1 << "Parameter                True values               Means                 Medians                 Modes\n";
	f1 << setiosflags(ios::scientific) << setiosflags(ios::right);
	for (int j = 0; j < nparamcom; j++) {
		//cout<<nomparam[j]<<"\n";
		f1 << nomparamO[j];
		for (int i = 0; i < 24 - (int)nomparamO[j].length(); i++) f1 << " ";
		for (int i = 0; i < ntest; i++) x[i] = enreg2[i].paramvv[j];
		mo[0] = cal_moyL(ntest, x);
		for (int i = 0; i < ntest; i++) x[i] = paramest[i][j].moy;
		mo[1] = cal_moyL(ntest, x);
		for (int i = 0; i < ntest; i++) x[i] = paramest[i][j].med;
		mo[2] = cal_moyL(ntest, x);
		for (int i = 0; i < ntest; i++) x[i] = paramest[i][j].mod;
		mo[3] = cal_moyL(ntest, x);
		f1 << setw(10) << setprecision(3) << mo[0];
		f1 << setw(24) << setprecision(3) << mo[1];
		f1 << setw(24) << setprecision(3) << mo[2];
		f1 << setw(22) << setprecision(3) << mo[3] << "\n";
		for (int i = 0; i < 35; i++) f1 << " ";
		f1 << setw(14) << setprecision(3) << "(" << paramestmoyS[j] << ")";
		f1 << setw(14) << setprecision(3) << "(" << paramestmedS[j] << ")\n";
		//f1<<setw(9)<<setprecision(3)<<"("<<paramestS[j].mod<<")\n";
	}
	f1 << "\n                                           Mean Relative Bias\n";
	f1 << "Parameter                     Means                    Medians                    Modes\n";
	//f1<<setiosflags(ios::fixed);
	for (int j = 0; j < nparamcom; j++) {
		//cout<<nomparam[j]<<"\n";
		f1 << nomparamO[j];
		for (int i = 0; i < 24 - (int)nomparamO[j].length(); i++) f1 << " ";
		f1 << fmLD(br_O[0][j], 7, 3) << " " << fmLDP(br_OS[0][j], 8, 3);
		f1 << fmLD(br_O[1][j], 18, 4) << " " << fmLDP(br_OS[1][j], 8, 3);;
		f1 << fmLD(br_O[2][j], 20, 4) << "\n";
		/*for (int i=0;i<34;i++) f1<<" ";
		f1<<fmLDP(br_OS[0][j],8,4);
		f1<<fmLDP(br_OS[1][j],15,4);
		f1<<fmLDP(br_OS[2][j],15,4)<<"\n";*/
	}
	f1 << "\n                            RRMISE              RMeanAD            Square root of mean square error/true value\n";
	f1 << "Parameter                                                           Mean                Median            Mode\n";
	for (int j = 0; j < nparamcom; j++) {
		//cout<<nomparam[j]<<"\n";
		f1 << nomparamO[j];
		for (int i = 0; i < 20 - (int)nomparamO[j].length(); i++) f1 << " ";
		f1 << fmLD(rrmise_O[j], 10, 3) << " " << fmLDP(rrmise_OS[j], 8, 3);
		f1 << fmLD(rmad_O[j], 11, 3) << " " << fmLDP(rmad_OS[j], 8, 3);
		f1 << fmLD(rmse_O[0][j], 11, 3) << " " << fmLDP(rmse_OS[0][j], 8, 3);
		f1 << fmLD(rmse_O[1][j], 11, 3) << " " << fmLDP(rmse_OS[1][j], 8, 3);
		f1 << fmLD(rmse_O[2][j], 11, 3) << "\n";
		/*for (int i=0;i<25;i++) f1<<" ";
		f1<<fmLDP(rrmise_OS[j],10,3);
		f1<<fmLDP(rmad_OS[j],18,3);
		f1<<fmLDP(rmse_OS[0][j],18,3);
		f1<<fmLDP(rmse_OS[1][j],19,3);
		f1<<fmLDP(rmse_OS[2][j],17,3)<<"\n";*/
	}
	f1 << "\n                                                                   Factor 2          Factor 2           Factor 2\n";
	f1 << "Parameter               50% Coverage        95% Coverage           (Mean)            (Median)            (Mode)  \n";
	for (int j = 0; j < nparamcom; j++) {
		//cout<<nomparam[j]<<"\n";
		f1 << nomparamO[j];
		for (int i = 0; i < 19 - (int)nomparamO[j].length(); i++) f1 << " ";
		f1 << fmLD(cov50_O[j], 10, 3) << " " << fmLDP(cov50_OS[j], 7, 3);
		f1 << fmLD(cov95_O[j], 12, 3) << " " << fmLDP(cov95_OS[j], 7, 3);
		f1 << fmLD(fac2_O[0][j], 12, 3) << " " << fmLDP(fac2_OS[0][j], 7, 3);
		f1 << fmLD(fac2_O[1][j], 11, 3) << " " << fmLDP(fac2_OS[1][j], 7, 3);
		f1 << fmLD(fac2_O[2][j], 15, 3) << "\n";
		/*for (int i=0;i<25;i++) f1<<" ";
		f1<<fmLDP(cov50_OS[j],10,3);
		f1<<fmLDP(cov95_OS[j],18,3);
		f1<<fmLDP(fac2_OS[0][j],18,3);
		f1<<fmLDP(fac2_OS[1][j],19,3);
		f1<<fmLDP(fac2_OS[2][j],17,3)<<"\n";*/
	}
	f1 << "\n                                          Median Relative Bias\n";
	f1 << "Parameter                 Means                    Medians                  Modes\n";
	for (int j = 0; j < nparamcom; j++) {
		//cout<<nomparam[j]<<"\n";
		f1 << nomparamO[j];
		for (int i = 0; i < 19 - (int)nomparamO[j].length(); i++) f1 << " ";
		f1 << fmLD(rmb_O[0][j], 8, 3) << " " << fmLDP(rmb_OS[0][j], 8, 3);
		f1 << fmLD(rmb_O[1][j], 18, 3) << " " << fmLDP(rmb_OS[1][j], 8, 3);
		f1 << fmLD(rmb_O[2][j], 18, 3) << "\n";
		/*for (int i=0;i<44;i++) f1<<" ";
			f1<<fmLDP(rmb_OS[0][j],8,3);
			f1<<fmLDP(rmb_OS[1][j],15,3);
			f1<<fmLDP(rmb_OS[2][j],15,3)<<"\n";*/
	}
	f1 << "\n                             RMedAD                  Median of the absolute error/true value\n";
	f1 << "Parameter                                          Means                    Medians                  Modes\n";
	for (int j = 0; j < nparamcom; j++) {
		//cout<<nomparam[j]<<"\n";
		f1 << nomparamO[j];
		for (int i = 0; i < 20 - (int)nomparamO[j].length(); i++)f1 << " ";
		f1 << fmLD(rmedad_O[j], 10, 3) << " " << fmLDP(rmedad_OS[j], 8, 3);
		f1 << fmLD(rmae_O[0][j], 15, 3) << " " << fmLDP(rmae_OS[0][j], 8, 3);
		f1 << fmLD(rmae_O[1][j], 15, 3) << " " << fmLDP(rmae_OS[1][j], 8, 3);
		f1 << fmLD(rmae_O[2][j], 19, 3) << "\n";
		/*for (int i=0;i<25;i++) f1<<" ";
		f1<<fmLDP(rmedad_OS[j],11,3);
		f1<<fmLDP(rmae_OS[0][j],16,3);
		f1<<fmLDP(rmae_OS[1][j],16,3);
		f1<<fmLDP(rmae_OS[2][j],16,3)<<"\n";*/
	}
	if (valinfinie) f1 << "\nN.B. Some statisics are noted '+inf' due to zero values of drawn parameters. Check prior minimum values.\n";
	if (valinfinie) cout << "VALEUR INFINIE\n";
	else cout << "PAS DE VALEUR INFINIE\n";
	f1 << "\n";
	valinfinie = false;
	f1.close();
}

/**
 *   Mise en forme des résultats pour les paramètres composites
 */
void ecriresC(int ntest, int nsel) {
	time_t rawtime;
	struct tm* timeinfo;
	time(&rawtime);
	timeinfo = localtime(&rawtime);
	string nomfiresult;
	long double *x, *mo;
	x = new long double[ntest];
	mo = new long double[4];
	nomfiresult = path + ident + "_bias_composite.txt";
	//strcpy(nomfiresult,path);
	//strcat(nomfiresult,ident);
	//strcat(nomfiresult,"_bias.txt");
	cout << "Les résultats sont dans " << nomfiresult << "\n";
	ofstream f1;
	f1.open(nomfiresult.c_str());
	f1 << "DIYABC :                 Bias and Mean Square Error Analysis                         " << asctime(timeinfo) << "\n";
	f1 << "Data file       : " << header.datafilename << "\n";
	f1 << "Reference table : " << rt.filename << "\n";
	switch (numtransf) {
		case 1: f1 << "No transformation of parameters\n";
			break;
		case 2: f1 << "Transformation LOG of parameters\n";
			break;
		case 3: f1 << "Transformation LOGIT of parameters\n";
			break;
		case 4: f1 << "Transformation LOG(TG) of parameters\n";
			break;
	}
	f1 << "Chosen scenario : " << rt.scenteste << "\n";
	f1 << "Number of simulated data sets : " << rt.nreclus << "\n";
	f1 << "Number of selected data sets  : " << nsel << "\n";
	f1 << "Results based on " << ntest << " test data sets\n";
	if (valinfinie) f1 << "Some statisics are noted '+inf' due to zero values of drawn parameters. Check prior minimum values.\n";
	f1 << "\n";
	f1 << "                                                                        Averages\n";
	f1 << "Parameter                True values               Means                 Medians                 Modes\n";
	f1 << setiosflags(ios::scientific) << setiosflags(ios::right);
	for (int j = 0; j < nparcompo; j++) {
		//cout<<nomparam[j]<<"\n";
		f1 << nomparamC[j];
		for (int i = 0; i < 24 - (int)nomparamC[j].length(); i++) f1 << " ";
		for (int i = 0; i < ntest; i++) x[i] = enreg2[i].paramvvC[j];
		mo[0] = cal_moyL(ntest, x);
		for (int i = 0; i < ntest; i++) x[i] = paramestcompo[i][j].moy;
		mo[1] = cal_moyL(ntest, x);
		for (int i = 0; i < ntest; i++) x[i] = paramestcompo[i][j].med;
		mo[2] = cal_moyL(ntest, x);
		for (int i = 0; i < ntest; i++) x[i] = paramestcompo[i][j].mod;
		mo[3] = cal_moyL(ntest, x);
		f1 << setw(10) << setprecision(3) << mo[0];
		f1 << setw(24) << setprecision(3) << mo[1];
		f1 << setw(24) << setprecision(3) << mo[2];
		f1 << setw(22) << setprecision(3) << mo[3] << "\n";
		for (int i = 0; i < 35; i++) f1 << " ";
		f1 << setw(14) << setprecision(3) << "(" << paramestcompoS[j].moy << ")";
		f1 << setw(14) << setprecision(3) << "(" << paramestcompoS[j].med << ")\n";
	}
	f1 << "\n                                           Mean Relative Bias\n";
	f1 << "Parameter                     Means                    Medians                    Modes\n";
	for (int j = 0; j < nparcompo; j++) {
		//cout<<nomparam[j]<<"\n";
		f1 << nomparamC[j];
		for (int i = 0; i < 24 - (int)nomparamC[j].length(); i++) f1 << " ";
		f1 << fmLD(br_C[0][j], 7, 3) << " " << fmLDP(br_CS[0][j], 8, 3);
		f1 << fmLD(br_C[1][j], 18, 4) << " " << fmLDP(br_CS[1][j], 8, 3);
		f1 << fmLD(br_C[2][j], 20, 4) << "\n";
		//for (int i=0;i<34;i++) f1<<" ";
		//f1<<fmLDP(br_CS[0][j],8,4);
		//f1<<fmLDP(br_CS[1][j],15,4);
		//f1<<fmLDP(br_CS[2][j],15,4)<<"\n";
	}
	f1 << "\n                            RRMISE              RMeanAD            Square root of mean square error/true value\n";
	f1 << "Parameter                                                           Mean                Median            Mode\n";
	for (int j = 0; j < nparcompo; j++) {
		//cout<<nomparam[j]<<"\n";
		f1 << nomparamC[j];
		for (int i = 0; i < 20 - (int)nomparamC[j].length(); i++) f1 << " ";
		f1 << fmLD(rrmise_C[j], 10, 3) << " " << fmLDP(rrmise_CS[j], 8, 3);
		f1 << fmLD(rmad_C[j], 11, 3) << " " << fmLDP(rmad_CS[j], 8, 3);
		f1 << fmLD(rmse_C[0][j], 11, 3) << " " << fmLDP(rmse_CS[0][j], 8, 3);
		f1 << fmLD(rmse_C[1][j], 11, 3) << " " << fmLDP(rmse_CS[1][j], 8, 3);
		f1 << fmLD(rmse_C[2][j], 11, 3) << "\n";
		/*for (int i=0;i<25;i++) f1<<" ";
		f1<<fmLDP(rrmise_CS[j],10,3);
		f1<<fmLDP(rmad_CS[j],18,3);
		f1<<fmLDP(rmse_CS[0][j],18,3);
		f1<<fmLDP(rmse_CS[1][j],19,3);
		f1<<fmLDP(rmse_CS[2][j],17,3)<<"\n";*/
	}
	f1 << "\n                                                                   Factor 2          Factor 2           Factor 2\n";
	f1 << "Parameter               50% Coverage        95% Coverage           (Mean)            (Median)            (Mode)  \n";
	for (int j = 0; j < nparcompo; j++) {
		//cout<<nomparam[j]<<"\n";
		f1 << nomparamC[j];
		for (int i = 0; i < 19 - (int)nomparamC[j].length(); i++) f1 << " ";
		f1 << fmLD(cov50_C[j], 10, 3) << " " << fmLDP(cov50_CS[j], 7, 3);
		f1 << fmLD(cov95_C[j], 12, 3) << " " << fmLDP(cov95_CS[j], 7, 3);
		f1 << fmLD(fac2_C[0][j], 12, 3) << " " << fmLDP(fac2_CS[0][j], 7, 3);
		f1 << fmLD(fac2_C[1][j], 11, 3) << " " << fmLDP(fac2_CS[1][j], 7, 3);
		f1 << fmLD(fac2_C[2][j], 15, 3) << "\n";
		/*for (int i=0;i<25;i++) f1<<" ";
		f1<<fmLDP(cov50_CS[j],10,3);
		f1<<fmLDP(cov95_CS[j],18,3);
		f1<<fmLDP(fac2_CS[0][j],18,3);33
		f1<<fmLDP(fac2_CS[1][j],19,3);
		f1<<fmLDP(fac2_CS[2][j],17,3)<<"\n";*/
	}
	f1 << "\n                                          Median Relative Bias\n";
	f1 << "Parameter                 Means                    Medians                  Modes\n";
	for (int j = 0; j < nparcompo; j++) {
		//cout<<nomparam[j]<<"\n";
		f1 << nomparamC[j];
		for (int i = 0; i < 19 - (int)nomparamC[j].length(); i++) f1 << " ";
		f1 << fmLD(rmb_C[0][j], 8, 3) << " " << fmLDP(rmb_CS[0][j], 8, 3);
		f1 << fmLD(rmb_C[1][j], 18, 3) << " " << fmLDP(rmb_CS[1][j], 8, 3);
		f1 << fmLD(rmb_C[2][j], 18, 3) << "\n";
		/*for (int i=0;i<44;i++) f1<<" ";
			f1<<fmLDP(rmb_CS[0][j],8,3);
			f1<<fmLDP(rmb_CS[1][j],15,3);
			f1<<fmLDP(rmb_CS[2][j],15,3)<<"\n";*/
	}
	f1 << "\n                             RMedAD                  Median of the absolute error/true value\n";
	f1 << "Parameter                                          Means                    Medians                  Modes\n";
	for (int j = 0; j < nparcompo; j++) {
		//cout<<nomparam[j]<<"\n";
		f1 << nomparamC[j];
		for (int i = 0; i < 20 - (int)nomparamC[j].length(); i++) f1 << " ";
		f1 << fmLD(rmedad_C[j], 10, 3) << " " << fmLDP(rmedad_CS[j], 8, 3);
		f1 << fmLD(rmae_C[0][j], 15, 3) << " " << fmLDP(rmae_CS[0][j], 8, 3);
		f1 << fmLD(rmae_C[1][j], 15, 3) << " " << fmLDP(rmae_CS[1][j], 8, 3);
		f1 << fmLD(rmae_C[2][j], 19, 3) << "\n";
		/*for         f1<<"Results based on "<<ntest<<" test data sets\n";
	if (valinfinie) f1<<"Some statisics are noted '+inf' due to zero values of drawn parameters. Check prior minimum values.\n";
	f1<<"\n";
(int i=0;i<25;i++) f1<<" ";
		f1<<fmLDP(rmedad_CS[j],11,3);
		f1<<fmLDP(rmae_CS[0][j],16,3);
		f1<<fmLDP(rmae_CS[1][j],16,3);
		f1<<fmLDP(rmae_CS[2][j],16,3)<<"\n";*/
	}
	if (valinfinie) f1 << "\nN.B. Some statisics are noted '+inf' due to zero values of drawn parameters. Check prior minimum values.\n";
	if (valinfinie) cout << "VALEUR INFINIE\n";
	else cout << "PAS DE VALEUR INFINIE\n";
	f1 << "\n";
	valinfinie = false;
	f1.close();

}

/**
 *   Mise en forme des résultats pour les paramètres scaled
 */
void ecriresS(int ntest, int nsel) {
	if (debuglevel == 11) cout << "debut de ecriresS\n";
	time_t rawtime;
	struct tm* timeinfo;
	time(&rawtime);
	timeinfo = localtime(&rawtime);
	string nomfiresult;
	long double *x, *mo;
	x = new long double[ntest];
	mo = new long double[4];
	nomfiresult = path + ident + "_bias_scaled.txt";
	//strcpy(nomfiresult,path);
	//strcat(nomfiresult,ident);
	//strcat(nomfiresult,"_bias.txt");
	cout << "Les résultats sont dans " << nomfiresult << "\n";
	ofstream f1;
	f1.open(nomfiresult.c_str());
	f1 << "DIYABC :                 Bias and Mean Square Error Analysis                         " << asctime(timeinfo) << "\n";
	f1 << "Data file       : " << header.datafilename << "\n";
	f1 << "Reference table : " << rt.filename << "\n";
	switch (numtransf) {
		case 1: f1 << "No transformation of parameters\n";
			break;
		case 2: f1 << "Transformation LOG of parameters\n";
			break;
		case 3: f1 << "Transformation LOGIT of parameters\n";
			break;
		case 4: f1 << "Transformation LOG(TG) of parameters\n";
			break;
	}
	f1 << "Chosen scenario : " << rt.scenteste << "\n";
	f1 << "Number of simulated data sets : " << rt.nreclus << "\n";
	f1 << "Number of selected data sets  : " << nsel << "\n";
	f1 << "Results based on " << ntest << " test data sets\n";
	if (valinfinie) f1 << "Some statisics are noted '+inf' due to zero values of drawn parameters. Check prior minimum values.\n";
	f1 << "\n";
	f1 << "                                                                        Averages\n";
	f1 << "Parameter                True values               Means                 Medians                 Modes\n";
	f1 << setiosflags(ios::scientific) << setiosflags(ios::right);
	for (int j = 0; j < nparscaled; j++) {
		//cout<<nomparam[j]<<"\n";
		f1 << nomparamS[j];
		for (int i = 0; i < 24 - (int)nomparamS[j].length(); i++) f1 << " ";
		for (int i = 0; i < ntest; i++) x[i] = enreg2[i].paramvvS[j];
		mo[0] = cal_moyL(ntest, x);
		for (int i = 0; i < ntest; i++) x[i] = paramestscaled[i][j].moy;
		mo[1] = cal_moyL(ntest, x);
		for (int i = 0; i < ntest; i++) x[i] = paramestscaled[i][j].med;
		mo[2] = cal_moyL(ntest, x);
		for (int i = 0; i < ntest; i++) x[i] = paramestscaled[i][j].mod;
		mo[3] = cal_moyL(ntest, x);
		f1 << setw(10) << setprecision(3) << mo[0];
		f1 << setw(24) << setprecision(3) << mo[1];
		f1 << setw(24) << setprecision(3) << mo[2];
		f1 << setw(22) << setprecision(3) << mo[3] << "\n";
		for (int i = 0; i < 35; i++) f1 << " ";
		f1 << setw(14) << setprecision(3) << "(" << paramestscaledS[j].moy << ")";
		f1 << setw(14) << setprecision(3) << "(" << paramestscaledS[j].med << ")\n";
	}
	f1 << "\n                                           Mean Relative Bias\n";
	f1 << "Parameter                     Means                    Medians                    Modes\n";
	for (int j = 0; j < nparscaled; j++) {
		//cout<<nomparam[j]<<"\n";
		f1 << nomparamS[j];
		for (int i = 0; i < 24 - (int)nomparamS[j].length(); i++) f1 << " ";
		f1 << fmLD(br_S[0][j], 7, 3) << " " << fmLDP(br_SS[0][j], 8, 3);
		f1 << fmLD(br_S[1][j], 18, 4) << " " << fmLDP(br_SS[1][j], 8, 3);
		f1 << fmLD(br_S[2][j], 18, 4) << "\n";
		/*for (int i=0;i<34;i++) f1<<" ";
		f1<<fmLDP(br_SS[0][j],8,4);
		f1<<fmLDP(br_SS[1][j],15,4);
		f1<<fmLDP(br_SS[2][j],15,4)<<"\n";*/
	}
	f1 << "\n                            RRMISE              RMeanAD            Square root of mean square error/true value\n";
	f1 << "Parameter                                                           Mean                Median            Mode\n";
	for (int j = 0; j < nparscaled; j++) {
		//cout<<nomparam[j]<<"\n";
		f1 << nomparamS[j];
		for (int i = 0; i < 20 - (int)nomparamS[j].length(); i++) f1 << " ";
		f1 << fmLD(rrmise_S[j], 10, 3) << " " << fmLDP(rrmise_SS[j], 8, 3);
		f1 << fmLD(rmad_S[j], 11, 3) << " " << fmLDP(rmad_SS[j], 8, 3);
		f1 << fmLD(rmse_S[0][j], 11, 3) << " " << fmLDP(rmse_SS[0][j], 8, 3);
		f1 << fmLD(rmse_S[1][j], 11, 3) << " " << fmLDP(rmse_SS[1][j], 8, 3);
		f1 << fmLD(rmse_S[2][j], 11, 3) << "\n";
		/*for (int i=0;i<25;i++) f1<<" ";
		f1<<fmLDP(rrmise_SS[j],10,3);
		f1<<fmLDP(rmad_SS[j],18,3);
		f1<<fmLDP(rmse_SS[0][j],18,3);
		f1<<fmLDP(rmse_SS[1][j],19,3);
		f1<<fmLDP(rmse_SS[2][j],17,3)<<"\n";*/
	}
	f1 << "\n                                                                   Factor 2          Factor 2           Factor 2\n";
	f1 << "Parameter               50% Coverage        95% Coverage          (Mean)            (Median)            (Mode)  \n";
	for (int j = 0; j < nparscaled; j++) {
		//cout<<nomparam[j]<<"\n";
		f1 << nomparamS[j];
		for (int i = 0; i < 19 - (int)nomparamS[j].length(); i++) f1 << " ";
		f1 << fmLD(cov50_S[j], 10, 3) << " " << fmLDP(cov50_SS[j], 7, 3);
		f1 << fmLD(cov95_S[j], 12, 3) << " " << fmLDP(cov95_SS[j], 7, 3);
		f1 << fmLD(fac2_S[0][j], 12, 3) << " " << fmLDP(fac2_SS[0][j], 7, 3);
		f1 << fmLD(fac2_S[1][j], 11, 3) << " " << fmLDP(fac2_SS[1][j], 7, 3);
		f1 << fmLD(fac2_S[2][j], 15, 3) << "\n";
		/*for (int i=0;i<25;i++) f1<<" ";
		f1<<fmLDP(cov50_SS[j],10,3);
		f1<<fmLDP(cov95_SS[j],18,3);
		f1<<fmLDP(fac2_SS[0][j],18,3);
		f1<<fmLDP(fac2_SS[1][j],19,3);
		f1<<fmLDP(fac2_SS[2][j],17,3)<<"\n";*/
	}
	f1 << "\n                                          Median Relative Bias\n";
	f1 << "Parameter                 Means                    Medians                  Modes\n";
	for (int j = 0; j < nparscaled; j++) {
		//cout<<nomparam[j]<<"\n";
		f1 << nomparamS[j];
		for (int i = 0; i < 19 - (int)nomparamS[j].length(); i++) f1 << " ";
		f1 << fmLD(rmb_S[0][j], 8, 3) << " " << fmLDP(rmb_SS[0][j], 8, 3);
		f1 << fmLD(rmb_S[1][j], 18, 3) << " " << fmLDP(rmb_SS[1][j], 8, 3);
		f1 << fmLD(rmb_S[2][j], 18, 3) << "\n";
		/*for (int i=0;i<44;i++) f1<<" ";
			f1<<fmLDP(rmb_SS[0][j],8,3);
			f1<<fmLDP(rmb_SS[1][j],15,3);
			f1<<fmLDP(rmb_SS[2][j],15,3)<<"\n";*/
	}
	f1 << "\n                             RMedAD                  Median of the absolute error/true value\n";
	f1 << "Parameter                                          Means                    Medians                  Modes\n";
	for (int j = 0; j < nparscaled; j++) {
		//cout<<nomparam[j]<<"\n";
		f1 << nomparamS[j];
		for (int i = 0; i < 20 - (int)nomparamS[j].length(); i++) f1 << " ";
		f1 << fmLD(rmedad_S[j], 10, 3) << " " << fmLDP(rmedad_SS[j], 8, 3);
		f1 << fmLD(rmae_S[0][j], 15, 3) << " " << fmLDP(rmae_SS[0][j], 8, 3);
		f1 << fmLD(rmae_S[1][j], 15, 3) << " " << fmLDP(rmae_SS[1][j], 8, 3);
		f1 << fmLD(rmae_S[2][j], 19, 3) << "\n";
		/*for (int i=0;i<25;i++) f1<<" ";
		f1<<fmLDP(rmedad_SS[j],11,3);
		f1<<fmLDP(rmae_SS[0][j],16,3);
		f1<<fmLDP(rmae_SS[1][j],16,3);
		f1<<fmLDP(rmae_SS[2][j],16,3)<<"\n";*/
	}
	if (valinfinie) f1 << "\nN.B. Some statisics are noted '+inf' due to zero values of drawn parameters. Check prior minimum values.\n";
	if (valinfinie) cout << "VALEUR INFINIE\n";
	else cout << "PAS DE VALEUR INFINIE\n";
	f1 << "\n";
	valinfinie = false;
	f1.close();

}

/*
 * Traitement des paramètres composites (microsat et séquences uniquement)
 */
void setcompo(int p) {
	long double pmut;
	int kk, qq, k = 0;
	for (int gr = 1; gr < header.ngroupes + 1; gr++) {
		if (groupe[gr].type == 0) {
			if (groupe[gr].priormutmoy.constant) {
				if (groupe[gr].priorsnimoy.constant) {
					pmut = (long double)(groupe[gr].mutmoy + groupe[gr].snimoy);
					for (int j = 0; j < npar; j++) {
						if (scenario[rt.scenteste - 1].histparam[numpar[0][j]].category < 2) {
							enreg2[p].paramvvC[k] = pmut * enreg2[p].paramvv[j];
							if (enreg2[p].paramvvC[k] < zeroplus) enreg2[p].paramvvC[k] = zeroplus;
							k++;
						}
					}
				}
				else {
					kk = 0;
					while (not ((header.mutparam[kk].groupe == gr)and (header.mutparam[kk].category == 2))) kk++;
					for (int j = 0; j < npar; j++) {
						if (scenario[rt.scenteste - 1].histparam[numpar[0][j]].category < 2) {
							pmut = (long double)groupe[gr].mutmoy + (long double)enreg2[p].paramvv[npar + kk];
							enreg2[p].paramvvC[k] = pmut * enreg2[p].paramvv[j];
							if (enreg2[p].paramvvC[k] < zeroplus) enreg2[p].paramvvC[k] = zeroplus;
							k++;
						}
					}
				}
			}
			else {
				if (groupe[gr].priorsnimoy.constant) {
					kk = 0;
					while (not ((header.mutparam[kk].groupe == gr)and (header.mutparam[kk].category == 0))) kk++;
					for (int j = 0; j < npar; j++) {
						if (scenario[rt.scenteste - 1].histparam[numpar[0][j]].category < 2) {
							pmut = (long double)enreg2[p].paramvv[npar + kk] + (long double)groupe[gr].snimoy;
							enreg2[p].paramvvC[k] = pmut * enreg2[p].paramvv[j];
							if (enreg2[p].paramvvC[k] < zeroplus) enreg2[p].paramvvC[k] = zeroplus;
							k++;
						}
					}
				}
				else {
					kk = 0;
					while (not ((header.mutparam[kk].groupe == gr)and (header.mutparam[kk].category == 0))) kk++;
					qq = 0;
					while (not ((header.mutparam[qq].groupe == gr)and (header.mutparam[qq].category == 2))) qq++;
					for (int j = 0; j < npar; j++) {
						if (scenario[rt.scenteste - 1].histparam[numpar[0][j]].category < 2) {
							pmut = (long double)enreg2[p].paramvv[npar + kk] + (long double)enreg2[p].paramvv[npar + qq];
							enreg2[p].paramvvC[k] = pmut * enreg2[p].paramvv[j];
							if (enreg2[p].paramvvC[k] < zeroplus) enreg2[p].paramvvC[k] = zeroplus;
							k++;
						}
					}
				}
			}
		}
		if (groupe[gr].type == 1) {
			if (groupe[gr].priormusmoy.constant) {
				pmut = (long double)groupe[gr].musmoy;
				for (int j = 0; j < npar; j++) {
					if (scenario[rt.scenteste - 1].histparam[numpar[0][j]].category < 2) {
						enreg2[p].paramvvC[k] = pmut * enreg2[p].paramvv[j];
						if (enreg2[p].paramvvC[k] < zeroplus) enreg2[p].paramvvC[k] = zeroplus;
						k++;
					}
				}
			}
			else {
				kk = 0;
				while (not ((header.mutparam[kk].groupe == gr)and (header.mutparam[kk].category == 3))) kk++;
				for (int j = 0; j < npar; j++) {
					if (scenario[rt.scenteste - 1].histparam[numpar[0][j]].category < 2) {
						pmut = enreg2[p].paramvv[npar + kk];
						enreg2[p].paramvvC[k] = pmut * enreg2[p].paramvv[j];
						if (enreg2[p].paramvvC[k] < zeroplus) enreg2[p].paramvvC[k] = zeroplus;
						k++;
					}
				}
			}
		}
	}
}

/**
 * Traitement des paramètres scaled 
 */
void setscaled(int p) {
	long double Ne;
	int k = 0, n, nNe = 0;
	Ne = 0.0;
	for (int j = 0; j < scenario[rt.scenchoisi[0] - 1].npop; j++) {
		for (int ievent = 0; ievent < scenario[rt.scenchoisi[0] - 1].nevent; ievent++) {
			if ((scenario[rt.scenchoisi[0] - 1].event[ievent].action == 'E')and (scenario[rt.scenchoisi[0] - 1].event[ievent].pop == j + 1)) {
				if (scenario[rt.scenchoisi[0] - 1].histparam[j].prior.constant) {
					n = 0;
					while (scenario[rt.scenchoisi[0] - 1].histparam[n].name != scenario[rt.scenchoisi[0] - 1].ne0[j].name) n++;
					Ne += (long double)scenario[rt.scenchoisi[0] - 1].histparam[n].prior.mini;
					//cout<<"j="<<j<< "   Ne="<<scenario[rt.scenchoisi[0]-1].histparam[n].prior.mini<<"\n";
				}
				else {
					//cout<<"  prior.variable"<< scenario[rt.scenchoisi[0]-1].ne0[j].name  <<"\n";
					n = 0;
					while (enreg2[p].name[n].compare(scenario[rt.scenchoisi[0] - 1].ne0[j].name) != 0) n++;
					Ne += (long double)enreg2[p].paramvv[n];
					nNe++;
					//cout<<"j="<<j<<"  "<< "   Ne="<< enreg2[p].paramvv[n] <<"\n";
				}
			}
		}
	}
	//cout<<"Ne="<<Ne<<"\n";
	Ne = Ne / (long double)nNe;
	cout << "Ne=" << Ne << "\n";
	k = 0;
	for (int j = 0; j < npar; j++) {
		if (scenario[rt.scenteste - 1].histparam[numpar[0][j]].category < 2) {
			enreg2[p].paramvvS[k] = enreg2[p].paramvv[j] / Ne;
			if (enreg2[p].paramvvS[k] < zeroplus) enreg2[p].paramvvS[k] = zeroplus;
			k++;
		}
	}
}

/**
 * Interprête la ligne de paramètres de l'option biais et lance les calculs correspondants
 */
void dobias(string opt, int seed) {
	cout << "debut de dobias\n";
	int iprog, nprog, k, nphistarOK, nenr = 10, nenrt = 0, nstatOK;
	int nrec = 0, nsel = 0, nrecpos = 0, ntest = 0, np, ng, npv, ncond, nt, *paordre, *paordreabs;
	string s, s0, s1, sg, entetelog, nomfitrace, nomfipar;
	vector<string> paname;
	float* stat_obs;
	long double** matC;
	bool posterior = false;
	bool firstime = true;
	string bidon;
	long double **phistar, **phistarcompo, **phistarscaled;
	mw.randinit(0, seed);
	progressfilename = path + ident + "_progress.txt";
	scurfile = path + "pseudo-observed_datasets_" + ident + ".txt";
	nomfitrace = path + ident + "_trace.txt";
	nomfipar = path + ident + "_param.txt";
	cout << scurfile << "\n";
	cout << "options : " << opt << "\n";
	vector<string> ss;
	splitwords(opt, ";", ss);
	int ns = ss.size();
	vector<string> ss1;
	for (int i = 0; i < ns; i++) { //cout<<ss[i]<<"\n";
		s0 = ss[i].substr(0, 2);
		s1 = ss[i].substr(2);
		if (s0 == "s:") {
			splitwords(s1, ",", ss1);
			rt.scenteste = atoi(ss1[0].c_str());
			nrecpos = 0;
			nrecpos += rt.nrecscen[rt.scenteste - 1];
			cout << "scenario choisi : " << rt.scenteste << "\n";
		}
		else if (s0 == "n:") {
			nrec = atoi(s1.c_str());
			if (nrec > nrecpos) nrec = nrecpos;
			cout << "nombre total de jeux de données considérés (pour le(s) scénario(s) choisi(s) )= " << nrec << "\n";
		}
		else if (s0 == "m:") {
			nsel = atoi(s1.c_str());
			cout << "nombre de jeux de données considérés pour la régression locale = " << nsel << "\n";
		}
		else if (s0 == "t:") {
			numtransf = atoi(s1.c_str());
			switch (numtransf) {
				case 1: cout << " pas de transformation des paramètres\n";
					break;
				case 2: cout << " transformation log des paramètres\n";
					break;
				case 3: cout << " transformation logit des paramètres\n";
					break;
				case 4: cout << " transformation log(tg) des paramètres\n";
					break;
			}
		}
		else if (s0 == "p:") {
			original = (s1.find("o") != string::npos);
			composite = (s1.find("c") != string::npos);
			scaled = (s1.find("s") != string::npos);
			if ((original)and (not composite)and (not scaled)) cout << "paramètres  originaux  ";
			if ((not original)and (composite)and (not scaled)) cout << "paramètres  composites  ";
			if ((not original)and (not composite)and (scaled)) cout << "paramètres  scaled  ";
			if ((original)and (composite)and (not scaled)) cout << "paramètres  originaux et composites ";
			if ((original)and (not composite)and (scaled)) cout << "paramètres  originaux et scaled ";
			if ((not original)and (composite)and (scaled)) cout << "paramètres  composites et scaled ";
			if ((original)and (composite)and (scaled)) cout << "paramètres  originaux, composites et scaled ";
			cout << "\n";
		}
		else if (s0 == "d:") {
			ntest = atoi(s1.c_str());
			cout << "nombre de jeux-tests à simuler = " << ntest << "\n";
		}
		else if (s0 == "h:") {
			splitwords(s1, " ", ss1);
			np = ss1.size();
			if (np < scenario[rt.scenteste - 1].nparam) {
				//cout<<"le nombre de paramètres transmis ("<<np<<") est incorrect. Le nombre attendu pour le scénario "<<rt.scenteste<<" est de "<<scenario[rt.scenteste-1].nparam<<"\n";
				cout << "the number of parameter transmitted (" << np << ") is incorrect. The expected number for scenario " << rt.scenteste << " is " << scenario[rt.scenteste - 1].nparam << "\n";
				exit(1);
			}
			ncond = np - scenario[rt.scenteste - 1].nparam;
			for (int j = 0; j < scenario[rt.scenteste - 1].nparam; j++) resethistparam(ss1[j]);
			if (ncond > 0) {
				cout << scenario[rt.scenteste - 1].nconditions << "\n";
				if (scenario[rt.scenteste - 1].nconditions != ncond) {
					//if (scenario[rt.scenteste-1].nconditions>0) delete []scenario[rt.scenteste-1].condition;
					scenario[rt.scenteste - 1].condition = vector<ConditionC>(ncond);
				}
				for (int j = 0; j < ncond; j++)
					scenario[rt.scenteste - 1].condition[j].readcondition(ss1[j + scenario[rt.scenteste - 1].nparam]);
			}
		}
		else if (s0 == "u:") {
			cout << s1 << "\n";
			splitwords(s1, "*", ss1);
			ng = ss1.size();
			if (ng != header.ngroupes) {
				cout << "le nombre de groupes transmis (" << ng << ") est incorrect. Le nombre attendu  est de " << header.ngroupes << "\n";
				//exit(1);
			}
			for (int j = 1; j <= ng; j++) resetmutparam(ss1[j - 1]);
		}
		else if (s0 == "po") {
			cout << "paramètres tirés dans les posteriors\n";
			posterior = true;
		}
	}
	rt.nscenchoisi = 1;
	rt.scenchoisi = new int[rt.nscenchoisi];
	rt.scenchoisi[0] = rt.scenteste;
	if (posterior) {
		//calcul des posteriors
		cout << "rt.nrec=" << rt.nrec << "\n";
		nstatOK = rt.cal_varstat();
		cout << "apres cal_varstat  nstatOK = " << nstatOK << "\n";
		cout << "nrec=" << nrec << "     nsel=" << nsel << "\n";
		rt.alloue_enrsel(nsel);
		rt.cal_dist(nrec, nsel, &header.stat_obs[0], false, false);
		cout << "apres cal_dist\n";
		det_numpar();
		cout << "apres det_numpar\n";
		rempli_mat(nsel, &header.stat_obs[0]);
		cout << "apres rempli_mat\n";
		if (not deltanul) matC = cal_matC(nsel);
		recalparamO(nsel);
		cout << "apres recalparam\n";
		if (not deltanul) {
			rempli_parsim(nsel, nparamcom);
			cout << "apres rempli_parsim(O)\n";
			local_regression(nsel, nparamcom, matC);
			cout << "apres local_regression\n";
		}
		phistar = new long double* [nsel];
		for (int i = 0; i < nsel; i++) phistar[i] = new long double[nparamcom];
		if (not deltanul) calphistarO(nsel, phistar);
		else copyphistar(nsel, nparamcom, phistar);
		cout << "apres calphistar\n";
		det_nomparam();
		cout << "apres det_nomparam\n";
		//savephistar(nsel,path,ident,phistar,phistarcompo,phistarscaled);                     cout<<"apres savephistar\n";
		phistarOK = new long double*[nsel];
		for (int i = 0; i < nsel; i++) phistarOK[i] = new long double[rt.nparam[rt.scenteste - 1]];
		cout << "scenario[rt.scenteste-1].nparam = " << scenario[rt.scenteste - 1].nparam << "\n";
		nphistarOK = detphistarOK(nsel, phistar);
		cout << "apres detphistarOK  nphistarOK=" << nphistarOK << "\n";
		cout << "   nphistarOK=" << nphistarOK << "   nstat=" << header.nstat << "\n";
		if (10 * nphistarOK < ntest) {
			cout << "Not enough suitable particles (" << nphistarOK << ")to perform model checking. Stopping computations." << endl;
			exit(1);
		}
		rt.desalloue_enrsel(nsel);
	}
	npv = rt.nparam[rt.scenteste - 1];
	enreg = new enregC[ntest];
	for (int p = 0; p < ntest; p++) {
		enreg[p].stat = vector<float>(header.nstat);
		enreg[p].param = vector<float>(npv);
		enreg[p].numscen = rt.scenteste;
	}
	//fpar.open(nomfipar.c_str());
	if (posterior) {
		cout << "avant dosimulphistar\n";
		while (nenrt < ntest) {
			ps.dosimulphistar(nenr, false, multithread, firstime, rt.scenteste, seed, nphistarOK);
			nenrt += nenr;
			firstime = false;
			if (nenrt > ntest) {
				nenr = ntest - (nenrt - nenr);
				nenrt = ntest;
			}
			cout << "nenrt = " << nenrt << "\r";
			fflush(stdout);
		}
		cout << "apres dosimulphistar\n";
	}
	else {
		cout << "avant dosimultabref\n";
		while (nenrt < ntest) {
			ps.dosimultabref(nenr, false, multithread, firstime, rt.scenteste, seed, 1);
			nenrt += nenr;
			firstime = false;
			if (nenrt > ntest) {
				nenr = ntest - (nenrt - nenr);
				nenrt = ntest;
			}
			cout << "nenrt = " << nenrt << "\r";
			fflush(stdout);
		}
		cout << "\napres dosimultabref\n";
	}
	//fpar.close();
	//header.entete=header.entetehist+header.entetemut0+header.entetestat;
	nprog = 10 * ntest + 6;
	iprog = 5;
	fprog.open(progressfilename.c_str());
	fprog << iprog << "   " << nprog << "\n";
	fprog.close();
	header.readHeader(headerfilename);
	cout << "apres readHeader\n";
	for (int p = 0; p < ntest; p++) {
		enreg[p].param.clear();
		enreg[p].stat.clear();
	}
	delete []enreg;
	det_numpar();
	cout << "naparmcom = " << nparamcom << "   nparcomp = " << nparcompo << "   nparscaled = " << nparscaled << "\n";
	cout << "header.nstat=" << header.nstat << "    rt.nstat=" << rt.nstat << "\n";
	enreg2 = new enreC[ntest];
	for (int p = 0; p < ntest; p++) {
		enreg2[p].stat = new double[header.nstat];
		enreg2[p].paramvv = new double[nparamcom];
		enreg2[p].name = new string[nparamcom];
		if (composite) enreg2[p].paramvvC = new double[nparcompo];
		if (scaled) enreg2[p].paramvvS = new double[nparscaled];
		enreg2[p].numscen = rt.scenchoisi[0];
	}
	cout << "courantfilename=" << scurfile << "\n";
	ifstream file(scurfile.c_str(), ios::in);
	getline(file, entetelog);
	cout << "entetelog=\n" << entetelog << "\n";
	splitwords(entetelog, " ", paname);
	nt = paname.size();
	paordreabs = new int[npv];
	paordre = new int[npv];
	for (int i = 0; i < npv; i++) {
		if (i < rt.nhistparam[rt.scenteste - 1]) {
			k = 0;
			while ((paname[k].compare(rt.histparam[rt.scenteste - 1][i].name) != 0)and (k < nt - 1)) k++;
			paordreabs[i] = k;
			cout << "i=" << i << "    k=" << k << "      ";
			cout << rt.histparam[rt.scenteste - 1][i].name << "   " << paname[k] << "\n";
		}
		else {
			k = 0;
			while ((paname[k].compare(rt.mutparam[i - rt.nhistparam[rt.scenteste - 1]].name) != 0)and (k < nt - 1)) k++;
			paordreabs[i] = k;
			cout << "i=" << i << "    k=" << k << "      ";
			cout << rt.mutparam[i - rt.nhistparam[rt.scenteste - 1]].name << "   " << paname[k] << "\n";
		}
		cout << "paordreabs[" << i << "]=" << paordreabs[i] << "\n";
	}
	for (int i = 0; i < npv; i++) {
		paordre[i] = i;
		if (i < rt.nhistparam[rt.scenteste - 1]) {
			paordre[i] = 0;
			for (int j = 0; j < npv; j++) if (paordreabs[i] > paordreabs[j]) paordre[i]++;
		}
		cout << "paordre[" << i << "]=" << paordre[i] << "\n";
	}
	for (int p = 0; p < ntest; p++) {
		getline(file, bidon);
		cout << "ligne " << p + 1 << " du fichier\n";
		cout << bidon << "\n";
		//cout<<bidon<<"\n";
		splitwords(bidon, " ", ss);
		ns = ss.size();
		cout << "ns=" << ns << "    npv=" << npv << "\n";
		enreg2[p].numscen = atoi(ss[0].c_str());
		cout << "bias.cpp   npv=" << npv << "\n";
		for (int i = 0; i < npv; i++) {
			enreg2[p].paramvv[i] = atof(ss[paordre[i] + 1].c_str());
			if (enreg2[p].paramvv[i] < zeroplus) enreg2[p].paramvv[i] = zeroplus;
			if (i < rt.nhistparam[rt.scenteste - 1]) enreg2[p].name[i] = rt.histparam[rt.scenteste - 1][i].name;
			else enreg2[p].name[i] = rt.mutparam[i - rt.nhistparam[rt.scenteste - 1]].name;
			if (p == 0)cout << "enreg2[p].name[" << i << "]=" << enreg2[p].name[i] << " = " << enreg2[p].paramvv[i] << "\n";
		}
		cout << "avant setcompo\n";
		if (composite) setcompo(p);
		cout << "avant setscaled\n";
		if (scaled) setscaled(p);
		cout << "apres setscaled\n";
		for (int i = 0; i < rt.nstat; i++) {
			cout << "i=" << i << "   i+1+npv=" << i + 1 + npv << "\n";
			cout << "ss[i+1+npv]=" << ss[i + 1 + npv] << "\n";
			enreg2[p].stat[i] = atof(ss[i + 1 + npv].c_str());

		}
		cout << "apres copie des stat\n";
	}
	file.close();
	nstatOK = rt.cal_varstat();
	cout << "apres cal_varstat  nstatOK = " << nstatOK << "\n";
	stat_obs = new float[rt.nstat];
	rt.alloue_enrsel(nsel);
	cout << "apres rt.alloue_enrsel\n";
	if (original) {
		initbiasO(ntest, nsel);
		initbiasOS(ntest, nsel);
		phistar = new long double* [nsel];
		for (int i = 0; i < nsel; i++) phistar[i] = new long double[nparamcom];
		lisimparO();
	}
	if (composite) {
		initbiasC(ntest, nsel); //PP: was nparcompo
		initbiasCS(ntest, nsel); //PP: was nparcompo
		phistarcompo = new long double* [nsel];
		for (int i = 0; i < nsel; i++) phistarcompo[i] = new long double[nparcompo];
		lisimparC();
	}
	if (scaled) {
		initbiasS(ntest, nsel); //PP: was nparcompo
		initbiasSS(ntest, nsel); //PP: was nparcompo
		phistarscaled = new long double* [nsel];
		for (int i = 0; i < nsel; i++) phistarscaled[i] = new long double[nparscaled];
		lisimparS();
	}
	ftrace.open(nomfitrace.c_str());
	numero = new int[nsel];
	for (int p = 0; p < ntest; p++) {
		cout << "\nanalysing data test " << p + 1 << " \n";
		for (int j = 0; j < rt.nstat; j++) stat_obs[j] = enreg2[p].stat[j];
		rt.cal_dist(nrec, nsel, stat_obs, false, false);
		if (debuglevel == 11) cout << "apres cal_dist\n";
		iprog += 6;
		fprog.open(progressfilename.c_str(), ios::out);
		fprog << iprog << "   " << nprog << "\n";
		fprog.close();
		if (p < 1) det_nomparam();
		rempli_mat(nsel, stat_obs);
		if (debuglevel == 11) cout << "apres rempli_mat\n";
		if (not deltanul) matC = cal_matC(nsel);
		if (debuglevel == 11) cout << "cal_matC\n";

		for (int kk = 0; kk < nsel; kk++) numero[kk] = mw.rand0(nsimpar);

		if (original) {
			recalparamO(nsel);
			if (debuglevel == 11) cout << "apres recalparamO\n";
			if (not deltanul) {
				rempli_parsim(nsel, nparamcom);
				if (debuglevel == 11) cout << "apres rempli_mat\n";
				local_regression(nsel, nparamcom, matC);
				if (debuglevel == 11) cout << "apres local_regression\n";
			}
			if (not deltanul) calphistarO(nsel, phistar);
			else copyphistar(nsel, nparamcom, phistar);
			if (debuglevel == 11) cout << "apres calphistarO\n";
			paramest[p] = calparstat(nsel, nparamcom, phistar);
			if (debuglevel == 11) cout << "apres calparstatO\n";
			for (int i = 0; i < nsel; i++) {
				for (int j = 0; j < nparamcom; j++) paretoil[i][j] = phistar[i][j];
				//for (int j=0;j<nparamcom;j++) cout<<"  "<<phistar[i][j]<<" ("<<enreg2[p].paramvv[j] <<")";
				//cout<<"\n";
			}
			biaisrelO(ntest, nsel, p);
			if (debuglevel == 11) cout << "apres biaisrelO\n";
			tracebiais(p);
		}
		if (composite) {
			recalparamC(nsel);
			if (debuglevel == 11) cout << "apres recalparamC\n";
			if (not deltanul) {
				rempli_parsim(nsel, nparcompo);
				if (debuglevel == 11) cout << "apres rempli_mat\n";
				local_regression(nsel, nparcompo, matC);
				if (debuglevel == 11) cout << "apres local_regression\n";
			}
			//phistarcompo = new long double*[nsel];
			//for (int i=0;i<nsel;i++) phistarcompo[i] = new long double[nparcompo];
			if (not deltanul) calphistarC(nsel, phistarcompo);
			else copyphistar(nsel, nparcompo, phistarcompo);
			if (debuglevel == 11) cout << "apres calphistarC\n";
			paramestcompo[p] = calparstat(nsel, nparcompo, phistarcompo);
			if (debuglevel == 11) cout << "apres calparstatC\n";
			for (int i = 0; i < nsel; i++) {
				for (int j = 0; j < nparcompo; j++) paretoilcompo[i][j] = phistarcompo[i][j];
			}
			//cout<<"avant delete phistarcompo\n";
			//for (int i=0;i<nsel;i++) delete []phistarcompo[i];delete phistarcompo;
			//cout<<"apres delete phistarcompo\n";
			if (debuglevel == 11) cout << "avant biaisrelC\n";
			biaisrelC(ntest, nsel, p);
			if (debuglevel == 11) cout << "après biaisrelC\n";

		}
		if (scaled) {
			recalparamS(nsel);
			if (debuglevel == 11) cout << "apres recalparamS\n";
			if (not deltanul) {
				rempli_parsim(nsel, nparscaled);
				if (debuglevel == 11) cout << "apres rempli_mat\n";
				local_regression(nsel, nparscaled, matC);
				if (debuglevel == 11) cout << "apres local_regression\n";
			}
			if (not deltanul) calphistarS(nsel, phistarscaled);
			else copyphistar(nsel, nparscaled, phistarcompo);
			if (debuglevel == 11) cout << "apres calphistarS\n";
			paramestscaled[p] = calparstat(nsel, nparscaled, phistarscaled);
			if (debuglevel == 11) cout << "apres calparstatS\n";
			for (int i = 0; i < nsel; i++) {
				for (int j = 0; j < nparscaled; j++) paretoilscaled[i][j] = phistarscaled[i][j];
				//for (int j=0;j<nparscaled;j++) cout<<"  "<<phistarscaled[i][j]<<" ("<<enreg2[p].paramvvS[j] <<")";
				//cout<<"\n";
			}
			biaisrelS(ntest, nsel, p);
			if (debuglevel == 11) cout << "apres biaisrelS\n";
		}
		if (not deltanul) {
			for (int i = 0; i < nstatOKsel + 1; i++) delete [] matC[i];
			delete [] matC;
			delete_mat(nsel);
		}
		iprog += 4;
		fprog.open(progressfilename.c_str());
		fprog << iprog << "   " << nprog << "\n";
		fprog.close();
	}
	ftrace.close();
	rt.desalloue_enrsel(nsel);
	//		for (int i=0;i<nsel;i++) delete []phistar[i];delete phistar;
	if (original) {
		for (int i = 0; i < nsimpar; i++) delete []simpar[i];
		delete [] simpar;
		finbiaisrelO(ntest);
		ecriresO(ntest, nsel);
	}
	if (composite) {
		finbiaisrelC(ntest);
		ecriresC(ntest, nsel);
	}
	if (scaled) {
		finbiaisrelS(ntest);
		if (debuglevel == 11) cout << "apres biaisrelS\n";
		ecriresS(ntest, nsel);
		if (debuglevel == 11) cout << "apres ecriresS\n";
	}
	iprog += 1;
	fprog.open(progressfilename.c_str());
	fprog << iprog << "   " << nprog << "\n";
	fprog.close();
}
