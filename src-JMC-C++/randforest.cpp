/*
* randforest.cpp
*
* Created on : 24 november 2014
*
*/

#include <string>
#include <vector>
#include <fstream>
#include <iostream>
#include <algorithm>
#include <iomanip>
#include <functional>
#include <stdlib.h>
#include <time.h>
#include <math.h>

#include "randomgenerator.hpp"
#include "particleset.hpp"
#include "mesutils.hpp"
#include "header.hpp"
#include "reftable.hpp"
#include "randforest.hpp"

#define LOG2 0.69314718055994530942;


using namespace std;
using namespace std::placeholders;

extern string ident, headerfilename;
extern bool multithread;
extern string progressfilename, path, ident;
extern string scurfile;
extern HeaderC header;
extern ReftableC rt;
extern int debuglevel;
extern enregC* enreg;
extern long double** phistarOK;
extern ParticleSetC ps;

RFC rf;
RFC* rfmin;

ofstream fout;
int nvoisins, nj, ndone, ntestes, njustes, nlim = 5;
double seuil = 0.01;
string nomrfmin = "foretmin.bin";
vector<string> nomstat;
vector<int> obs_estscen;
vector<vector<int>> sim_estscen;
resAFD afd;
bool LD = false, flagk = false, flago = false;
TreeC ctree;
vector<vector<int>> nimportance;
vector<vector<double>> importance;
vector<double> obs_estim;

bool operator <(const VMC& lhs, const VMC& rhs) {
	return lhs.x < rhs.x;
}

bool operator <(const VMD& lhs, const VMD& rhs) {
	return lhs.x > rhs.x;
}

/**
* Definition de l'operateur = pour une instance de la classe NodeC
*/
NodeRC& NodeRC::operator=(NodeRC const& source) {
	if (this == &source) return *this;
	int imax;
	this->pere = source.pere;
	this->filsG = source.filsG;
	this->filsD = source.filsD;
	this->nvar = source.nvar;
	this->nsets = source.nsets;
	this->npassages = source.npassages;
	this->nsetG = source.nsetG;
	this->nsetD = source.nsetD;
	this->model = source.model;
	this->imax = source.imax;
	this->cutval = source.cutval;
	this->disval = source.disval;
	this->delta = source.delta;
	this->terminal = source.terminal;
	if (not this->indvar.empty()) this->indvar.clear();
	if (not source.indvar.empty()) {
		this->indvar = vector<int>(source.indvar.size());
		imax = (int)source.indvar.size();
		for (int i = 0; i < imax; i++) this->indvar[i] = source.indvar[i];
	}
	if (not this->numset.empty()) this->numset.clear();
	if (not source.numset.empty()) {
		this->numset = vector<int>(source.numset.size());
		imax = (int)source.numset.size();
		for (int i = 0; i < imax; i++) this->numset[i] = source.numset[i];
	}
	if (not this->numsetG.empty()) this->numsetG.clear();
	if (not source.numsetG.empty()) {
		this->numsetG = vector<int>(source.numsetG.size());
		imax = (int)source.numsetG.size();
		for (int i = 0; i < imax; i++) this->numsetG[i] = source.numsetG[i];
	}
	if (not this->numsetD.empty()) this->numsetD.clear();
	if (not source.numsetD.empty()) {
		this->numsetD = vector<int>(source.numsetD.size());
		imax = (int)source.numsetD.size();
		for (int i = 0; i < imax; i++) this->numsetD[i] = source.numsetD[i];
	}
	return *this;
}


/**
* Definition de l'operateur = pour une instance de la classe TreeC
*/
TreeC& TreeC::operator=(TreeC const& source) {
	if (this == &source) return *this;
	int imax;
	this->nnodes = source.nnodes;
	this->nsets = source.nsets;
	this->nvar = source.nvar;
	this->fin = source.fin;
	/// ON NE COPIE PAS LE GENERATEUR DE NOMBRES ALEATOIRES 
	if (not this->numset.empty()) this->numset.clear();
	if (not source.numset.empty()) {
		this->numset = vector<int>(source.numset.size());
		imax = (int)source.numset.size();
		for (int i = 0; i < imax; i++) this->numset[i] = source.numset[i];
	}
	/*if (not this->indsel.empty()) this->indsel.clear();
	if (not source.indsel.empty()) {
	this->indsel = vector<int>(source.indsel.size());
	imax=(int)source.indsel.size();
	for (int i=0;i<imax;i++) this->indsel[i] = source.indsel[i];
	}*/
	if (not this->score.empty()) this->score.clear();
	if (not source.score.empty()) {
		this->score = vector<int>(source.score.size());
		imax = (int)source.score.size();
		for (int i = 0; i < imax; i++) this->score[i] = source.score[i];
	}
	if (not this->index.empty()) this->index.clear();
	if (not source.index.empty()) {
		this->index = vector<int>(source.index.size());
		imax = (int)source.index.size();
		for (int i = 0; i < imax; i++) this->index[i] = source.index[i];
	}
	if (not this->node.empty()) this->node.clear();
	if (not source.node.empty()) {
		this->node = vector<NodeRC>(source.node.size());
		imax = (int)source.node.size();
		for (int i = 0; i < imax; i++) this->node[i] = source.node[i];
	}
	if (not this->varused.empty()) this->varused.clear();
	if (not source.varused.empty()) {
		this->varused = vector<bool>(source.varused.size());
		imax = (int)source.varused.size();
		for (int i = 0; i < imax; i++) this->varused[i] = source.varused[i];
	}
	return *this;
}

/**
* Definition de l'operateur = pour une instance de la classe RFC
*/
RFC& RFC::operator=(RFC const& source) {
	if (this == &source) return *this;
	int imax, jmax;
	this->ntrees = source.ntrees;
	this->ntot = source.ntot;
	this->nsets = source.nsets;
	this->nstat = source.nstat;
	this->nmodel = source.nmodel;
	this->nvar = source.nvar;
	this->nsel = source.nsel;
	this->nstatclass = source.nstatclass;
	this->nbootsamp = source.nbootsamp;
	if (not this->model.empty()) this->model.clear();
	if (not source.model.empty()) {
		this->model = vector<int>(source.model.size());
		imax = (int)source.model.size();
		for (int i = 0; i < imax; i++) this->model[i] = source.model[i];
	}
	if (not this->vote.empty()) {
		imax = (int)this->vote.size();
		for (int i = 0; i < imax; i++) {
			if (not this->vote[i].empty()) this->vote[i].clear();
		}
		this->vote.clear();
	}
	if (not source.vote.empty()) {
		this->vote = vector<vector<double>>(source.vote.size());
		imax = (int)source.vote.size();
		for (int i = 0; i < imax; i++) {
			if (not source.vote[i].empty()) {
				this->vote[i] = vector<double>(source.vote[i].size());
				jmax = (int)source.vote[i].size();
				for (int j = 0; j < jmax; j++) this->vote[i][j] = source.vote[i][j];
			}
		}
	}
	if (not this->stat.empty()) {
		imax = (int)this->stat.size();
		for (int i = 0; i < imax; i++) {
			if (not this->stat[i].empty()) this->stat[i].clear();
		}
		this->stat.clear();
	}
	if (not source.stat.empty()) {
		this->stat = vector<vector<double>>(source.stat.size());
		imax = (int)source.stat.size();
		for (int i = 0; i < imax; i++) {
			if (not source.stat[i].empty()) {
				this->stat[i] = vector<double>(source.stat[i].size());
				jmax = (int)source.stat[i].size();
				for (int j = 0; j < jmax; j++) this->stat[i][j] = source.stat[i][j];
			}
		}
	}
	if (not this->bootsamp.empty()) this->bootsamp.clear();
	if (not source.bootsamp.empty()) {
		this->bootsamp = vector<int>(source.bootsamp.size());
		imax = (int)source.bootsamp.size();
		for (int i = 0; i < imax; i++) this->bootsamp[i] = source.bootsamp[i];
	}
	if (not this->varimp.empty()) this->varimp.clear();
	if (not source.varimp.empty()) {
		this->varimp = vector<double>(source.varimp.size());
		imax = (int)source.varimp.size();
		for (int i = 0; i < imax; i++) this->varimp[i] = source.varimp[i];
	}
	if (not this->tree.empty()) this->tree.clear();
	if (not source.tree.empty()) {
		this->tree = vector<TreeC>(source.tree.size());
		imax = (int)source.tree.size();
		for (int i = 0; i < imax; i++) this->tree[i] = source.tree[i];
	}
	if (not this->statobs.empty()) this->statobs.clear();
	if (not source.statobs.empty()) {
		this->statobs = vector<double>(source.statobs.size());
		imax = (int)source.statobs.size();
		for (int i = 0; i < imax; i++) this->statobs[i] = source.statobs[i];
	}
	if (not this->statname.empty()) this->statname.clear();
	if (not source.statname.empty()) {
		this->statname = vector<string>(source.statname.size());
		imax = (int)source.statname.size();
		for (int i = 0; i < imax; i++) this->statname[i] = source.statname[i];
	}
	return *this;
}

int NodeRC::regle3(vector<VMC> vm, vector<int> b, MwcGen& mw) {
	bool ident;
	int modmax, nident, bmax;
	vector<double> xmin, xmax;
	xmin = vector<double>(this->nvar);
	xmax = vector<double>(this->nvar);
	for (int i = 0; i < this->nvar; i++) {
		for (int j = 0; j < nsets; j++) vm[j].ind = rf.model[numset[j]];
		vm[0].x = rf.stat[numset[0]][indvar[i]];
		xmax[i] = xmin[i] = vm[0].x;
		for (int j = 1; j < nsets; j++) {
			vm[j].x = rf.stat[numset[j]][indvar[i]];
			if (vm[j].x > xmax[i]) xmax[i] = vm[j].x;
			if (vm[j].x < xmin[i]) xmin[i] = vm[j].x;
		}
		ident = (xmax[i] - xmin[i]) < 1E-16 * (xmax[i] + xmin[i]);
		if (not ident) {
			xmin.clear();
			xmax.clear();
			return -1; //il existe au moins deux individus aux stat différentes
		}
	}
	modmax = 0;
	bmax = b[0];
	for (int j = 1; j < rf.nmodel; j++)
		if (b[j] > bmax) {
			bmax = b[j];
			modmax = j;
		}
	nident = 0;
	for (int j = 0; j < rf.nmodel; j++) if (b[j] == bmax) nident++;
	if (nident < 2) { //
		xmin.clear();
		xmax.clear();
		return modmax; //tous les individus ont les mêmes valeurs de stat et le modèle modmax est majoritaire
	}
	else {
		int kk = mw.rand1(nident);
		int kkk = 0;
		for (int j = 0; j < rf.nmodel; j++) {
			if (b[j] == bmax) {
				kkk++;
				if (kkk == kk) modmax = j;
			}
		}
		xmin.clear();
		xmax.clear();
		return modmax;
	}
}

double NodeRC::calGinimin(vector<VMC>& vm, double& cutvalc, vector<int> modfreq) const {
	//cout<<"calGinimin 0\n";
	double epsilon, gini, res, cg, cd;
	int gg = 0, dd = this->nsets, j = 0, imax;
	vector<int> ga, dr;
	ga = vector<int>(rf.nmodel, 0);
	dr = vector<int>(rf.nmodel);
	for (int i = 0; i < rf.nmodel; i++) dr[i] = modfreq[i];
	vector<double> v;
	v.resize(0);
	//cout<<"calGinimin 1\n";
	epsilon = vm[this->nsets / 2].x * 1E-10;
	for (int i = 1; i < this->nsets; i++) {
		if (vm[i].x - vm[i - 1].x > epsilon) v.push_back(0.5 * (vm[i].x + vm[i - 1].x));
	}
	//cout<<"calGinimin v.size="<<v.size()<<"\n";
	gini = 1.0;
	imax = (int)v.size();
	for (int i = 0; i < imax; i++) {
		while (vm[j].x < v[i]) {
			ga[vm[j].ind]++;
			gg++;
			dr[vm[j].ind]--;
			dd--;
			j++;
		}
		//cout<<"calGinimin j="<<j<<"\n";
		res = 0.0;
		for (int m = 0; m < rf.nmodel; m++) {
			if (gg > 0) cg = (double)ga[m] / (double)gg; else cg = 0.0;
			if (dd > 0) cd = (double)dr[m] / (double)dd; else cd = 0.0;
			res += cg * (1.0 - cg) * (double)gg + cd * (1.0 - cd) * (double)dd;
		}
		res = res / (double)(gg + dd);
		//cout<<"calGinimin res="<<res<<"\n";
		if (res < gini) {
			gini = res;
			cutvalc = v[i];
		}
		//cout<<"calGinimin gini="<<gini<<"\n";
	}
	ga.clear();
	dr.clear();
	v.clear();
	return gini;
}

double NodeRC::calvarmin(vector<VMC>& vm, double& cutvalc) const {
	//cout<<"calGinimin 0\n";
	double epsilon, vv, va, vg, vd, sx2g, sx2d, sxg, sxd, d2, d;
	int ng, nd, j = 0, imax;
	vector<double> v;
	v.resize(0);
	epsilon = vm[this->nsets / 2].x * 1E-10;
	for (int i = 1; i < this->nsets; i++) {
		if (vm[i].x - vm[i - 1].x > epsilon) v.push_back(0.5 * (vm[i].x + vm[i - 1].x));
	}
	vv = 10.0;
	imax = (int)v.size();
	sx2g = 0.0;
	sxg = 0.0;
	ng = 0;
	sx2d = 0.0;
	sxd = 0.0;
	nd = this->nsets;
	for (int i = 1; i < this->nsets; i++) {
		sx2d += (double)(vm[i].ind * vm[i].ind);
		sxd += (double)vm[i].ind;
	}
	for (int i = 0; i < imax; i++) {
		//cout<<"v["<<i<<"]="<<v[i]<<"\n";
		//cout<<"ng="<<ng<<"  sx2g="<<sx2g<<"  sxg="<<sxg<<"    nd="<<nd<<"  sx2d="<<sx2d<<"  sxd="<<sxd<<"\n";
		//                      d2=0.0;d=0.0;
		while (vm[j].x < v[i]) {
			d = (double)(vm[j].ind);
			d2 = (double)(vm[j].ind * vm[j].ind);
			sx2d -= d2;
			sxd -= d;
			sx2g += d2;
			sxg += d;
			ng++;
			nd--;
			j++;
		}
		vg = sx2g - sxg * sxg / (double)ng;
		vd = sx2d - sxd * sxd / (double)nd;
		va = (vg + vd) / (double)(this->nsets);
		//cout<<"va="<<va<<"     i="<<i<<"    j="<<j<<"\n";
		if (va < vv) {
			vv = va;
			cutvalc = v[i];
		}
		//cout<<"   calvarmin="<<vv<<"\n";
	}
	v.clear();//exit(1);
	return vv;
}

double NodeRC::calGini(vector<VMC>& vm, double cutval) const {
	double cg, cd, gini = 0.0;
	int gg = 0, dd = 0;
	vector<int> ga, dr;
	ga = vector<int>(rf.nmodel, 0);
	dr = vector<int>(rf.nmodel, 0);
	gini = 1.0;
	for (int i = 0; i < this->nsets; i++) {
		if (vm[i].x < cutval) {
			ga[vm[i].ind]++;
			gg++;
		}
		else {
			dr[vm[i].ind]++;
			dd++;
		}
	}
	for (int m = 0; m < rf.nmodel; m++) {
		if (gg > 0) cg = (double)ga[m] / (double)gg; else cg = 0.0;
		if (dd > 0) cd = (double)dr[m] / (double)dd; else cd = 0.0;
		gini += cg * (1.0 - cg) * (double)gg + cd * (1.0 - cd) * (double)dd;
	}
	gini /= (double)(gg + dd);
	ga.clear();
	dr.clear();
	return gini;
}

double NodeRC::calvarmoy(vector<VMC>& vm, double cutval) const {
	double sx2g = 0, sxg = 0, sx2d = 0, sxd = 0, vg, vd = 0.0, v;
	int ng = 0, nd = 0;
	for (int i = 0; i < this->nsets; i++) {
		if (vm[i].x < cutval) {
			ng++;
			sx2g += (double)(vm[i].ind * vm[i].ind);
			sxg += (double)(vm[i].ind);
		}
		else {
			nd++;
			sx2d += (double)(vm[i].ind * vm[i].ind);
			sxd += (double)(vm[i].ind);
		}
	}
	vg = (sx2g - sxg * sxg / (double)ng);
	vd = (sx2d - sxd * sxd / (double)nd);
	v = (vg + vd) / (double)this->nsets;
	// cout << "      calvarmoy=" << v << "   vg=" << vg << "  vd=" << vd << "    nsets=" << this->nsets << "\n";
	return v;
}


/**
*  teste si le noeud est terminal :
*  retourne -1 si le noeud n'est pas terminal
*  retourne le numéro du modèle si le noeud est terminal
*/
int NodeRC::getdisval(MwcGen& mw) {
	double gini, ginimin, cutvalmin, dfa, c, deltamax, cutvalc = 0.0;
	int modmax, freqmax, ii;
	vector<VMC> vm;
	vm = vector<VMC>(this->nsets);
	for (int j = 0; j < this->nsets; j++) vm[j].ind = rf.model[this->numset[j]];
	//calcul du nombre d'individus de chaque modèle
	vector<int> modfreq;
	modfreq = vector<int>(rf.nmodel, 0);
	for (int j = 0; j < this->nsets; j++) modfreq[vm[j].ind]++;
	modmax = 0;
	for (int k = 1; k < rf.nmodel; k++) if (modfreq[k] > modfreq[modmax]) modmax = k;
	if (modfreq[modmax] == this->nsets) { //tous les individus sont du même modèle
		vm.clear();
		modfreq.clear();
		return modmax;
	}
	modmax = regle3(vm, modfreq, mw);
	if (modmax != -1) { //tous les individus ont les mêmes valeurs de stat
		vm.clear();
		modfreq.clear();
		return modmax;
	}
	//calcul du Gini du noeud avant split
	dfa = 0.0;
	for (int m = 0; m < rf.nmodel; m++) {
		c = (double)modfreq[m] / (double)this->nsets;
		dfa += c * (1.0 - c);
	}
	//calcul du nombre d'individus du modèle le plus fréquent à ce noeud
	freqmax = modfreq[0];
	for (int k = 1; k < rf.nmodel; k++) if (modfreq[k] > freqmax) freqmax = modfreq[k];
	if (freqmax == 1) { //il y a au plus 1 seul individu par modèle
		do {
			ii = mw.rand0(this->nvar);
			for (int j = 0; j < this->nsets; j++) {
				vm[j].x = rf.stat[this->numset[j]][indvar[ii]];
				vm[j].ind = rf.model[this->numset[j]];
			}
			sort(&vm[0], &vm[0] + nsets);
			deltamax = fabs(vm[this->nsets - 1].x - vm[0].x);
		} while (deltamax < fabs(vm[this->nsets / 2].x * 1E-10));
		cutvalmin = vm[this->nsets / 2].x;
		this->imax = indvar[ii];
		for (int j = 0; j < this->nsets; j++) {
			vm[j].x = rf.stat[this->numset[j]][this->imax];
			vm[j].ind = rf.model[this->numset[j]];
		}
		ginimin = calGini(vm, cutvalmin);
	}
	else { //il y a plusieurs individus d'un modèle donné
		for (int i = 0; i < this->nvar; i++) {
			for (int j = 0; j < this->nsets; j++) {
				vm[j].x = rf.stat[this->numset[j]][indvar[i]];
				vm[j].ind = rf.model[this->numset[j]];
			}
			sort(&vm[0], &vm[0] + nsets);
			gini = calGinimin(vm, cutvalc, modfreq);
			if (i == 0) {
				this->imax = indvar[i];
				cutvalmin = cutvalc;
				ginimin = gini;
			}
			else {
				if (gini < ginimin) {
					this->imax = indvar[i];
					cutvalmin = cutvalc;
					ginimin = gini;
				}
			}
		}
		for (int j = 0; j < this->nsets; j++) {
			vm[j].x = rf.stat[this->numset[j]][imax];
			vm[j].ind = rf.model[this->numset[j]];
		}
	}
	this->delta = (dfa - ginimin) * this->nsets;
	this->disval = ginimin;
	this->cutval = cutvalmin;
	this->numsetG.resize(0);
	this->numsetD.resize(0);
	for (int j = 0; j < this->nsets; j++) {
		if (vm[j].x < this->cutval) this->numsetG.push_back(this->numset[j]);
		else this->numsetD.push_back(this->numset[j]);
	}
	this->nsetG = this->numsetG.size();
	this->nsetD = this->numsetD.size();
	vm.clear();
	modfreq.clear();
	return -1;
}

/**
*  teste si le noeud est terminal :
*  retourne -1 si le noeud n'est pas terminal
*  retourne le numéro du modèle si le noeud est terminal
*/
double NodeRC::getdisval2(MwcGen& mw) {
	bool terminal = false;
	double va, vamin, cutvalmin, deltamax, cutvalc = 0.0;
	int modmax, freqmax, ii;
	double freqmoy = 0.0;
	vector<VMC> vm = vector<VMC>(this->nsets);
	for (int j = 0; j < this->nsets; j++) {
		vm[j].ind = rf.model[this->numset[j]];
		freqmoy += (double)rf.model[this->numset[j]];
	}
	freqmoy /= (double)this->nsets;
	//calcul du nombre d'individus de chaque modèle
	vector<int> modfreq = vector<int>(rf.nmodel, 0);
	for (int j = 0; j < this->nsets; j++) modfreq[vm[j].ind]++;
	modmax = 0;
	for (int k = 1; k < rf.nmodel; k++) if (modfreq[k] > modfreq[modmax]) modmax = k;
	if (this->nsets <= nlim) { //la limite de 5 datasets est atteinte
		vm.clear();
		modfreq.clear();
		return freqmoy;
	}
	if (modfreq[modmax] == this->nsets) { //tous les datasets sont du même modèle
		vm.clear();
		modfreq.clear();
		return freqmoy;
	}
	//  modmax = regle3(vm, modfreq, mw);
	//  if (modmax != -1) { //tous les datasets ont les mêmes valeurs de stat
	//      vm.clear();
	//cout << "Erreur regle 3 activee en regression, impossible" << endl;
	//exit(1);
	//      modfreq.clear();
	//      return freqmoy;
	//  }
	//calcul du Gini du noeud avant split
	//dfa=0.0;for (int m=0;m<rf.nmodel;m++) {c=(double)modfreq[m]/(double)this->nsets;dfa +=c*(1.0-c);}
	//calcul du nombre d'individus du modèle le plus fréquent à ce noeud
	freqmax = modfreq[0];
	for (int k = 1; k < rf.nmodel; k++) if (modfreq[k] > freqmax) freqmax = modfreq[k];
	if (freqmax == 1) { //il y a au plus 1 seul individu par modèle
		do {
			ii = mw.rand0(this->nvar);
			for (int j = 0; j < this->nsets; j++) {
				vm[j].x = rf.stat[this->numset[j]][indvar[ii]];
				vm[j].ind = rf.model[this->numset[j]];
			}
			sort(&vm[0], &vm[0] + nsets);
			deltamax = fabs(vm[this->nsets - 1].x - vm[0].x);
		} while (deltamax < fabs(vm[this->nsets / 2].x * 1E-10));
		cutvalmin = vm[this->nsets / 2].x;
		//ginimin=calGini(vm,cutvalmin);
		this->imax = indvar[ii];
		for (int j = 0; j < this->nsets; j++) {
			vm[j].x = rf.stat[this->numset[j]][this->imax];
			vm[j].ind = rf.model[this->numset[j]];
		}
		vamin = calvarmoy(vm, cutvalmin);
	}
	else { //il y a plusieurs individus d'un modèle donné
		double ns = static_cast<double>(this->nsets);
		if (6.15*ns*log(ns) < static_cast<double>(rf.nsets)) {
			//if (true) {
			for (int i = 0; i < this->nvar; i++) {
				for (int j = 0; j < this->nsets; j++) {
					vm[j].x = rf.stat[this->numset[j]][indvar[i]];
					vm[j].ind = rf.model[this->numset[j]];
				}
				sort(&vm[0], &vm[0] + nsets);
				va = calvarmin(vm, cutvalc);
				if (i == 0) {
					this->imax = indvar[i];
					cutvalmin = cutvalc;
					vamin = va;
				}
				else {
					if (va < vamin) {
						this->imax = indvar[i];
						cutvalmin = cutvalc;
						vamin = va;
					}
				}
			}
			for (int j = 0; j < this->nsets; j++) {
				vm[j].x = rf.stat[this->numset[j]][imax];
				vm[j].ind = rf.model[this->numset[j]];
			}
		}
		else {
			vector<int> isIn = vector<int>(rf.nsets, 0);
			for (int i = 0; i < this->nvar; i++) {
				for (int j = 0; j < this->nsets; j++) {
					isIn[this->numset[j]]++;
				}
				int j = 0;
				int realind;
				int current_var_ind = indvar[i];
				vector<int>& current_stat_vec = rf.stat_sorted_ind[current_var_ind];
				for (int jj = 0; jj < rf.nsets; jj++) {
					realind = current_stat_vec[jj];
					for (int jji = 0; jji < isIn[realind]; jji++) {
						vm[j].x = rf.stat[realind][current_var_ind];
						vm[j].ind = rf.model[realind];
						j++;
					}
				}
				va = calvarmin(vm, cutvalc);
				if (i == 0) {
					imax = indvar[i];
					cutvalmin = cutvalc;
					vamin = va;
				}
				else {
					if (va < vamin) {
						imax = indvar[i];
						cutvalmin = cutvalc;
						vamin = va;
					}
				}
				for (int jj = 0; jj < this->nsets; jj++) {
					isIn[this->numset[jj]] = 0;
				}

			}
			for (int j = 0; j < this->nsets; j++) {
				vm[j].x = rf.stat[this->numset[j]][imax];
				vm[j].ind = rf.model[this->numset[j]];
			}
		}
	}
	//this->delta=(dfa-ginimin)*this->nsets;        
	this->disval = vamin;
	this->cutval = cutvalmin;
	this->numsetG.resize(0);
	this->numsetD.resize(0);
	for (int j = 0; j < this->nsets; j++) {
		if (vm[j].x < this->cutval) this->numsetG.push_back(this->numset[j]);
		else this->numsetD.push_back(this->numset[j]);
	}
	this->nsetG = this->numsetG.size();
	this->nsetD = this->numsetD.size();
	double result = -1;
	vm.clear();
	modfreq.clear();
	//    if ((this->nsetG < 5)or (this->nsetD < 5)) return freqmoy;
	return -1;
}

int TreeC::infermodel(const vector<double>& stat) {
	int k = 0;
	do {
		if (stat[node[k].imax] < node[k].cutval) k = node[k].filsG; else k = node[k].filsD;
	} while (not node[k].terminal);
	return node[k].model;
}

void TreeC::initree() {
	this->nsets = rf.nsel;
	this->nvar = rf.nvar;
	this->nnodes = 2 * this->nsets - 1;
	this->node = vector<NodeRC>(2 * this->nsets);
	this->numset = vector<int>(this->nsets);
	//this->indsel = vector <int>(this->nsets);
	this->score = vector<int>(rf.nstat);
	this->mw.samplewith(rf.nsets, this->nsets, this->numset/*this->indsel*/);
	this->sim_participe = vector<bool>(rf.nsets);
	for (int i = 0; i < rf.nsets; i++) this->sim_participe[i] = false;
	for (int i = 0; i < this->nsets; i++) this->sim_participe[this->numset[i]] = true;
	//for (int j=0;j<this->nsets;j++) this->numset[j] = this->indsel[j];
	this->varused = vector<bool>(rf.nstat, false);
}

void TreeC::deletree() {
	if (not this->numset.empty()) this->numset.clear();
	if (not this->sim_participe.empty()) this->sim_participe.clear();
	if (not this->score.empty()) this->score.clear();
	if (not this->index.empty()) this->index.clear();
	if (not this->varused.empty()) this->varused.clear();
	if (not this->node.empty()) {
		for (int m = 0; m < this->nnodes; m++) {
			if (not this->node[m].indvar.empty()) this->node[m].indvar.clear();
			if (not this->node[m].numset.empty()) this->node[m].numset.clear();
			if (not this->node[m].numsetG.empty()) this->node[m].numsetG.clear();
			if (not this->node[m].numsetD.empty()) this->node[m].numsetD.clear();
		}
		this->node.clear();
	}
}


void TreeC::estim() {
	//cout<<"debut ESTIM\n";
	obs_estscen[this->infermodel(rf.statobs)]++;
	//cout<<"avant la boucle\n";
	for (int j = 0; j < rf.nsets; j++) {
		//cout<<"this->sim_participe[j]="<<this->sim_participe[j]<<"\n";
		if (not this->sim_participe[j]) {
			sim_estscen[j][this->infermodel(rf.stat[j])]++;
		}
	}
	//cout<<"fin ESTIM\n";
}

double TreeC::inferobs(vector<double>& stat) {
	int k = 0;
	do {
		if (stat[node[k].imax] < node[k].cutval) k = node[k].filsG; else k = node[k].filsD;
	} while (not node[k].terminal);
	return node[k].modmoy;
}

void TreeC::buildtree1(int seed, int i, int rep) {
	//cout<<"debut BUILDTREE\n";
	int k, kk;
	this->mw.randinit(seed + i, 3 * (i + seed));
	this->initree();
	//cout<<"apres initree de tree "<<i<<"\n";fflush(stdout);
	this->node[0].nsets = this->nsets;
	this->node[0].nvar = this->nvar;
	this->node[0].npassages = 0;
	this->node[0].pere = 0;
	this->node[0].indvar = vector<int>(this->nvar);
	this->mw.resample(rf.nstat, this->nvar, this->node[0].indvar);
	this->node[0].numset = vector<int>(this->nsets);
	for (int j = 0; j < this->nsets; j++) this->node[0].numset[j] = this->numset[j];
	//for (int m=0;m<10;m++) cout<<this->node[0].numset[m]<<"  ";cout<<"\n";
	k = 0;
	this->fin = false;//cout<<"avant la boucle while\n";
	while (not this->fin) {
		//cout<<k<<"\r";fflush(stdout);
		//cout<<" \nAVANT LE NOEUD "<<k<<"\n";
		//cout<<"Noeud k="<<k<<"   ";if (this->node[k].terminal) cout<<"terminal\n";else cout<<"not terminal\n";
		this->node[k].model = this->node[k].getdisval(this->mw);
		this->node[k].terminal = (this->node[k].model != -1);
		if (not this->node[k].terminal) { //cout<<"noeud non terminal\n";
			this->varused[this->node[k].imax] = true;//cout<<"Noeud k="<<k<<"   imax="<<this->node[k].imax<<"\n";
			importance[rep][this->node[k].imax] += this->node[k].delta;
			nimportance[rep][this->node[k].imax]++;
			this->node[k].filsG = k + 1;
			if (k == 2 * rf.nsel) {
				cout << "dépassement du nombre de noeuds 1 dans l'arbre " << i << "   kk=" << kk << "\n";
				cout << "k=" << k << "\n";
				for (kk = 0; kk <= 30; kk++) cout << "node=" << kk << "  pere=" << this->node[kk].pere << "  filsG=" << this->node[kk].filsG << "  filsD=" << this->node[kk].filsD << "  n=" << this->node[kk].npassages << "\n";
				exit(1);
			}
			k++;//cout<<"création du noeud "<<k<<" dans la descente sur "<<2*rf.nsel-1<<"\n";
			this->node[k].npassages = 0;
			this->node[k].pere = k - 1;
			this->node[k].nsets = this->node[this->node[k].pere].nsetG;
			this->node[k].nvar = this->nvar;
			this->node[k].numset = vector<int>(this->node[k].nsets);
			for (int j = 0; j < this->node[k].nsets; j++) this->node[k].numset[j] = this->node[this->node[k].pere].numsetG[j];
			this->node[k].indvar = vector<int>(this->nvar);
			this->mw.resample(rf.nstat, this->nvar, this->node[k].indvar);

			this->node[k].npassages++;
			//cout<<"noeud k="<<k<<"   npassages="<<this->node[k].npassages<<"\n";
		}
		else {//cout<<"noeud terminal\n";
			  //this->node[k].model = rf.model[this->node[k].numset[0]];
			kk = k;
			do {
				kk = this->node[kk].pere;
				//cout<<"on remonte d'un cran  kk="<<kk<<"   npassages="<<this->node[kk].npassages<<"\n";
			} while ((this->node[kk].npassages == 2) and (kk != 0));
			this->fin = ((kk == 0) and (this->node[kk].npassages == 1));
			if (not this->fin) {
				if (k == 2 * rf.nsel) {
					cout << "dépassement du nombre de noeuds 2 dans l'arbre " << i << "   kk=" << kk << "\n";
					cout << "k=" << k << "\n";
					for (kk = 0; kk <= k; kk++) {
						cout << "node=" << kk << "  pere=" << this->node[kk].pere;
						if (not this->node[kk].terminal) {
							cout << "  filsG=" << this->node[kk].filsG << " (" << this->node[kk].nsetG << ")";
							cout << "  filsD=" << this->node[kk].filsD << " (" << this->node[kk].nsetD << ")  n=" << this->node[kk].npassages << "\n";
						}
						else cout << "  terminal\n";
					}
					exit(1);
				}
				k++;//cout<<"création du noeud "<<k<<" apres remontée sur "<<2*nsel-1<<"\n";
				this->node[k].npassages = 0;
				this->node[k].pere = kk;
				this->node[kk].filsD = k;
				this->node[k].nsets = this->node[this->node[k].pere].nsetD;
				this->node[k].nvar = this->nvar;
				this->node[k].numset = vector<int>(this->node[k].nsets);
				for (int j = 0; j < this->node[k].nsets; j++) this->node[k].numset[j] = this->node[this->node[k].pere].numsetD[j];
				this->node[k].indvar = vector<int>(this->nvar);
				this->mw.resample(rf.nstat, this->nvar, this->node[k].indvar);
				this->node[kk].npassages++;
				this->node[k].npassages++;
			}
		}
	}
	this->nnodes = k + 1;
	//this->ecris(rt,i);
	for (int m = 0; m <= k; m++) {
		if (not this->node[m].terminal) {
			this->node[m].numsetD.clear();
			this->node[m].numsetG.clear();
			this->node[m].numset.clear();
			this->node[m].indvar.clear();
		}
	}
	ndone++;
	cout << "   construction de l'arbre " << ndone << "\r";
	fflush(stdout);
	//cout<<"fin de tree "<<i+1<<"   "<<tree[i].nnodes<<" noeuds\n";fflush(stdout);
}

void TreeC::buildtree2(int seed, int i, int rep) {
	//cout<<"debut BUILDTREE\n";
	int k, kk;
	this->mw.randinit(seed + i, 5 * (i + seed));
	this->initree();
	//cout<<"apres initree de tree "<<i<<"\n";fflush(stdout);
	this->node[0].nsets = this->nsets;
	this->node[0].nvar = this->nvar;
	this->node[0].npassages = 0;
	this->node[0].pere = 0;
	this->node[0].indvar = vector<int>(this->nvar);
	this->mw.resample(rf.nstat, this->nvar, this->node[0].indvar);
	this->node[0].numset = vector<int>(this->nsets);
	for (int j = 0; j < this->nsets; j++) this->node[0].numset[j] = this->numset[j];
	//for (int m=0;m<10;m++) cout<<this->node[0].numset[m]<<"  ";cout<<"\n";
	k = 0;
	this->fin = false;//cout<<"avant la boucle while\n";
	while (not this->fin) {
		//cout<<k<<"\r";fflush(stdout);
		//cout<<" \nAVANT LE NOEUD "<<k<<"\n";
		//cout<<"Noeud k="<<k<<"   ";if (this->node[k].terminal) cout<<"terminal\n";else cout<<"not terminal\n";
		this->node[k].modmoy = this->node[k].getdisval2(this->mw);
		//cout<<"model "<<this->node[k].model<<"\n";
		this->node[k].terminal = (this->node[k].modmoy > -0.5);
		if (not this->node[k].terminal) { //cout<<"noeud non terminal\n";
										  //this->varused[this->node[k].imax] = true;cout<<"Noeud k="<<k<<"   imax="<<this->node[k].imax<<"\n";
			this->node[k].filsG = k + 1;
			if (k == 2 * rf.nsel) {
				cout << "dépassement du nombre de noeuds 1 dans l'arbre " << i << "   kk=" << kk << "\n";
				cout << "k=" << k << "\n";
				for (kk = 0; kk <= 30; kk++) cout << "node=" << kk << "  pere=" << this->node[kk].pere << "  filsG=" << this->node[kk].filsG << "  filsD=" << this->node[kk].filsD << "  n=" << this->node[kk].npassages << "\n";
				exit(1);
			}
			k++;//cout<<"création du noeud "<<k<<" dans la descente sur "<<2*this->nsets-1<<"\n";
			this->node[k].npassages = 0;
			this->node[k].pere = k - 1;
			this->node[k].nsets = this->node[this->node[k].pere].nsetG;
			this->node[k].nvar = this->nvar;
			this->node[k].numset = vector<int>(this->node[k].nsets);
			for (int j = 0; j < this->node[k].nsets; j++) this->node[k].numset[j] = this->node[this->node[k].pere].numsetG[j];
			this->node[k].indvar = vector<int>(this->nvar);
			this->mw.resample(rf.nstat, this->nvar, this->node[k].indvar);

			this->node[k].npassages++;
			//cout<<"noeud k="<<k<<"   npassages="<<this->node[k].npassages<<"\n";
		}
		else {//cout<<"noeud terminal\n";
			  //this->node[k].model = rf.model[this->node[k].numset[0]];
			kk = k;
			do {
				kk = this->node[kk].pere;
				//cout<<"on remonte d'un cran  kk="<<kk<<"   npassages="<<this->node[kk].npassages<<"\n";
			} while ((this->node[kk].npassages == 2) and (kk != 0));
			this->fin = ((kk == 0) and (this->node[kk].npassages == 1));
			if (not this->fin) {
				if (k == 2 * rf.nsel) {
					cout << "dépassement du nombre de noeuds 2 dans l'arbre " << i << "   kk=" << kk << "\n";
					cout << "k=" << k << "\n";
					for (kk = 0; kk <= k; kk++) {
						cout << "node=" << kk << "  pere=" << this->node[kk].pere;
						if (not this->node[kk].terminal) {
							cout << "  filsG=" << this->node[kk].filsG << " (" << this->node[kk].nsetG << ")";
							cout << "  filsD=" << this->node[kk].filsD << " (" << this->node[kk].nsetD << ")  n=" << this->node[kk].npassages << "\n";
						}
						else cout << "  terminal\n";
					}
					exit(1);
				}
				k++;//cout<<"création du noeud "<<k<<" apres remontée sur "<<2*this->nsets-1<<"\n";
				this->node[k].npassages = 0;
				this->node[k].pere = kk;
				this->node[kk].filsD = k;
				this->node[k].nsets = this->node[this->node[k].pere].nsetD;
				this->node[k].nvar = this->nvar;
				this->node[k].numset = vector<int>(this->node[k].nsets);
				for (int j = 0; j < this->node[k].nsets; j++) this->node[k].numset[j] = this->node[this->node[k].pere].numsetD[j];
				this->node[k].indvar = vector<int>(this->nvar);
				this->mw.resample(rf.nstat, this->nvar, this->node[k].indvar);
				this->node[kk].npassages++;
				this->node[k].npassages++;
			}

		}
	}
	this->nnodes = k + 1;
	//this->ecris(rt,i);
	for (int m = 0; m <= k; m++) {
		if (not this->node[m].terminal) {
			this->node[m].numsetD.clear();
			this->node[m].numsetG.clear();
			this->node[m].numset.clear();
			this->node[m].indvar.clear();
		}

	}
	ndone++;
	cout << "   construction de l'arbre " << ndone << "\r";
	fflush(stdout);
}

void var_importance3(int rep) {
	//cout<<"\nRecherche des stats les plus informatives\n";
	int ns;
	ns = rf.nstatclass;
	if (ns > rf.nstat) ns = rf.nstat;
	vector<VMD> vd;
	vd = vector<VMD>(rf.nstat);
	for (int i = 0; i < rf.nstat; i++) {
		if (nimportance[rep][i] > 0) vd[i].x = importance[rep][i] / rf.ntrees /*/(double)nimportance[i]*/; else vd[i].x = 0.0;
		vd[i].name = nomstat[i];//cout<<"vd["<<i<<"].name = "<<vd[i].name<<"\n";
		for (int j = vd[i].name.length(); j < 12; j++) vd[i].name += " ";
	}
	//cout<<"avant le sort\n";
	sort(&vd[0], &vd[0] + rf.nstat);
	if (ns == rf.nstatclass) cout << " Classement des " << rf.nstatclass << " stats les plus informatives:\n";
	else cout << "\n\n Classement des stats selon leur informativité:\n";
	for (int i = 0; i < ns; i++) {
		cout << fixed << setw(4) << i + 1 << "  " << vd[i].name << "   " << fixed << setw(10) << setprecision(2) << vd[i].x;
		cout << "   (" << fixed << setw(6) << setprecision(2) << vd[i].x / vd[0].x * 100.0 << ")\n";
	}
	fout << "\nMost informative summary statistics\n";
	for (int i = 0; i < ns; i++) {
		fout << fixed << setw(4) << i + 1 << "  " << vd[i].name << "   " << fixed << setw(10) << setprecision(2) << vd[i].x;
		fout << "   (" << fixed << setw(6) << setprecision(2) << vd[i].x / vd[0].x * 100.0 << ")\n";
	}
}

bool RFC::sort_stat(int var, int i, int j) {
	return this->stat[i][var] < this->stat[j][var];
}

void RFC::readstat(bool LD) {
	cout << "\nLecture des données\n";
	int nscenOK = 0, nparamax = 0, iscen, bidon;
	enregC enr;
	bool scenOK;
	this->nmodel = rt.nscenchoisi;
	cout << "nmodel=" << this->nmodel << "\n";
	this->model = vector<int>(this->nsets);
	//for (int i=0;i<this->nmodel;i++) this->model[i]=rt.scenchoisi[i];
	this->nstat = rt.nstat;
	cout << "nstat = " << this->nstat << "   nsets=" << this->nsets << "\n";
	for (int i = 0; i < rt.nscen; i++) if (rt.nparam[i] > nparamax) nparamax = rt.nparam[i];
	enr.param = vector<float>(nparamax);
	enr.stat = vector<float>(this->nstat);
	this->stat = vector<vector<double>>(this->nsets);
	//LECTURE DE LA TABLE DE REFERENCE
	rt.filename = path + "reftableRF.bin";
	rt.openfile2();
	bidon = 0;
	cout << "ouverture de " << rt.filename << "\n";
	while ((nscenOK < this->nsets) and (bidon == 0)) {
		bidon = rt.readrecord(&enr);//cout<<"nscenOK="<<nscenOK<<"   bidon="<<bidon<<"\n";
		if (bidon == 0) {
			scenOK = false;
			iscen = 0;//cout<<"numscen="<<enr.numscen<<"    iscen="<<iscen<<"\n";
			while ((not scenOK) and (iscen < this->nmodel)) {
				scenOK = (enr.numscen == rt.scenchoisi[iscen]);
				iscen++;
			}
			if (scenOK) {
				this->stat[nscenOK] = vector<double>(this->nstat);
				for (int j = 0; j < this->nstat; j++) this->stat[nscenOK][j] = (double)enr.stat[j];
				iscen--;
				this->model[nscenOK] = iscen;// enr.numscen-1;
											 //if (nscenOK<10) cout<<"scenOK=true   enr.numscen="<<this->model[nscenOK]<<"   nscenOK="<<nscenOK<<"\n";
				nscenOK++;
			}
		}
	}
	rt.closefile();
	this->nsets = nscenOK;
	cout << "nsets = " << this->nsets << "\n";
	//LECTURE DES STATS OBSERVEES
	this->statobs = vector<double>(this->nstat);
	for (int j = 0; j < this->nstat; j++) this->statobs[j] = (double)header.stat_obs[j];

	this->statname = vector<string>(this->nstat);
	for (int i = 0; i < this->nstat; i++) this->statname[i] = header.statname[i];

	//AJOUT CONDITIONNEL DES LD
	if (LD) {
		cout << "\nAnalyse Discriminante pour l'ajout des LD...\n";
		long double *w, **X;
		int* mod;
		mod = new int[this->nsets];
		for (int i = 0; i < this->nsets; i++) mod[i] = this->model[i];
		vector<long double> statpiv;
		w = new long double[this->nsets];
		for (int i = 0; i < this->nsets; i++) w[i] = 1.0;
		X = new long double*[this->nsets];
		for (int i = 0; i < this->nsets; i++) {
			X[i] = new long double[this->nstat];
			for (int j = 0; j < this->nstat; j++) X[i][j] = (long double) this->stat[i][j];
		}
		afd = AFD(this->nsets, this->nstat, mod, w, X, 1.0);
		//cout<<"afd.nlambda = "<<afd.nlambda<<"\n";
		cout << "Fin de l'analyse discriminante\n";
		for (int i = 0; i < this->nsets; i++) {
			this->stat[i].resize(this->nstat + afd.nlambda - 1);
		}
		//calcul des LD sur les jeux simulés
		statpiv = vector<long double>(afd.nlambda);
		for (int i = 0; i < this->nsets; i++) {
			for (int j = 0; j < afd.nlambda - 1; j++) {
				statpiv[j] = 0.0;
				for (int k = 0; k < this->nstat; k++) statpiv[j] += (X[i][k] - afd.moy[k]) * afd.vectprop[k][j];
				this->stat[i][nstat + j] = (double)statpiv[j];
			}
		}
		//delete des vecteurs et tableaux
		for (int i = 0; i < this->nsets; i++) delete[]X[i];
		delete[]X;
		delete[]w;
		delete[]mod;
		//calcul des LD sur le jeu observé
		this->statobs.resize(this->nstat + afd.nlambda - 1);
		for (int j = 0; j < afd.nlambda - 1; j++) {
			statpiv[j] = 0.0;
			for (int k = 0; k < this->nstat; k++) statpiv[j] += ((long double)this->statobs[k] - afd.moy[k]) * afd.vectprop[k][j];
			this->statobs[nstat + j] = (double)statpiv[j];
		}
		statpiv.clear();
		//ajout des noms de LD
		this->statname = vector<string>(this->nstat);
		for (int i = 0; i < this->nstat; i++) this->statname[i] = header.statname[i];
		this->statname.resize(this->nstat + afd.nlambda - 1);
		for (int i = 0; i < afd.nlambda - 1; i++) this->statname[this->nstat + i] = "LD" + IntToString(i + 1);
		this->nstat += afd.nlambda - 1;
		cout << "Fin de l'ajout des composantes LD\n";
	}
	else cout << "Pas d'ajout des LD\n";
	nomstat = vector<string>(this->nstat);
	for (int i = 0; i < this->nstat; i++) nomstat[i] = this->statname[i];
	cout << "FIN de readstat nstat=" << this->nstat << "\n";
	cout << "Début du tri des stats" << endl;
	vector<int> stat_range = vector<int>(this->nsets);
	for (int i = 0; i < this->nsets; i++) stat_range[i] = i;
	this->stat_sorted_ind = vector<vector<int>>(this->nstat, vector<int>(stat_range));
	for (int i = 0; i < this->nstat; i++) {
		sort(this->stat_sorted_ind[i].begin(), this->stat_sorted_ind[i].end(), bind(&RFC::sort_stat, this, i, _1, _2));
	}
}

void dorandfor(string opt, int seed) {
	cout << "\nDébut de dorandfor\n";
	double duree;
	clock_t debut, debut_buildtree2;
	debut = clock();
	time_t rawtime;
	time(&rawtime);
	struct tm* timeinfo;
	timeinfo = localtime(&rawtime);
	int nrecpos, repmin, mod;
	double ppe, ppemin;
	string nomfiresult, s0, s1;
	nomfiresult = path + ident + "_randomforest.txt";
	progressfilename = path + ident + "_progress.txt";
	cout << "options : " << opt << "\n";
	vector<string> ss;
	vector<string> ss1;
	splitwords(opt, ";", ss);
	rf.ntrees = 500;
	int sssize = (int)ss.size();
	for (int i = 0; i < sssize; i++) {
		s0 = ss[i].substr(0, 2);
		s1 = ss[i].substr(2);
		if (s0 == "s:") {
			splitwords(s1, ",", ss1);
			rt.nscenchoisi = ss1.size();
			rt.scenchoisi = new int[rt.nscenchoisi];
			for (int j = 0; j < rt.nscenchoisi; j++) rt.scenchoisi[j] = atoi(ss1[j].c_str());
			nrecpos = 0;
			for (int j = 0; j < rt.nscenchoisi; j++) nrecpos += rt.nrecscen[rt.scenchoisi[j] - 1];
			cout << "scenarios à tester : ";
			for (int j = 0; j < rt.nscenchoisi; j++) {
				cout << rt.scenchoisi[j];
				if (j < rt.nscenchoisi - 1) cout << ",";
			}
			cout << "\n";
			ss1.clear();
		}
		else if (s0 == "n:") {
			rf.nsets = atoi(s1.c_str());
			cout << "nombre total de jeux de données considérés (tous scénarios confondus)= " << rf.nsets << "\n";
		}
		else if (s0 == "d:") {
			LD = true;
			cout << "ajout des linear discriminant scores\n";
		}
		else if (s0 == "t:") {
			rf.ntrees = atoi(s1.c_str());
			cout << "nombre d'arbres ' = " << rf.ntrees << "\n";
		}
		else if (s0 == "b:") {
			splitwords(s1, ",", ss1);
			rf.nbootsamp = ss1.size();
			rf.bootsamp = vector<int>(rf.nbootsamp);
			for (int j = 0; j < rf.nbootsamp; j++) rf.bootsamp[j] = atoi(ss1[j].c_str());
			cout << "taille des échantillons pour bootstrap : ";
			for (int j = 0; j < rf.nbootsamp; j++) {
				cout << rf.bootsamp[j];
				if (j < rf.nbootsamp - 1) cout << ",";
			}
			cout << "\n";
			ss1.clear();
		}
		else if (s0 == "k:") {
			rf.nvar = atoi(s1.c_str());
			flagk = true;
		}
		else if (s0 == "o:") {
			rf.nstatclass = atoi(s1.c_str());
			cout << "affichage des " << rf.nstatclass << " statistiques les plus efficaces\n";
			flago = true;
		}
	}
	cout << "fin de l'analyse des options\n\n";
	rf.readstat(LD);
	if (not flagk) rf.nvar = (int)sqrt(rf.nstat);
	if (not flago) rf.nstatclass = 30;
	fout.open(nomfiresult.c_str());
	fout << "DIYABC :                 Search of the best scenario through Random Forest                         " << asctime(timeinfo) << "\n";
	fout << "Data file       : " << header.datafilename << "\n";
	fout << "Reference table : " << rt.filename << "\n";
	fout << "Total number of simulated data sets : " << rf.nsets << "\n";
	if (LD) fout << "Linear discriminant scores have been added to summary statistics\n";
	fout << "The forest includes " << rf.ntrees << " trees\n";
	fout << "Calibration has been performed with trees containing ";
	for (int j = 0; j < rf.nbootsamp; j++) {
		fout << rf.bootsamp[j];
		if (j < rf.nbootsamp - 1) fout << ",";
	}
	fout << " data sets\n";
	fout << "Total number of summary statistics : " << rf.nstat << "    Number drawn at random at each node : " << rf.nvar << "\n\n";
	fout << "*************************Calibration**************************\n\n";
	ppemin = 1.0;
	importance = vector<vector<double>>(rf.nbootsamp);
	nimportance = vector<vector<int>>(rf.nbootsamp);
	obs_estscen = vector<int>(rt.nscenchoisi);
	sim_estscen = vector<vector<int>>(rf.nsets);
	rf.bienestime = vector<bool>(rf.nsets);
	vector<int> bienvu;
	bienvu = vector<int>(rf.nsets);
	obs_estim = vector<double>(rf.ntrees);
	for (int i = 0; i < rf.nsets; i++) sim_estscen[i] = vector<int>(rt.nscenchoisi);
	for (int rep = 0; rep < rf.nbootsamp; rep++) {
		rf.nsel = rf.bootsamp[rep];
		for (int j = 0; j < rt.nscenchoisi; j++) {
			obs_estscen[j] = 0;
			for (int i = 0; i < rf.nsets; i++) sim_estscen[i][j] = 0;
		}
		cout << "\nrep=" << rep + 1 << "    nsel=" << rf.nsel << "\n";
		cout << rf.nstat << " statistiques résumées en tout. Tirage de " << rf.nvar << " stat à chaque noeud\n";
		cout << rf.ntrees << " arbres à construire avec " << rf.nsel << " enregistrements par arbre\n";
		//fout<<rf.nstat<<" statistiques résumées en tout. Tirage de "<<rf.nvar<<" stat à chaque noeud\n";
		//fout<<rf.ntrees<<" arbres à construire avec "<<rf.nsel<<" enregistrements par arbre\n";
		fout << "-----------Trees with " << rf.bootsamp[rep] << " data sets-------------\n";
		importance[rep] = vector<double>(rf.nstat);
		nimportance[rep] = vector<int>(rf.nstat);
		for (int i = 0; i < rf.nstat; i++) {
			importance[rep][i] = 0.0;
			nimportance[rep][i] = 0;
		}
		ndone = 0;
#pragma omp parallel for shared(ndone,seed,rep,obs_estscen,sim_estscen,importance,nimportance) private(ctree) if(multithread)
		for (int i = 0; i < rf.ntrees; i++) {
			ctree.buildtree1(seed, i, rep);
			ctree.estim();
			ctree.deletree();
		}
		cout << "\nscenario     ";
		for (int i = 0; i < rt.nscenchoisi; i++) cout << "  " << rt.scenchoisi[i] << "  ";
		cout << "\nvotes      ";
		for (int i = 0; i < rt.nscenchoisi; i++) cout << setw(5) << obs_estscen[i];
		cout << "\n";
		ntestes = 0;
		njustes = 0;
		for (int i = 0; i < rf.nsets; i++) {
			mod = 0;
			for (int j = 1; j < rt.nscenchoisi; j++) if (sim_estscen[i][j] > sim_estscen[i][mod]) mod = j;
			rf.bienestime[i] = (rf.model[i] == mod);
			if (rf.bienestime[i]) njustes++;
		}
		ppe = (double)(rf.nsets - njustes) / (double)(rf.nsets);
		cout << "prior predictive error = " << ppe << "\n";
		if (ppe <= ppemin - 0.001) {
			ppemin = ppe;
			repmin = rep;
			for (int i = 0; i < rf.nsets; i++) if (rf.bienestime[i]) bienvu[i] = 1; else bienvu[i] = 0;
		}
	}
	cout << "\nRésultat de la calibration avec " << rf.nsets << " : nsel=" << rf.bootsamp[repmin] << "  (prior predictive error = " << ppemin << ")\n\n";
	fout << "\nResults of calibration with  " << rf.nsets << " overall data sets : \n";
	fout << "Minimum prior predictive error = " << ppemin << " obtained with trees of " << rf.bootsamp[repmin] << " data sets\n";
	var_importance3(repmin);
	for (int i = 0; i < rf.nsets; i++) rf.model[i] = bienvu[i];
	rf.nmodel = 2;
	rf.nvar = rf.nstat / 3;
	rf.nsel = rf.bootsamp[repmin];
	ndone = 0;
	//cout<<"avant omp parallel\n";
	debut_buildtree2 = clock();
#pragma omp parallel for shared(ndone,seed,obs_estscen,sim_estscen) private(ctree) if(multithread)
	for (int i = 0; i < rf.ntrees; i++) {
		ctree.buildtree2(seed, i, repmin);
		obs_estim[i] = ctree.inferobs(rf.statobs);
		ctree.deletree();
	}
	cout << endl << "Durée buildtree2 " << TimeToStr(walltime(debut_buildtree2)) << endl;
	double som = 0.0;
	for (int i = 0; i < rf.ntrees; i++) som += obs_estim[i];
	cout << "\nsom=" << som << "\n";
	som /= rf.ntrees;
	som = 1.0 - som;
	cout.precision(5);
	cout << "\n\nlocal error = " << som << "\n";
	fout.precision(5);
	fout << "\n\nlocal error = " << som << "\n";
	duree = walltime(debut);
	cout << "\n\ndurée totale =" << TimeToStr(duree) << "\n";
	fout << "\n\nTotal duration =" << TimeToStr(duree) << "\n";
	fout.close();
}

