/*
 * randforest.cpp
 * 
 * Created on : 24 november 2014
 * 
 */

#include <fstream>
#include <iostream>
#include <string>
#include <cstring>
#include "randforest.h"
#include "mesutils.h"
#include "matrices.h"
#include <stdlib.h>
#include "randomgenerator.h"
#include <algorithm>
#include <vector>
#include <math.h>
#include <iomanip>
#include "particleset.h"
#include <time.h>

#define LOG2 0.69314718055994530942;


using namespace std;

extern string ident, headerfilename;
extern bool multithread;
extern string progressfilename, path, ident;
extern string scurfile;
extern HeaderC header;
extern ReftableC rt;
extern int debuglevel;
extern enregC *enreg;
extern long double **phistarOK;
extern ParticleSetC ps;

RFC rf;
RFC *rfmin;

ofstream fout;
int nvoisins, nj;
double seuil=0.01;
string nomrfmin="foretmin.bin";
vector <string> nomstat;
resAFD afd;
bool LD=false,flagk=false,flago=false;


bool operator < (const VMC & lhs,const VMC & rhs) {
	return lhs.x < rhs.x;
}

bool operator < (const VMD & lhs,const VMD & rhs) {
	return lhs.x > rhs.x;
}

/**
 * Definition de l'operateur = pour une instance de la classe NodeC
 */
  NodeRC & NodeRC::operator= (NodeRC const & source) {
	if (this== &source) return *this;
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
		imax=(int)source.indvar.size();
		for (int i=0;i<imax;i++) this->indvar[i] = source.indvar[i];
	}
	if (not this->numset.empty()) this->numset.clear();
	if (not source.numset.empty()) {
		this->numset = vector<int>(source.numset.size());
		imax=(int)source.numset.size();
		for (int i=0;i<imax;i++) this->numset[i] = source.numset[i];
	}
	if (not this->numsetG.empty()) this->numsetG.clear();
	if (not source.numsetG.empty()) {
		this->numsetG = vector<int>(source.numsetG.size());
		imax=(int)source.numsetG.size();
		for (int i=0;i<imax;i++) this->numsetG[i] = source.numsetG[i];
	}
	if (not this->numsetD.empty()) this->numsetD.clear();
	if (not source.numsetD.empty()) {
		this->numsetD = vector<int>(source.numsetD.size());
		imax=(int)source.numsetD.size();
		for (int i=0;i<imax;i++) this->numsetD[i] = source.numsetD[i];
	}
	return *this;
}


/**
 * Definition de l'operateur = pour une instance de la classe TreeC
 */
  TreeC & TreeC::operator= (TreeC const & source) {
	if (this== &source) return *this;
	int imax;
	this->nnodes = source.nnodes;
	this->nsets = source.nsets;
	this->nvar = source.nvar;
	this->fin = source.fin;
	/// ON NE COPIE PAS LE GENERATEUR DE NOMBRES ALEATOIRES 
	if (not this->numset.empty()) this->numset.clear();
	if (not source.numset.empty()) {
		this->numset = vector<int>(source.numset.size());
		imax=(int)source.numset.size();
		for (int i=0;i<imax;i++) this->numset[i] = source.numset[i];
	}
	if (not this->indsel.empty()) this->indsel.clear();
	if (not source.indsel.empty()) {
		this->indsel = vector<int>(source.indsel.size());
		imax=(int)source.indsel.size();
		for (int i=0;i<imax;i++) this->indsel[i] = source.indsel[i];
	}
	if (not this->score.empty()) this->score.clear();
	if (not source.score.empty()) {
		this->score = vector<int>(source.score.size());
		imax=(int)source.score.size();
		for (int i=0;i<imax;i++) this->score[i] = source.score[i];
	}
	if (not this->index.empty()) this->index.clear();
	if (not source.index.empty()) {
		this->index = vector<int>(source.index.size());
		imax=(int)source.index.size();
		for (int i=0;i<imax;i++) this->index[i] = source.index[i];
	}
	if (not this->node.empty()) this->node.clear();
	if (not source.node.empty()) {
		this->node = vector<NodeRC>(source.node.size());
		imax=(int)source.node.size();
		for (int i=0;i<imax;i++) this->node[i] = source.node[i];
	}
	if (not this->varused.empty()) this->varused.clear();
	if (not source.varused.empty()) {
		this->varused = vector<bool>(source.varused.size());
		imax=(int)source.varused.size();
		for (int i=0;i<imax;i++) this->varused[i] = source.varused[i];
	}
	return *this;
}

/**
 * Definition de l'operateur = pour une instance de la classe RFC
 */
  RFC & RFC::operator= (RFC const & source) {
	if (this== &source) return *this;
	int imax,jmax;
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
		imax=(int)source.model.size();
		for (int i=0;i<imax;i++) this->model[i] = source.model[i];
	}
	if (not this->indsel.empty()) this->indsel.clear();
	if (not source.indsel.empty()) {
		this->indsel = vector<int>(source.indsel.size());
		imax=(int)source.indsel.size();
		for (int i=0;i<imax;i++) this->indsel[i] = source.indsel[i];
	}
	if (not this->nimportance.empty()) {
		imax=(int)this->nimportance.size();
		for (int i=0;i<imax;i++) {
			if (not this->nimportance[i].empty()) this->nimportance[i].clear();
		}
		this->nimportance.clear();
	}
	if (not source.nimportance.empty()) {
		this->nimportance = vector < vector  <int> >(source.nimportance.size());
		imax=(int)source.nimportance.size();
		for (int i=0;i<imax;i++) {
			if (not source.nimportance[i].empty()) {
				this->nimportance[i] = vector<int>(source.nimportance[i].size());
				jmax=(int)source.nimportance[i].size();
				for (int j=0;j<jmax;j++) this->nimportance[i][j] = source.nimportance[i][j];
			}
		}
	}
	if (not this->vote.empty()) {
		imax=(int)this->vote.size();
		for (int i=0;i<imax;i++) {
			if (not this->vote[i].empty()) this->vote[i].clear();
		}
		this->vote.clear();
	}
	if (not source.vote.empty()) {
		this->vote = vector < vector <double> >(source.vote.size());
		imax=(int)source.vote.size();
		for (int i=0;i<imax;i++) {
			if (not source.vote[i].empty()) {
				this->vote[i] = vector<double>(source.vote[i].size());
				jmax=(int)source.vote[i].size();
				for (int j=0;j<jmax;j++) this->vote[i][j] = source.vote[i][j];
			}
		}
	}
	if (not this->stat.empty()) {
		imax=(int)this->stat.size();
		for (int i=0;i<imax;i++) {
			if (not this->stat[i].empty()) this->stat[i].clear();
		}
		this->stat.clear();
	}
	if (not source.stat.empty()) {
		this->stat = vector < vector <double> >(source.stat.size());
		imax=(int)source.stat.size();
		for (int i=0;i<imax;i++) {
			if (not source.stat[i].empty()) {
				this->stat[i] = vector<double>(source.stat[i].size());
				jmax=(int)source.stat[i].size();
				for (int j=0;j<jmax;j++) this->stat[i][j] = source.stat[i][j];
			}
		}
	}
	if (not this->importance.empty()) {
		imax=(int)this->importance.size();
		for (int i=0;i<imax;i++) {
			if (not this->importance[i].empty()) this->importance[i].clear();
		}
		this->importance.clear();
	}
	if (not source.importance.empty()) {
		this->importance = vector < vector <double> >(source.importance.size());
		imax=(int)source.importance.size();
		for (int i=0;i<imax;i++) {
			if (not source.importance[i].empty()) {
				this->importance[i] = vector<double>(source.importance[i].size());
				jmax=(int)source.importance[i].size();
				for (int j=0;j<jmax;j++) this->importance[i][j] = source.importance[i][j];
			}
		}
	}
	if (not this->bootsamp.empty()) this->bootsamp.clear();
	if (not source.bootsamp.empty()) {
		this->bootsamp = vector<int>(source.bootsamp.size());
		imax=(int)source.bootsamp.size();
		for (int i=0;i<imax;i++) this->bootsamp[i] = source.bootsamp[i];
	}
	if (not this->varimp.empty()) this->varimp.clear();
	if (not source.varimp.empty()) {
		this->varimp = vector<double>(source.varimp.size());
		imax=(int)source.varimp.size();
		for (int i=0;i<imax;i++) this->varimp[i] = source.varimp[i];
	}
	if (not this->tree.empty()) this->tree.clear();
	if (not source.tree.empty()) {
		this->tree = vector<TreeC>(source.tree.size());
		imax=(int)source.tree.size();
		for (int i=0;i<imax;i++) this->tree[i] = source.tree[i];
	}
	if (not this->statobs.empty()) this->statobs.clear();
	if (not source.statobs.empty()) {
		this->statobs = vector<double>(source.statobs.size());
		imax=(int)source.statobs.size();
		for (int i=0;i<imax;i++) this->statobs[i] = source.statobs[i];
	}
	if (not this->statname.empty()) this->statname.clear();
	if (not source.statname.empty()) {
		this->statname = vector<string>(source.statname.size());
		imax=(int)source.statname.size();
		for (int i=0;i<imax;i++) this->statname[i] = source.statname[i];
	}
	return *this;
}
	
/*
	double NodeRC::caldisval(int nscen, int nsets, const vector<VMC>&vm, double val, vector <int>& nn) {
		double res,cg,cd,ggg,ddd,gd;
		vector <int>ga,dr;
		int gg,dd;
		ga = vector <int>(nscen);
		dr = vector <int>(nscen);
		for (int m=0;m<nscen;m++) {ga[m]=0;dr[m]=0;}
		gg=0;dd=0;
		for (int k=0;k<nsets;k++) {
			if (vm[k].x<val) {ga[vm[k].ind]++;gg++;}
			if (vm[k].x==val) {
				if (gg>dd) {dr[vm[k].ind]++;dd++;}
				else {ga[vm[k].ind]++;gg++;}
			}
			if (vm[k].x>val){dr[vm[k].ind]++;dd++;}
		}
		res=0.0;
		if (gg*dd==0) return 100.0;
		ggg=(double)gg;ddd=(double)dd;gd=ggg+ddd;
		for (int m=0;m<nscen;m++) {
			if (gg>0) cg= (double)ga[m]/ggg;
			if (dd>0) cd= (double)dr[m]/ddd;
			if (debuglevel==41) cout<<"caldisval nsets="<<nsets<<"   val="<<val<<"   gg="<<gg<<"    dd="<<dd<<"   ga["<<m<<"]="<<ga[m]<<"     dr["<<m<<"]="<<dr[m]<<"\n";
			res +=cg*(1.0-cg)*ggg+cd*(1.0-cd)*ddd; 
		}
		res=res/gd;
		nn[0]=gg;nn[1]=dd;
		if (res != res) {
			cout<<"val="<<val<<"\n";
			for (int k=0;k<nsets;k++) cout<<"("<<vm[k].ind<<")  "<<vm[k].x<<"   ";cout<<"\n";
			for (int m=0;m<nscen;m++) cout<<"ga["<<m<<"]="<<ga[m]<<"     dr["<<m<<"]="<<dr[m]<<"\n";
			cout<<"gg="<<gg<<"    dd="<<dd<<"\n";
			exit(1);
		}
		ga.clear();dr.clear();
		return res;
	}
*/	
	int NodeRC::regle3(vector <VMC> vm, vector <int>b, MwcGen& mw){
		bool ident;
		int modmax,nident,bmax;
		vector <double> xmin,xmax;
		xmin = vector <double>(this->nvar);
		xmax = vector <double>(this->nvar);
		for (int i=0;i<this->nvar;i++) {
			for (int j=0;j<nsets;j++) vm[j].ind=rf.model[numset[j]];
			vm[0].x=rf.stat[numset[0]][indvar[i]];xmax[i]=xmin[i]=vm[0].x;
			for (int j=1;j<nsets;j++){
				vm[j].x=rf.stat[numset[j]][indvar[i]];
				if (vm[j].x>xmax[i]) xmax[i]=vm[j].x;
				if (vm[j].x<xmin[i]) xmin[i]=vm[j].x;
			}
			ident=(xmax[i]-xmin[i])<1E-16*(xmax[i]+xmin[i]);
			if (not ident) {
				xmin.clear();xmax.clear();
				return -1;	//il existe au moins deux individus aux stat différentes
			}
		}
		modmax=0;bmax=b[0];
		for (int j=1;j<rf.nmodel;j++) if (b[j]>bmax) {bmax=b[j];modmax=j;}
		nident=0;for (int j=0;j<rf.nmodel;j++) if (b[j]==bmax) nident++;
		if (nident<2) { //
			xmin.clear();xmax.clear();
			return modmax; //tous les individus ont les mêmes valeurs de stat et le modèle modmax est majoritaire
		} else {
				int kk=mw.rand1(nident);
				int kkk=0;
				for (int j=0;j<rf.nmodel;j++) {
					if (b[j]==bmax) {
						kkk++;
						if(kkk==kk) modmax=j;
					}
				}
			xmin.clear();xmax.clear();
			return modmax;
		}
	} 
	
	double NodeRC::calGinimin(vector <VMC>& vm, double &cutvalc, vector <int> modfreq) {
		//cout<<"calGinimin 0\n";
		double epsilon,gini,res,cg,cd;
		int gg=0,dd=this->nsets,j=0,imax;
		vector <int> ga,dr;
		ga = vector <int>(rf.nmodel,0);
		dr = vector <int>(rf.nmodel);
		for (int i=0;i<rf.nmodel;i++) dr[i]=modfreq[i];
		vector <double> v;
		v.resize(0);
		//cout<<"calGinimin 1\n";
		epsilon=vm[this->nsets/2].x*1E-10;
		for (int i=1;i<this->nsets;i++) {
			if (vm[i].x-vm[i-1].x>epsilon) v.push_back(0.5*(vm[i].x+vm[i-1].x));
		}
		//cout<<"calGinimin v.size="<<v.size()<<"\n";
		gini=1.0;
		imax=(int)v.size();
		for (int i=0;i<imax;i++) {
			while (vm[j].x<v[i]) {
				ga[vm[j].ind]++;gg++;
				dr[vm[j].ind]--;dd--;
				j++;
			}
			//cout<<"calGinimin j="<<j<<"\n";
			res=0.0;
			for (int m=0;m<rf.nmodel;m++) {
				if (gg>0) cg= (double)ga[m]/(double)gg; else cg=0.0;
				if (dd>0) cd= (double)dr[m]/(double)dd; else cd=0.0;
				res +=cg*(1.0-cg)*(double)gg+cd*(1.0-cd)*(double)dd; 
			}
			res=res/(double)(gg+dd);
			//cout<<"calGinimin res="<<res<<"\n";
			if (res<gini) {
				gini=res;
				cutvalc=v[i];
			}
			//cout<<"calGinimin gini="<<gini<<"\n";
		}
		ga.clear();dr.clear();v.clear();
		return gini;
	}
	
	double NodeRC::calGini(vector <VMC>& vm, double cutval) {
		double cg,cd,gini=0.0;
		int gg=0,dd=0;
		vector <int> ga,dr;
		ga = vector <int>(rf.nmodel,0);
		dr = vector <int>(rf.nmodel,0);
		gini=1.0;
		for (int i=0;i<this->nsets;i++) {
			if (vm[i].x<cutval) {ga[vm[i].ind]++;gg++;}
			else	{dr[vm[i].ind]++;dd++;}
		}
		for (int m=0;m<rf.nmodel;m++) {
			if (gg>0) cg= (double)ga[m]/(double)gg; else cg=0.0;
			if (dd>0) cd= (double)dr[m]/(double)dd; else cd=0.0;
			gini +=cg*(1.0-cg)*(double)gg+cd*(1.0-cd)*(double)dd; 
		}
		gini /=(double)(gg+dd);
		ga.clear();dr.clear();
		return gini;
	}
	
	
	
	
	
/**
 *  teste si le noeud est terminal :
 *  retourne -1 si le noeud n'est pas terminal
 *  retourne le numéro du modèle si le noeud est terminal
 */
	int NodeRC::getdisval(MwcGen& mw) {
		double gini,ginimin,cutvalmin,dfa,c,deltamax,cutvalc=0.0;
		int modmax,freqmax,ii;
		vector <VMC> vm;
		vm = vector <VMC>(this->nsets);
		for (int j=0;j<this->nsets;j++) vm[j].ind=rf.model[this->numset[j]];
//calcul du nombre d'individus de chaque modèle
		vector <int> modfreq;
		modfreq = vector <int>(rf.nmodel,0);
		for (int j=0;j<this->nsets;j++) modfreq[vm[j].ind]++;
		modmax=0;for (int k=1;k<rf.nmodel;k++) if (modfreq[k]>modfreq[modmax]) modmax=k;
		if (modfreq[modmax]==this->nsets) { //tous les individus sont du même modèle
			vm.clear();
			modfreq.clear();
			return modmax;
		}
		modmax=regle3(vm,modfreq,mw);
		if (modmax !=-1) { //tous les individus ont les mêmes valeurs de stat
			vm.clear();
			modfreq.clear();
			return modmax;
		}
//calcul du Gini du noeud avant split
		dfa=0.0;for (int m=0;m<rf.nmodel;m++) {c=(double)modfreq[m]/(double)this->nsets;dfa +=c*(1.0-c);}
//calcul du nombre d'individus du modèle le plus fréquent à ce noeud
		freqmax=modfreq[0];for (int k=1;k<rf.nmodel;k++) if (modfreq[k]>freqmax) freqmax=modfreq[k];
		if (freqmax==1) {  //il y a au plus 1 seul individu par modèle
			do {
				ii=mw.rand0(this->nvar);
				for (int j=0;j<this->nsets;j++){
					vm[j].x=rf.stat[this->numset[j]][indvar[ii]];
					vm[j].ind=rf.model[this->numset[j]];
				}
				sort(&vm[0],&vm[nsets]);
				deltamax=fabs(vm[this->nsets-1].x-vm[0].x);
			} while (deltamax<fabs(vm[this->nsets/2].x*1E-10));
			cutvalmin=vm[this->nsets/2].x;
			//ginimin=calGini(vm,cutvalmin);
			this->imax=indvar[ii];
/*			for (int i=0;i<this->nvar;i++) {
				for (int j=0;j<this->nsets;j++){
					vm[j].x=rf.stat[this->numset[j]][indvar[i]];
					vm[j].ind=rf.model[this->numset[j]];
				}
				sort(&vm[0],&vm[nsets]);
				if (i==0) {
					deltamax=vm[nsets-1].x-vm[0].x;
					imax=indvar[i];
					cutvalmin=0.5*(vm[nsets-1].x+vm[0].x);
				} else {
					if (deltamax<vm[nsets-1].x-vm[0].x) {
						deltamax=vm[nsets-1].x-vm[0].x;
						imax=indvar[i];
						cutvalmin=0.5*(vm[nsets-1].x+vm[0].x);
					}
				}
			}*/
			for (int j=0;j<this->nsets;j++){
				vm[j].x=rf.stat[this->numset[j]][this->imax];
				vm[j].ind=rf.model[this->numset[j]];
			}
			ginimin=calGini(vm,cutvalmin);
		} else { //il y a plusieurs individus d'un modèle donné
			for (int i=0;i<this->nvar;i++){
				for (int j=0;j<this->nsets;j++){
					vm[j].x=rf.stat[this->numset[j]][indvar[i]];
					vm[j].ind=rf.model[this->numset[j]];
				}
				sort(&vm[0],&vm[nsets]);
				gini=calGinimin(vm,cutvalc,modfreq);
				if (i==0) {
					this->imax=indvar[i];
					cutvalmin=cutvalc;
					ginimin=gini;
				}else {
					if (gini<ginimin) {
						this->imax=indvar[i];
						cutvalmin=cutvalc;
						ginimin=gini;
					}
				}
			}
			for (int j=0;j<this->nsets;j++){
				vm[j].x=rf.stat[this->numset[j]][imax];
				vm[j].ind=rf.model[this->numset[j]];
			}
		}
		this->delta=(dfa-ginimin)*this->nsets;	
		this->disval=ginimin;	
		this->cutval=cutvalmin;
		this->numsetG.resize(0);this->numsetD.resize(0);
		for (int j=0;j<this->nsets;j++) {
			if (vm[j].x<this->cutval) this->numsetG.push_back(this->numset[j]);
			else this->numsetD.push_back(this->numset[j]);
		}
		this->nsetG=this->numsetG.size();
		this->nsetD=this->numsetD.size();
		vm.clear();modfreq.clear();
		return -1;
	}
	
	// calcul la règle au noeud courant
/*	bool NodeRC::getdisval(MwcGen& mw) {
		//double *dis,*cut,
		double c,d,da,db,dfa,epsilon,disG,disD,deltamax,xmin,xmax,xmed,fxmin,fxmax,fxmed,r=0.618,aa,bb;
		vector <VMC> vm;
		vector <int> nn,nnm;
		vector <double> f,b;
		vector <int> a;
		int k,kk,n,n1,m;
		bool ecr=true;
		disval=100.0;
		vm = vector <VMC>(rf.nsets);
		nn = vector <int>(2);
		nnm= vector <int>(2);
		f = vector <double>(rf.nmodel);
		a = vector <int>(rf.nmodel,0);
		b = vector <double>(rf.nmodel);
		
		for (int j=0;j<nsets;j++) vm[j].ind=rf.model[this->numset[j]];
		for (int j=0;j<nsets;j++) a[vm[j].ind]++;
		n=0;for (int m=0;m<rf.nmodel;m++) if (a[m]>0) n++;
		//dd=0.0;for (int m=0;m<rf.nmodel;m++) if ((a[m]>0)and(a[m]<nsets)) {c=(double)a[m]/(double)nsets;dd-=c*log(c)/LOG2;}
		dfa=0.0;for (int m=0;m<rf.nmodel;m++) {c=(double)a[m]/(double)nsets;dfa +=c*(1.0-c);}
		if (debuglevel==41) cout<<"dfa="<<dfa<<"\n";
		for (int m=0;m<rf.nmodel;m++) if (a[m]>0) b[m]=1.0/(double)a[m]; else b[m]=0.0;
		sort(&a[0],&a[rf.nmodel]);
		n1=(int)a[rf.nmodel-1];
		if (debuglevel==41) {for (int m=0;m<rf.nmodel;m++) cout<<"b["<<m<<"]="<<b[m]<<"  ";cout<<"    n1="<<n1<<"\n";}
		if (nsets==n1) {a.clear();b.clear();f.clear();vm.clear();return true;}//cout<<"nsets="<<nsets<<"     n1="<<n1<<"\n";ecr=false;


		if (n1>1){ 
			for (int i=0;i<nvar;i++) {
				for (int j=0;j<nsets;j++) vm[j].ind=rf.model[numset[j]];
				vm[0].x=rf.stat[numset[0]][indvar[i]];xmax=xmin=vm[0].x;
				for (int j=1;j<nsets;j++){
					vm[j].x=rf.stat[numset[j]][indvar[i]];
					if (vm[j].x>xmax) xmax=vm[j].x;
					if (vm[j].x<xmin) xmin=vm[j].x;
				}
				aa=xmin;bb=xmax;epsilon=0.5*(aa+bb)/10000.0;
				//cout<<"\naa="<<aa<<"   bb="<<bb<<"\n";
				xmin=r*aa+(1.0-r)*bb;
				xmax=aa+bb-xmin;
				fxmin=caldisval(rf.nmodel, nsets, vm, xmin,nn);
				fxmax=caldisval(rf.nmodel, nsets, vm, xmax,nn);
				if (debuglevel==41) cout<<"aa="<<aa<<"   bb="<<bb<<"   xmin="<<xmin<<"   xmax="<<xmax<<"   fxmin="<<fxmin<<"   fxmax="<<fxmax<<"\n";
				m=0;
				while ((fabs(bb-aa)>fabs(epsilon))and(m<100)) {
					//cout<<"avant  "<<fabs(bb-aa)<<"   "<<epsilon<<"\n";
					if (fxmin<=fxmax) {
						bb=xmax;
						xmax=xmin;
						xmin=r*aa+(1-r)*bb;
						fxmax=fxmin;
						fxmin=caldisval(rf.nmodel, nsets, vm, xmin,nn);
					} else {
						aa=xmin;
						xmin=xmax;
						xmax=(1-r)*aa+r*bb;
						fxmin=fxmax;
						fxmax=caldisval(rf.nmodel, nsets, vm, xmax,nn);
					}
					//cout<<"aa="<<aa<<"   bb="<<bb<<"   xmin="<<xmin<<"   xmax="<<xmax<<"   fxmin="<<fxmin<<"   fxmax="<<fxmax<<"\n";
					//cout<<"apres  "<<fabs(bb-aa)<<"   "<<epsilon<<"\n";
					m++;
				}
				if (m>99) {
					cout<<"m="<<m<<"\n";
					cout<<"aa="<<aa<<"   bb="<<bb<<"   xmin="<<xmin<<"   xmax="<<xmax<<"   fxmin="<<fxmin<<"   fxmax="<<fxmax<<"\n";
					cout<<"apres  "<<fabs(bb-aa)<<"   "<<epsilon<<"\n";
					
				}
				c=0.5*(aa+bb);
				d=caldisval(rf.nmodel, nsets, vm,c,nn);
				da=caldisval(rf.nmodel, nsets, vm,aa,nn);
				db=caldisval(rf.nmodel, nsets, vm,bb,nn);
				if ((da<d)and(da<db)) {d=da;c=aa;}
				else if ((db<d)and(db<da)){d=db;c=bb;}
				d=caldisval(rf.nmodel, nsets, vm,c,nn);
				if (debuglevel==41) cout<<"aa="<<aa<<"   bb="<<bb<<"   c="<<c<<"   d="<<d<<"    epsilon="<<epsilon<<"\n";
				//cout<<"\nvar="<<indvar[i]<<"   c="<<c<<"   d="<<d<<"    disval="<<disval<<"  nn[0]="<<nn[0]<<"  nn[1]="<<nn[1]<<"\n";
				if (disval>d) {
					disval=d;
					cutval=c;
					imax=indvar[i];
					delta = nsets*(dfa - d);
					nnm[0]=nn[0];
					nnm[1]=nn[1];
				}
				//cout<<"var="<<indvar[i]<<"   c="<<c<<"   d="<<d<<"    disval="<<disval<<"\n";
				if (debuglevel==41) {cout<<"var "<<indvar[i]<<"     ";for (int m=0;m<rf.nmodel;m++) cout<<"f["<<m<<"]="<<f[m]<<"   ";}
				if (debuglevel==41) {cout<<"    d="<<d<<"   dis="<<disval<<"   cut="<<cutval<<"\n";}
			}
		} else {
			//cout<<"     ATTENTION n1=1\n";
			disval=0.0;
			for (int i=0;i<nvar;i++) {
				for (int j=0;j<nsets;j++){vm[j].x=rf.stat[numset[j]][indvar[i]];}
				sort(&vm[0],&vm[nsets]);
				if (i==0) {
					deltamax=vm[nsets-1].x-vm[0].x;
					imax=indvar[i];
					cutval=0.5*(vm[nsets-1].x+vm[0].x);
				}
				else {
					if (deltamax<vm[nsets-1].x-vm[0].x) {
						deltamax=vm[nsets-1].x-vm[0].x;
						imax=indvar[i];
						cutval=0.5*(vm[nsets-1].x+vm[0].x);
					}
				}
			}
		}
		if (debuglevel==42) cout<<"imax="<<imax<<"   disval="<<disval<<"   cutval="<<cutval<<"    nvar="<<nvar<<"    nsets="<<nsets<<"\n-----------------------------------------------------\n";
		nsetG=0;nsetD=0;
		for (int j=0;j<nsets;j++) {
			if (rf.stat[numset[j]][imax]<cutval) nsetG++;
			if (rf.stat[numset[j]][imax]==cutval) {
				if (nsetD>=nsetG) nsetG++; else nsetD++;
			}
			if (rf.stat[numset[j]][imax]>cutval) nsetD++;
		}
		if ((nsetG<1)or(nsetD<1)){ 
			cout<<"nsets="<<nsets<<"    nsetG="<<nsetG<<"    nsetD="<<nsetD<<"   disval="<<disval<<"   cutval="<<cutval<<"   imax="<<imax<<"\n";
			//for (int j=0;j<nsets;j++) cout<<stat[numset[j]][imax]<<"   ";cout<<"\n";
		}
		numsetG = vector <int>(nsetG);
		numsetD = vector <int>(nsetD);
		nsetG=0;nsetD=0;
		//for (int j=0;j<nsets;j++) cout<<"("<<rt.model[numset[j]]<<"):"<<stat[numset[j]][imax]<<"  ";cout<<"\n";
		for (int j=0;j<nsets;j++) {
			//cout<<"indiv "<<numset[j]<<"   x="<<stat[numset[j]][imax]<<"  scen "<<rt.model[numset[j]];
			if (rf.stat[numset[j]][imax]<cutval) {
				numsetG[nsetG] = numset[j];//cout<<" -> à gauche\n";
				nsetG++;
			}  
			if (rf.stat[numset[j]][imax]==cutval) {
				if (nsetD>=nsetG){
					numsetG[nsetG] = numset[j];//cout<<" -> à gauche\n";
					nsetG++;
				} else {
					numsetD[nsetD] = numset[j];//cout<<" -> à droite\n";
					nsetD++;
				}
			} 
			if (rf.stat[numset[j]][imax]>cutval) {
					numsetD[nsetD] = numset[j];//cout<<" -> à droite\n";
					nsetD++;
			}
		}
		if (n1>1){
			if ((nnm[0]!=nsetG)or(nnm[1]!=nsetD)) {
				cout<<"imax="<<imax<<"\n";
				cout<<"nsetG="<<nsetG<<"   nnm[0]="<<nnm[0]<<"\n";
				cout<<"nsetD="<<nsetD<<"   nnm[1]="<<nnm[1]<<"\n";
				for (int j=0;j<nsets;j++) cout<<rf.stat[numset[j]][imax]<<" ("<< rf.model[numset[j]]<<")  ";cout<<" \n";
				exit(1);
			}
		}
		a.clear();b.clear();f.clear();vm.clear();
		return false;
	}
*/
	void TreeC::ecris(int num) {
		cout<<"\nEcriture de l'arbre "<<num<<"         "<<nnodes<<" noeuds\n";
		for (int i=0;i<nnodes;i++) {
			cout<<"Noeud "<<i<<"   pere "<<node[i].pere;
			if (not node[i].terminal )	{cout<<"   ("<<header.statname[node[i].imax]<<"   filsG "<<node[i].filsG<<"   filsD "<<node[i].filsD<<"   nsets="<<node[i].nsets;
				cout<<"   disval="<<node[i].disval<<"\n";
			}
			else cout<<"   TERMINAL   scenario "<<node[i].model<<"\n";
			for (int j=0;j<node[i].nsets;j++)cout<<"("<<rf.model[node[i].numset[j]]<<"):"<<node[i].numset[j]<<"  ";cout<<"\n";
		}
		cout<<"Fin de l'arbre "<<num<<"\n";
	}
	
	int TreeC::infermodel(const vector <double>& stat) {
		int k=0;
		do {
			if (stat[node[k].imax]<node[k].cutval) k=node[k].filsG; else k=node[k].filsD;
		} while (not node[k].terminal);
		return node[k].model;
	}
	
	int TreeC::infermodel2(double* stat) {
		int k=0;
		do {
			if (stat[node[k].imax]<node[k].cutval) k=node[k].filsG; else k=node[k].filsD;
		} while (not node[k].terminal);
		return node[k].model;
	}
	
	int TreeC::calprofondeur(const vector <double>& stat) {
		int k=0,n=0;
		do {
			if (stat[node[k].imax]<node[k].cutval) k=node[k].filsG; else k=node[k].filsD;
			n++;
		} while (not node[k].terminal);
		return n;
	}
	
	
	int TreeC::calprox(int n0, const vector <double>& stat,const vector <double>& statobs) {
	//cout<<"calprox  n0="<<n0<<"\n";
		int k=0,n=0;
		bool fin=false;
		while ((not fin)and(n<n0)) {
			if ((statobs[node[k].imax]<node[k].cutval) and (stat[node[k].imax]<node[k].cutval)) k=node[k].filsG;
			else if ((statobs[node[k].imax]>=node[k].cutval) and (stat[node[k].imax]>=node[k].cutval)) k=node[k].filsD;
			else fin=true;
			if (not fin) n++;//cout<<"n="<<n<<"\n";
		}
		return (n0-n);
	}
	
	void TreeC::initree(int seed, int i) {
		this->nsets = rf.nsel;
		this->nvar = rf.nvar;
		this->nnodes = 2*this->nsets-1;
		this->node = vector <NodeRC>(2*this->nsets);
		this->mw.randinit(seed+i,3*(i+seed));
		this->numset = vector <int>(this->nsets);
		this->indsel = vector <int>(this->nsets);
		this->score = vector <int>(rf.nstat);
		this->mw.samplewith(rf.nsets,this->nsets,this->indsel);
		for (int j=0;j<this->nsets;j++) this->numset[j] = this->indsel[j];
		this->varused = vector <bool>(rf.nstat,false);
	}
	
	void TreeC::deletree() {
		if (not this->numset.empty()) this->numset.clear();
		if (not this->indsel.empty()) this->indsel.clear();
		if (not this->score.empty()) this->score.clear();
		if (not this->index.empty()) this->index.clear();
		if (not this->varused.empty()) this->varused.clear();
		if (not this->node.empty()) {
			for (int m=0;m<this->nnodes;m++) {
				if (not this->node[m].indvar.empty()) this->node[m].indvar.clear();
				if (not this->node[m].numset.empty()) this->node[m].numset.clear();
				if (not this->node[m].numsetG.empty()) this->node[m].numsetG.clear();
				if (not this->node[m].numsetD.empty()) this->node[m].numsetD.clear();
			}
			this->node.clear();
		}
	}
	
	void RFC::dimimportance() {
		this->importance = vector< vector <double> > (this->nbootsamp);
		this->nimportance = vector< vector <int> > (this->nbootsamp);
	}
	/*void RFC::writetree(string filename, const ReftableC& rt) {
		
		
	}*/
	
	void RFC::growtrees(int seed, int rep) {
		//cout<<"debut de growtrees("<<rep<<")\n";
		int k,kk,ndone=0;
		importance[rep] = vector <double>(this->nstat);
		nimportance[rep] = vector <int>(this->nstat);
		for (int i=0;i<this->nstat;i++) {importance[rep][i]=0.0;nimportance[rep][i]=0;}
		cout<<ntrees<<" arbres à construire dans growtrees\n";
		tree = vector <TreeC>(this->ntrees);
#pragma omp parallel for shared(ndone) private(k,kk) if(multithread)
		for (int i=0;i<this->ntrees;i++) {
			tree[i].initree(seed,i);
			//cout<<"debut de tree "<<i<<"\n";fflush(stdout);
			tree[i].node[0].nsets=tree[i].nsets;
			tree[i].node[0].nvar=tree[i].nvar;
			tree[i].node[0].npassages = 0;
			tree[i].node[0].pere = 0;
			tree[i].node[0].indvar = vector <int> (tree[i].nvar);
			tree[i].mw.resample(this->nstat,tree[i].nvar,tree[i].node[0].indvar);
			tree[i].node[0].numset = vector <int>(tree[i].nsets);
			for (int j=0;j<tree[i].nsets;j++) tree[i].node[0].numset[j]=tree[i].numset[j];
			//for (int m=0;m<10;m++) cout<<tree[i].node[0].numset[m]<<"  ";cout<<"\n";
			k=0;
			tree[i].fin=false;//cout<<"avant la boucle while\n";
			while (not tree[i].fin) {
				//cout<<k<<"\r";fflush(stdout);
				//cout<<" \nAVANT LE NOEUD "<<k<<"\n";
				//cout<<"Noeud k="<<k<<"   ";if (tree[i].node[k].terminal) cout<<"terminal\n";else cout<<"not terminal\n";
				tree[i].node[k].model=tree[i].node[k].getdisval(tree[i].mw);
				tree[i].node[k].terminal = (tree[i].node[k].model!=-1);
				if (not tree[i].node[k].terminal){ //cout<<"noeud non terminal\n";
					tree[i].varused[tree[i].node[k].imax] = true;//cout<<"Noeud k="<<k<<"   imax="<<tree[i].node[k].imax<<"\n";
					importance[rep][tree[i].node[k].imax] +=tree[i].node[k].delta;
					nimportance[rep][tree[i].node[k].imax]++;
					tree[i].node[k].filsG=k+1;
					if (k==2*nsel) {
						cout<<"dépassement du nombre de noeuds 1 dans l'arbre "<<i<<"   kk="<<kk<<"\n";
						cout<<"k="<<k<<"\n";
						for (kk=0;kk<=30;kk++) cout<<"node="<<kk<<"  pere="<<tree[i].node[kk].pere<<"  filsG="<<tree[i].node[kk].filsG<<"  filsD="<<tree[i].node[kk].filsD<<"  n="<<tree[i].node[kk].npassages<<"\n";
						exit(1);
					}
					k++;//cout<<"création du noeud "<<k<<" dans la descente sur "<<2*nsel-1<<"\n";
					tree[i].node[k].npassages=0;
					tree[i].node[k].pere=k-1;
					tree[i].node[k].nsets=tree[i].node[tree[i].node[k].pere].nsetG;
					tree[i].node[k].nvar=tree[i].nvar;
					tree[i].node[k].numset = vector <int>(tree[i].node[k].nsets);
					for (int j=0;j<tree[i].node[k].nsets;j++) tree[i].node[k].numset[j]=tree[i].node[tree[i].node[k].pere].numsetG[j];
					tree[i].node[k].indvar = vector <int>(tree[i].nvar);
					tree[i].mw.resample(this->nstat,tree[i].nvar,tree[i].node[k].indvar);
					
					tree[i].node[k].npassages++;
					//cout<<"noeud k="<<k<<"   npassages="<<tree[i].node[k].npassages<<"\n";
					
				} else {//cout<<"noeud terminal\n";
					//tree[i].node[k].model = rf.model[tree[i].node[k].numset[0]];
					kk=k;
					do {
						kk=tree[i].node[kk].pere;
						//cout<<"on remonte d'un cran  kk="<<kk<<"   npassages="<<tree[i].node[kk].npassages<<"\n";
					} while ((tree[i].node[kk].npassages==2)and(kk!=0));
					tree[i].fin= ((kk==0)and(tree[i].node[kk].npassages==1));
					if (not tree[i].fin){
						if (k==2*nsel){
							cout<<"dépassement du nombre de noeuds 2 dans l'arbre "<<i<<"   kk="<<kk<<"\n";
							cout<<"k="<<k<<"\n";
							for (kk=0;kk<=k;kk++) { 
								cout<<"node="<<kk<<"  pere="<<tree[i].node[kk].pere;
								if (not tree[i].node[kk].terminal){ 
									cout<<"  filsG="<<tree[i].node[kk].filsG<<" ("<<tree[i].node[kk].nsetG<<")";
									cout<<"  filsD="<<tree[i].node[kk].filsD<<" ("<<tree[i].node[kk].nsetD<<")  n="<<tree[i].node[kk].npassages<<"\n";
								}else cout<<"  terminal\n";
							}
							exit(1);
						}
						k++;//cout<<"création du noeud "<<k<<" apres remontée sur "<<2*nsel-1<<"\n";
						tree[i].node[k].npassages=0;
						tree[i].node[k].pere=kk;
						tree[i].node[kk].filsD=k;
						tree[i].node[k].nsets=tree[i].node[tree[i].node[k].pere].nsetD;
						tree[i].node[k].nvar=tree[i].nvar;
						tree[i].node[k].numset = vector <int>(tree[i].node[k].nsets);
						for (int j=0;j<tree[i].node[k].nsets;j++) tree[i].node[k].numset[j]=tree[i].node[tree[i].node[k].pere].numsetD[j];
						tree[i].node[k].indvar = vector <int>(tree[i].nvar);
						tree[i].mw.resample(this->nstat,tree[i].nvar,tree[i].node[k].indvar);
						tree[i].node[kk].npassages++;
						tree[i].node[k].npassages++;
					}
				}
			}
			tree[i].nnodes=k+1;
			//tree[i].ecris(rt,i);
			for (int m=0;m<=k;m++) {
				if (not tree[i].node[m].terminal) {
					tree[i].node[m].numsetD.clear();
					tree[i].node[m].numsetG.clear();
					tree[i].node[m].numset.clear();
					tree[i].node[m].indvar.clear();
				}
			}
			ndone++;
			cout<<"construction de l'arbre "<<ndone<<"\r";fflush(stdout);
			//cout<<"fin de tree "<<i+1<<"   "<<tree[i].nnodes<<" noeuds\n";fflush(stdout);
		}
		cout<<"\n";
	}
	
	int RFC::bestmodel(int nscen,vector <double>& vote,const vector <double>& stat) {
		int infmodel,best,i;
		for (i=0;i<nscen;i++) vote[i]=0.0;
		for (i=0;i<ntrees;i++) {
			infmodel=tree[i].infermodel(stat);
			vote[infmodel]+=1.0;
		}
		best=0;
		for (i=1;i<nscen;i++) if (vote[i]>vote[best]) best=i;
		return best;
	}
	
	int RFC::bestmodel2(int k, int nscen,const vector <double>& stat) {
		int infmodel,best,i;
		vector <int>vote;
		bool participe;
		vote = vector <int>(nscen);
		for (i=0;i<nscen;i++) vote[i]=0;
		for (i=0;i<ntrees;i++) {
			participe=false;
			for (int j=0;j<tree[i].nsets;j++) {
				participe=(k==tree[i].numset[j]);
				if (participe) break;
			}
			if (not participe){
				infmodel=tree[i].infermodel(stat);
				vote[infmodel]++;
			}
		}
		best=0;
		for (i=1;i<nscen;i++) if (vote[i]>vote[best]) best=i;
		vote.clear();
		return best;
	}
	
	int RFC::bestmodel3(int k, int nscen,double* stat) {
		int infmodel,best,i;
		vector <int>vote;
		bool participe;
		vote = vector <int>(nscen);
		for (i=0;i<nscen;i++) vote[i]=0;
		for (i=0;i<ntrees;i++) {
			participe=false;
			for (int j=0;j<tree[i].nsets;j++) {
				participe=(k==tree[i].numset[j]);
				if (participe) break;
			}
			if (not participe){
				infmodel=tree[i].infermodel2(stat);
				vote[infmodel]++;
			}
		}
		best=0;
		for (i=1;i<nscen;i++) if (vote[i]>vote[best]) best=i;
		vote.clear();
		return best;
	}
	
	void RFC::infermodel() {
		int best;
		vector <double> vote;
		vote = vector <double>(rt.nscen);
		best=bestmodel(rt.nscen,vote,this->statobs);
		cout<<"Estimation du modèle pour les données observées\nScénario            ";
		for (int i=0;i<rt.nscen;i++) cout<<"   "<<i+1<<"    ";cout<<"\n"<<"Proportion de votes ";
		for (int i=0;i<rt.nscen;i++) cout<< std::fixed<<setw(4)<<setprecision(3)<<vote[i]<<"  "; cout<<"\n";
		cout<<"Le scénario "<<best+1<<" l'emporte\n";
		fout<<"\nInferred scenario for observed data\nScenario              ";
		for (int i=0;i<rt.nscen;i++) fout<<"  "<<i+1<<"   ";fout<<"\n"<<"Proportion of votes ";
		for (int i=0;i<rt.nscen;i++) fout<< std::fixed<<setw(4)<<setprecision(3)<<vote[i]<<"  "; fout<<"\n";
		fout<<"Scenario "<<best+1<<" is the winner\n";
		fout<<"\nScenario  number of trees\n";
		for (int i=0;i<rt.nscen;i++) fout<<"   "<<i+1<<"        "<<(int)vote[i]<<"\n";fout<<"\n";
		vote.clear();
	}
	
	double RFC::training_accuracy() {
		cout<<"début de training accuracy\n";
		double sfv,ppe;
		int sum,sumt=0,correct=0,ndone=0,nta=10000,dta;
		vector <vector <int> >vote;
		vector <vector <double> >fvote;
		vector <int>truemodel;
		vector <int>infmodel;
		vote  = vector <vector<int> >(this->nmodel);
		fvote = vector <vector<double> >(this->nmodel);
		truemodel = vector <int>(nta);
		infmodel = vector <int>(nta);
		dta=nsets-nta;
		for (int i=0;i<rt.nscen;i++) {
			vote[i] = vector <int>(this->nmodel);
			fvote[i] = vector <double>(this->nmodel);
			for (int j=0;j<rt.nscen;j++) vote[i][j]=0;
		}
		//cout<<"avant la boucle multithread\n";
#pragma omp parallel for shared(truemodel,infmodel,ndone) if(multithread)
		for (int i=0;i<nta;i++) {
			truemodel[i] = this->model[i+dta];
			infmodel[i] = bestmodel2(i+dta,this->nmodel,this->stat[i+dta]);
			ndone++;
			if (((ndone+1) % 100)==0) cout<<"test de l'individu "<<ndone+1<<" / "<<nta<<"                      \r";fflush(stdout);
			//cout<<"trumodel["<<i<<"] = "<<truemodel[i]<<"     infmodel["<<i<<"] = "<<infmodel[i]<<"\n";
		}
		cout<<"\n";
		//cout<<"après la boucle multithread\n";
		for (int i=0;i<nta;i++) vote[truemodel[i]][infmodel[i]]++;
		cout<<"\nTraining accuracy\n";
		//cout<<"En données brutes\n True\\Estim.\n";
		//cout<<"    ";for (int i=0;i<rt.nscen;i++) cout<<"    "<<i+1<<"   ";cout<<"\n";
		fout<<"\nTraining accuracy\n";
		fout<<"Raw data\n True\\Estim.\n";
		fout<<"    ";for (int i=0;i<this->nmodel;i++) fout<<"    "<<i+1<<"   ";fout<<"\n";
		/*for (int i=0;i<rt.nscen;i++) {
			cout<< std::fixed<<setw(3)<<i+1<<" ";
			for (int j=0;j<rt.nscen;j++) cout<<setw(6)<<setprecision(0)<<vote[i][j]<<"  ";cout<<"\n";
		}*/
		for (int i=0;i<this->nmodel;i++) {
			fout<<std::fixed<<setw(3)<<i+1<<" "<< std::fixed;
			for (int j=0;j<rt.nscen;j++) fout<<setw(6)<<setprecision(0)<<vote[i][j]<<"  ";fout<<"\n";
		}
		for (int i=0;i<this->nmodel;i++) {
			sum=0;for (int j=0;j<this->nmodel;j++) sum +=vote[i][j];
			sumt +=sum;
			correct +=vote[i][i];
			for (int j=0;j<this->nmodel;j++) fvote[i][j] =(double)vote[i][j]/(double)sum;
		}
		//cout<<"Somme totale = "<<sumt<<"   (attendus : "<<nta<<")\n";
		//cout<<"\nEn données normalisées (somme de chaque ligne=1)\n True\\Estim.\n";
		//cout<<"   ";for (int i=0;i<rt.nscen;i++) cout<<"    "<<i+1<<"  ";
		//cout<<" class.error\n";
		for (int i=0;i<this->nmodel;i++) {
			//cout<<std::fixed<<setw(3)<<i+1<<"  "<< std::fixed;
			//for (int j=0;j<rt.nscen;j++) cout<<setw(5)<<setprecision(3)<<fvote[i][j]<<"  ";
			sfv=0.0;for (int j=0;j<rt.nscen;j++) if (j!=i) sfv+=fvote[i][j];
			//cout<<setw(8)<<setprecision(4)<<sfv<<"\n";
		}
		fout<<"Sum = "<<sumt<<"   (expected : "<<nta<<")\n";
		fout<<"\nNormalized data (sum of each line = 1)\n True\\Estim.\n";
		fout<<"   ";for (int i=0;i<this->nmodel;i++) fout<<"    "<<i+1<<"  ";fout<<" class.error\n";
		for (int i=0;i<this->nmodel;i++) {
			fout<<std::fixed<<setw(3)<<i+1<<"  "<< std::fixed;
			for (int j=0;j<rt.nscen;j++) fout<<setw(5)<<setprecision(3)<<fvote[i][j]<<"  ";
			sfv=0.0;for (int j=0;j<this->nmodel;j++) if (j!=i) sfv+=fvote[i][j];
			fout<<setw(8)<<setprecision(4)<<sfv<<"\n";
		}
		ppe = (double)(sumt-correct)/(double)sumt;
		cout<<"Prior predictive error : "<<fixed<<setw(6)<<setprecision(3)<<ppe<<"\n";
		fout<<"\nPrior predictive error : "<<fixed<<setw(6)<<setprecision(3)<<ppe<<"\n";
		for (int i=0;i<rt.nscen;i++) vote[i].clear();
		for (int i=0;i<rt.nscen;i++) fvote[i].clear();
		vote.clear();
		fvote.clear();
		truemodel.clear();
		infmodel.clear();
		return ppe;
	}
	
	void var_importance3(int rep) {
		//cout<<"\nRecherche des stats les plus informatives\n";
		int ns;
		ns=rfmin->nstatclass;if (ns>rfmin->nstat) ns=rfmin->nstat;
		vector <VMD> vd;
		vd =  vector <VMD>(rfmin->nstat);
		for (int i=0;i<rfmin->nstat;i++) {
			if (rfmin->nimportance[rep][i]>0) vd[i].x=rfmin->importance[rep][i]/rfmin->ntrees /*/(double)nimportance[i]*/;else vd[i].x=0.0;
			vd[i].name = nomstat[i];//cout<<"vd["<<i<<"].name = "<<vd[i].name<<"\n";
			for (int j=vd[i].name.length();j<12;j++) vd[i].name +=" ";
		}
		//cout<<"avant le sort\n";
		sort(&vd[0],&vd[rfmin->nstat]);
		if (ns==rfmin->nstatclass) cout<<" Classement des "<<rfmin->nstatclass<<" stats les plus informatives:\n";
		else cout<<"\n\n Classement des stats selon leur informativité:\n";
		for (int i=0;i<ns;i++) {
			cout<<fixed<<setw(4)<<i+1<<"  "<<vd[i].name<<"   "<<fixed<<setw(10)<<setprecision(2)<<vd[i].x;
			cout<<"   ("<<fixed<<setw(6)<<setprecision(2)<<vd[i].x/vd[0].x*100.0<<")\n";
		}
		fout<<"\nMost informative summary statistics\n";
		for (int i=0;i<ns;i++) {
			fout<<fixed<<setw(4)<<i+1<<"  "<<vd[i].name<<"   "<<fixed<<setw(10)<<setprecision(2)<<vd[i].x;
			fout<<"   ("<<fixed<<setw(6)<<setprecision(2)<<vd[i].x/vd[0].x*100.0<<")\n";
		}
	}

	void RFC::readstat(bool LD) {
		cout<<"\nLecture des données\n";
		int nscenOK=0,nparamax=0,iscen,bidon;
        enregC enr;
		bool scenOK;
		this->nmodel=rt.nscenchoisi;cout<<"nmodel="<<this->nmodel<<"\n";
		this->model = vector <int>(this->nsets);
		//for (int i=0;i<this->nmodel;i++) this->model[i]=rt.scenchoisi[i];
		this->nstat=rt.nstat;cout<<"nstat = "<<this->nstat<<"   nsets="<<this->nsets<<"\n";
        for (int i=0;i<rt.nscen;i++) if (rt.nparam[i]>nparamax) nparamax=rt.nparam[i];
        enr.param = vector <float>(nparamax);
        enr.stat  = vector <float>(this->nstat);
		this->stat = vector <vector <double> >(this->nsets);
//LECTURE DE LA TABLE DE REFERENCE
		rt.filename=path+"reftableRF.bin";
        rt.openfile2();bidon=0;cout<<"ouverture de "<<rt.filename<<"\n";
		while ((nscenOK<this->nsets)and(bidon==0)) {
			bidon=rt.readrecord(&enr);//cout<<"nscenOK="<<nscenOK<<"   bidon="<<bidon<<"\n";
			if (bidon==0){
				scenOK=false;iscen=0;//cout<<"numscen="<<enr.numscen<<"    iscen="<<iscen<<"\n";
				while((not scenOK)and(iscen<this->nmodel)) {
					scenOK=(enr.numscen==rt.scenchoisi[iscen]);
					iscen++;
				}
				if (scenOK) {
					this->stat[nscenOK] = vector <double>(this->nstat);
					for (int j=0;j<this->nstat;j++) this->stat[nscenOK][j] = (double)enr.stat[j]; 
					this->model[nscenOK] = enr.numscen-1;
					//if (nscenOK<10) cout<<"scenOK=true   enr.numscen="<<this->model[nscenOK]<<"   nscenOK="<<nscenOK<<"\n";
					nscenOK++;	
				}
			}
		}
		rt.closefile();
		this->nsets=nscenOK;
		cout<<"nsets = "<<this->nsets<<"\n";
//LECTURE DES STATS OBSERVEES
		this->statobs = vector <double>(this->nstat);
		for (int j=0;j<this->nstat;j++) this->statobs[j] = (double)header.stat_obs[j];

		this->statname = vector <string>(this->nstat);
		for (int i=0;i<this->nstat;i++) this->statname[i]=header.statname[i];
		
//AJOUT CONDITIONNEL DES LD
		if (LD) {
			cout<<"\nAnalyse Discriminante pour l'ajout des LD...\n";
			long double *w,**X;
			int *mod;
			mod = new int[this->nsets];for (int i=0;i<this->nsets;i++) mod[i]=this->model[i];
			vector <long double>statpiv;
			w = new long double[this->nsets];
			for (int i=0;i<this->nsets;i++) w[i] = 1.0;
			X = new long double*[this->nsets];
			for (int i=0;i<this->nsets;i++) {
				X[i] = new long double[this->nstat];
				for (int j=0;j<this->nstat;j++) X[i][j] = (long double) this->stat[i][j];
			}
			afd = AFD(this->nsets,this->nstat,mod,w,X,1.0);
			//cout<<"afd.nlambda = "<<afd.nlambda<<"\n";
			cout<<"Fin de l'analyse discriminante\n";
			for (int i=0;i<this->nsets;i++){this->stat[i].resize(this->nstat+afd.nlambda-1);}
		//calcul des LD sur les jeux simulés
			statpiv = vector <long double>(afd.nlambda);
			for (int i=0;i<this->nsets;i++) {
				for (int j=0;j<afd.nlambda-1;j++){
					statpiv[j] =0.0;
					for (int k=0;k<this->nstat;k++) statpiv[j] += (X[i][k] - afd.moy[k])*afd.vectprop[k][j];
					this->stat[i][nstat+j] = (double)statpiv[j];
				}
			}
		//delete des vecteurs et tableaux
			for (int i=0;i<this->nsets;i++) delete []X[i];delete []X;
			delete []w;delete []mod;
		//calcul des LD sur le jeu observé
			this->statobs.resize(this->nstat+afd.nlambda-1);
			for (int j=0;j<afd.nlambda-1;j++){
				statpiv[j] =0.0;
				for (int k=0;k<this->nstat;k++) statpiv[j] += ((long double)this->statobs[k] - afd.moy[k])*afd.vectprop[k][j];
				this->statobs[nstat+j] = (double)statpiv[j];
			}
			statpiv.clear();
		//ajout des noms de LD
			this->statname = vector <string>(this->nstat);
			for (int i=0;i<this->nstat;i++) this->statname[i]=header.statname[i];
			this->statname.resize(this->nstat+afd.nlambda-1);
			for (int i=0;i<afd.nlambda-1;i++) this->statname[this->nstat+i]="LD"+IntToString(i+1);
			this->nstat+=afd.nlambda-1;
			cout<<"Fin de l'ajout des composantes LD\n";
		} else cout<<"Pas d'ajout des LD\n";
		nomstat = vector<string>(this->nstat);
		for (int i=0;i<this->nstat;i++) nomstat[i]=this->statname[i];
		cout<<"FIN de readstat nstat="<<this->nstat<<"\n";
	}

	void calposterior(int seed) {
		cout<<"\nCalcul de la Posterior Predictive Error\n";
		int k,ii,debut=0,nt,np;
		nvoisins = (int)((double)rfmin->nsel*seuil); 
		nj=10000/nvoisins;cout<<"nvoisins="<<nvoisins<<"   nj="<<nj<<"    total="<<nj*nvoisins<<"\n";	
		//Calcul de la profondeur de l'observation dans chaque arbre
		//vector <int> profobs;
		vector <long double>statpiv;
		/*profobs = vector<int>(rfmin->ntrees);
		for (int i=0;i<rfmin->ntrees;i++) {
			profobs[i]=rfmin->tree[i].calprofondeur(rfmin->statobs);
			//cout<<"Profondeur de l'observé dans l'arbre "<<i+1<<" = "<<profobs[i]<<"\n";
		}*/
/*cout<<"vérification des datasets de l'arbre 0 nsets = "<<rfmin->tree[0].nsets<<"\n";
for (k=0;k<rfmin->tree[0].nsets;k++) cout<<rfmin->tree[0].numset[k]<<"  ";
cout<<"\n\n";*/
		
		cout<<"Calcul de la proximité des jeux simulés  rfmin->nsel="<<rfmin->nsel<<"    rfmin->nsel="<<rfmin->nsel<<"\n";
	//Calcul de la proximité des jeux simulés
		VMC *prox;
		prox = new VMC[rfmin->nsel];
		bool participe,fin;
		for (int i=0;i<rfmin->nsel;i++) {
			prox[i].ind = i;
			prox[i].x=0;
			nt=0;np=0;
			for (int j=0;j<rfmin->ntrees;j++) {
				participe=false;
				k=i;
				ii=rfmin->tree[j].numset[k];
				participe=(i==ii);
				fin=false;
				while ((not fin)and((not participe))) {
					if(i>ii){ii=rfmin->tree[j].numset[++k];fin=(i<ii);}
					else {ii=rfmin->tree[j].numset[--k];fin=i>ii;}
					participe=(i==ii);
					/*cout<<"i="<<i<<"   k="<<k<<"   ii="<<ii;
					if (fin) cout<<"   fin";
					if (participe)  cout<<"   participe";
					cout<<"\n";*/
				}
				/*
				for (k=0;k<rfmin->tree[j].nsets;k++) {
					participe=(i==rfmin->tree[j].numset[k]);
					if (participe) break;
				}*/
				if (not participe) {
					nt++;
					k=0;fin=false;
					while ((not rfmin->tree[j].node[k].terminal)and(not fin)) {
						ii=rfmin->tree[j].node[k].imax;
						if ((rfmin->statobs[ii]<rfmin->tree[j].node[k].cutval) and (rfmin->stat[i][ii]<rfmin->tree[j].node[k].cutval)) k=rfmin->tree[j].node[k].filsG;
						else if ((rfmin->statobs[ii]>=rfmin->tree[j].node[k].cutval) and (rfmin->stat[i][ii]>=rfmin->tree[j].node[k].cutval)) k=rfmin->tree[j].node[k].filsD;
						else fin=true;
						if ((not fin)and(rfmin->tree[j].node[k].terminal)) np++; 
					}
				}
			}
			//cout<<"dataset "<<i<<"    np="<<np<<"   sur nt="<<nt<<setw(4)<<"   prox="<<(double)np/(double)nt<<"\r";fflush(stdout);
			//if ((double)np/(double)nt>0.001) cout<<"\n";
			prox[i].x = 1.0-(double)np/(double)nt;
			//prox[i].x += (double)rfmin->tree[j].calprox(profobs[j],rfmin->stat[i],rfmin->statobs);
			cout<<"prox["<<i<<"].ind = "<<prox[i].ind<<"     x = "<<prox[i].x<<"\r";
		}
		cout<<"\nClassement des plus proches voisins\n";
		sort(&prox[0],&prox[rfmin->nsel]);
		for (int i=0;i<rfmin->nsel;i++) prox[i].x=1.0-prox[i].x;
		cout<<"\n";
		for (int i=0;i<10;i++) cout<<"prox["<<i<<"].ind = "<<prox[i].ind<<"     x = "<<prox[i].x<<"\n";
		cout<<"..................\n";
		for (int i=nvoisins-11;i<nvoisins;i++) cout<<"prox["<<i<<"].ind = "<<prox[i].ind<<"     x = "<<prox[i].x<<"\n";
		cout<<"\n";
		double proxlim=prox[nvoisins-1].x;
		int nv=0;
		for (int p=0;p<nvoisins;p++) if (prox[p].x>proxlim) nv++;
		enreg = new enregC[rfmin->nsel];
		rt.openfile2();
		int nphistarOK;
		//double **stata;
		double **statb;
		vector <vector <long double> >phistarOKK;
		for (int p=0;p<rfmin->nsel;p++) {
			//cout<<"lecture de l'enregistrement "<<p;
			enreg[p].stat = vector <float>(rfmin->nstat);
			enreg[p].param = vector <float>(rt.nparamax);
			enreg[p].numscen = rt.readparam(enreg[p].param);
			//cout<<"   scénario "<<enreg[p].numscen<<"\n";
		}
		rt.closefile(); 
		cout<<"fin de la lecture des "<<rfmin->nsel<<" enregistrements de la reftable\n";
		statb = new double* [nj];for (int p=0;p<nj;p++) statb[p] = new double [rfmin->nstat];
		cout<<"nombre de colonnes de statb :"<<rfmin->nstat<<"\n";
		if (LD) {statpiv = vector <long double>(header.nstat);cout<<"nombre de colonnes de statpiv :"<<header.nstat<<"\n";}
		int nerror=0,kk,nerrorscen,scenestim;
		vector <int> nscenestim;
		nscenestim = vector <int>(rt.nscen);
		phistarOKK=vector <vector <long double> >(nj);for (int i=0;i<nj;i++) phistarOKK[i] = vector <long double> (rt.nparamax);
		cout<<"                                                ";for (int j=0;j<rt.nscen;j++) cout<<setw(5)<<j+1;cout<<"\n";
		fout<<"\nSelection of "<<nv<<" nearest neighbors\n";
		fout<<"prox[closest] = "<<prox[0].x<<"   prox[500-th closest]"<<prox[499].x<<"\n";
		fout<<"Simulating "<<nj*nv<<" datasets ("<<nj<<" datasets per neighbor)\n";
		fout<<" True\\Estim.\n        ";
		for (int j=0;j<rt.nscen;j++) fout<<setw(4)<<j+1<<"   ";fout<<"\n";
		for (int iscen=0;iscen<rt.nscen;iscen++) {
			nphistarOK=0;nerrorscen=0;kk=0;
			for (int j=0;j<rt.nscen;j++) nscenestim[j]=0;
			for (int p=0;p<nv;p++) if (enreg[prox[p].ind].numscen-1==iscen) nphistarOK++;
			if (nphistarOK>0){
				for (int p=0;p<nv;p++) if (enreg[prox[p].ind].numscen-1==iscen) {
					for (int i=0;i<nj;i++) {
						//cout<<"remplissage des paramètres de phistarOK["<<i<<"] pour le dataset "<<prox[p].ind<<"\n";
						for (int j=0;j<rt.nparam[iscen];j++) {
							phistarOKK[i][j]=(long double)enreg[prox[p].ind].param[j];
							//cout<<setw(15)<<setprecision(9)<<"     param["<<j<<"] = "<<phistarOKK[i][j]<<"\n";
						}
					}
					//cout<<"avant dosimulstat2\n";
					ps.dosimulstat2(debut,nj,phistarOKK,false,multithread,iscen+1,seed,statb);
					kk +=nj;
					cout<<"Sur "<<setw(4)<<kk<<" particules (scenario "<<iscen+1<<")   ";
					if (debut==0) debut = 10;
					//cout<<"apres dosimulstat2\n";
					if (LD) {
						for (int i=0;i<nj;i++) {
							for (int j=0;j<afd.nlambda-1;j++){
								statpiv[j]=0.0;
								for (int k=0;k<header.nstat;k++) statpiv[j] += (statb[i][k] - afd.moy[k])*afd.vectprop[k][j];
								statb[i][header.nstat+j] = (double)statpiv[j];
							}
						}
					}
					for (int i=0;i<nj;i++) {
						scenestim=rfmin->bestmodel3(prox[p].ind,rt.nscen,statb[i]);
						//for(int j=0;j<rfmin->nstat;j++) cout<<"  "<<statb[i][j];cout<<"\n";
						//cout<<"i="<<i<<"   scenestim="<<scenestim<<"\n";
						nscenestim[scenestim]++;
						if (scenestim != iscen) {nerror++;nerrorscen++;}
					}
					cout<<setw(4)<<nerrorscen<<" errors :";
					for (int j=0;j<rt.nscen;j++) cout<<setw(5)<<nscenestim[j];
					cout<<"    dataset "<<prox[p].ind;
					cout<<"\r";fflush(stdout);
				}
				cout<<"\n";
				fout<<"  "<<iscen+1<<"    ";
				for (int j=0;j<rt.nscen;j++) fout<<setw(5)<<nscenestim[j]<<"  ";fout<<"\n"; 
			}
			
		}
		
		double pperror;
		pperror = (double)nerror/(double)(nv*nj);
		cout<<"\n"<<nerror<<" errors over "<<nv*nj<<" datasets ("<<nv<<")\n\n";
		cout<<"\nPosterior predictive error = "<<pperror<<"\n";
		fout<<"\nPosterior predictive error = "<<pperror<<"\n";
	}
	
	void TreeC::ecrifich(ofstream& foret) {
		int n=this->node.size();		
		foret.write((char*)&n,sizeof(int));
		for (int i=0;i<n;i++) {
			foret.write((char*)&(this->node[i].cutval),sizeof(double));
			foret.write((char*)&(this->node[i].filsG),sizeof(int));
			foret.write((char*)&(this->node[i].filsD),sizeof(int));
			foret.write((char*)&(this->node[i].terminal),sizeof(bool));
			if (this->node[i].terminal) foret.write((char*)&(this->node[i].model),sizeof(int));
			foret.write((char*)&(this->node[i].imax),sizeof(int));
		}
		int ns=this->numset.size();
		foret.write((char*)&ns,sizeof(int));
		for (int i=0;i<ns;i++) {foret.write((char*)&(this->numset[i]),sizeof(int));}
		
	}
	
	void TreeC::lifich(ifstream& foret) {
		int n,m,ns;
		double x;
		bool t;
		foret.read((char*)&m,sizeof(int));
		this->node = vector<NodeRC>(m);
		for (int i=0;i<m;i++) {
			foret.read((char*)&x,sizeof(double));this->node[i].cutval=x;
			foret.read((char*)&n,sizeof(int));this->node[i].filsG=n;
			foret.read((char*)&n,sizeof(int));this->node[i].filsD=n;
			foret.read((char*)&t,sizeof(bool));this->node[i].terminal=t;
			if (this->node[i].terminal) foret.read((char*)&n,sizeof(int));this->node[i].model=n;
			foret.read((char*)&n,sizeof(int));this->node[i].imax=n;
		}
		foret.read((char*)&ns,sizeof(int));
		this->nsets = ns;
		this->numset = vector<int>(ns);
		for (int i=0;i<ns;i++) {foret.read((char*)&n,sizeof(int));this->numset[i]=n;}
		sort(&this->numset[0],&this->numset[nsets]);
	}
	
	void RFC::ecrifich(string nomfi) {
		int n,m;
		ofstream foret(nomfi.c_str(),ios::binary);
		foret.write((char*)&(this->ntrees),sizeof(int));
		foret.write((char*)&(this->ntot),sizeof(int));
		foret.write((char*)&(this->nstat),sizeof(int));
		foret.write((char*)&(this->nmodel),sizeof(int));
		foret.write((char*)&(this->nvar),sizeof(int));
		foret.write((char*)&(this->nsel),sizeof(int));
		foret.write((char*)&(this->nstatclass),sizeof(int));
		foret.write((char*)&(this->nbootsamp),sizeof(int));
		
		n=this->statobs.size();
		foret.write((char*)&n,sizeof(int));
		if (n>0) for (int i=0;i<n;i++) foret.write((char*)&(this->statobs[i]),sizeof(double));		
		
		n=this->stat.size();
		foret.write((char*)&n,sizeof(int));
		if (n>0) for (int i=0;i<n;i++) {
			m=this->stat[i].size();
			foret.write((char*)&m,sizeof(int));
			if (m>0) for (int j=0;j<m;j++) foret.write((char*)&(this->stat[i][j]),sizeof(double));
		}
		
		n=this->tree.size();
		foret.write((char*)&n,sizeof(int));
		if (n>0) for (int i=0;i<n;i++) this->tree[i].ecrifich(foret);
		
		n=this->nimportance.size();
		foret.write((char*)&n,sizeof(int));
		if (n>0) for (int i=0;i<n;i++) {
			m=this->nimportance[i].size();
			foret.write((char*)&m,sizeof(int));
			if (m>0) for (int j=0;j<m;j++) foret.write((char*)&(this->nimportance[i][j]),sizeof(int));
		}

		n=this->importance.size();
		foret.write((char*)&n,sizeof(int));
		if (n>0) for (int i=0;i<n;i++) {
			m=this->importance[i].size();
			foret.write((char*)&m,sizeof(int));
			if (m>0) for (int j=0;j<m;j++) foret.write((char*)&(this->importance[i][j]),sizeof(double));
		}

		foret.close();
	}
	
	
	void lifich(string nomfi) {
		cout<<"\ndebut de lifich\n";
		int n,m,k;
		double x;
		ifstream foret(nomfi.c_str(),ios::binary);
		rfmin = new RFC;
		foret.read((char*)&n,sizeof(int));rfmin->ntrees=n;
		foret.read((char*)&n,sizeof(int));rfmin->ntot=n;
		foret.read((char*)&n,sizeof(int));rfmin->nstat=n;
		foret.read((char*)&n,sizeof(int));rfmin->nmodel=n;
		foret.read((char*)&n,sizeof(int));rfmin->nvar=n;
		foret.read((char*)&n,sizeof(int));rfmin->nsel=n;
		foret.read((char*)&n,sizeof(int));rfmin->nstatclass=n;
		foret.read((char*)&n,sizeof(int));rfmin->nbootsamp=n;
		
		foret.read((char*)&n,sizeof(int));
		rfmin->statobs = vector<double>(n);
		for (int i=0;i<n;i++) {foret.read((char*)&x,sizeof(double));rfmin->statobs[i]=x;}
		
		foret.read((char*)&n,sizeof(int));
		rfmin->stat = vector<vector< double > >(n);
		for (int i=0;i<n;i++) {
			foret.read((char*)&m,sizeof(int));
			rfmin->stat[i] = vector<double>(m);
			for (int j=0;j<m;j++) {foret.read((char*)&x,sizeof(double));rfmin->stat[i][j]=x;}
		}
		
		foret.read((char*)&n,sizeof(int));
		rfmin->tree = vector<TreeC>(n);
		for (int i=0;i<n;i++) {rfmin->tree[i].lifich(foret);cout<<"lecture de l'arbre "<<i+1<<"\r";fflush(stdout);}
		
		foret.read((char*)&n,sizeof(int));
		rfmin->nimportance = vector<vector< int > >(n);
		for (int i=0;i<n;i++) {
			foret.read((char*)&m,sizeof(int));
			rfmin->nimportance[i] = vector<int>(m);
			for (int j=0;j<m;j++) {foret.read((char*)&k,sizeof(int));rfmin->nimportance[i][j]=k;}
		}
		
		foret.read((char*)&n,sizeof(int));
		rfmin->importance = vector<vector< double > >(n);
		for (int i=0;i<n;i++) {
			foret.read((char*)&m,sizeof(int));
			rfmin->importance[i] = vector<double>(m);
			for (int j=0;j<m;j++) {foret.read((char*)&x,sizeof(double));rfmin->importance[i][j]=x;}
		}
		
		
		foret.close();
		cout<<"\nfin de lifich\n";
	}
	
	void dorandfor(string opt,int seed) {
		cout<<"\nDébut de dorandfor\n";
		double duree; clock_t debut;
		debut = clock();
        time_t rawtime;
        time ( &rawtime );
        struct tm * timeinfo;
        timeinfo = localtime ( &rawtime );
		int nrecpos,repmin;
		double ppe,ppemin;
		string nomfiresult,s0,s1;
		nomfiresult = path + ident + "_randomforest.txt";
		progressfilename = path + ident + "_progress.txt";
		cout<<"options : "<<opt<<"\n";
		vector<string> ss;
		vector<string> ss1;  
		splitwords(opt,";",ss);
		rf.ntrees = 500;
		int sssize=(int)ss.size();
		for (int i=0;i<sssize;i++) {
            s0=ss[i].substr(0,2);
            s1=ss[i].substr(2);
			if (s0=="s:"){
				splitwords(s1,",",ss1);rt.nscenchoisi=ss1.size();
                rt.scenchoisi = new int[rt.nscenchoisi];
                for (int j=0;j<rt.nscenchoisi;j++) rt.scenchoisi[j] = atoi(ss1[j].c_str());
                nrecpos=0;for (int j=0;j<rt.nscenchoisi;j++) nrecpos +=rt.nrecscen[rt.scenchoisi[j]-1];
				cout <<"scenarios à tester : ";
                for (int j=0;j<rt.nscenchoisi;j++) {cout<<rt.scenchoisi[j]; if (j<rt.nscenchoisi-1) cout <<",";}cout<<"\n";
				ss1.clear();
			} else if (s0=="n:") {
                rf.nsets=atoi(s1.c_str());
                cout<<"nombre total de jeux de données considérés (tous scénarios confondus)= "<<rf.nsets<<"\n";
            } else if (s0=="d:") {
                LD=true;
                cout<<"ajout des linear discriminant scores\n";
            } else if (s0=="t:") {
                rf.ntrees =atoi(s1.c_str());
                cout<<"nombre d'arbres ' = "<<rf.ntrees<<"\n";
            } else if (s0=="b:") {
				splitwords(s1,",",ss1);rf.nbootsamp=ss1.size();
                rf.bootsamp = vector <int>(rf.nbootsamp);
                for (int j=0;j<rf.nbootsamp;j++) rf.bootsamp[j] = atoi(ss1[j].c_str());
				cout <<"taille des échantillons pour bootstrap : ";
                for (int j=0;j<rf.nbootsamp;j++) {cout<<rf.bootsamp[j]; if (j<rf.nbootsamp-1) cout <<",";}cout<<"\n";
				ss1.clear();
            } else if (s0=="k:") {
				rf.nvar=atoi(s1.c_str());
				flagk=true;
            } else if (s0=="o:") {
				rf.nstatclass=atoi(s1.c_str());
				cout<<"affichage des "<<rf.nstatclass<<" statistiques les plus efficaces\n";
				flago=true;
			}
				
		}
		cout<<"fin de l'analyse des options\n\n";
		rf.readstat(LD);
		if (not flagk)rf.nvar=(int)sqrt(rf.nstat);
		if (not flago)rf.nstatclass=30;
		fout.open(nomfiresult.c_str());
        fout<<"DIYABC :                 Search of the best scenario through Random Forest                         "<<asctime(timeinfo)<<"\n";
        fout<<"Data file       : "<<header.datafilename<<"\n";
        fout<<"Reference table : "<<rt.filename<<"\n";
        fout<<"Total number of simulated data sets : "<<rf.nsets<<"\n";
        if (LD) fout<<"Linear discriminant scores have been added to summary statistics\n";
        fout<<"The forest includes "<<rf.ntrees<<" trees\n";
        fout<<"Calibration has been performed with trees containing ";
        for (int j=0;j<rf.nbootsamp;j++) {fout<<rf.bootsamp[j]; if (j<rf.nbootsamp-1) fout <<",";} fout<<" data sets\n";
        fout<<"Total number of summary statistics : "<<rf.nstat<<"    Number drawn at random at each node : "<<rf.nvar<<"\n\n";
		fout<<"*************************Calibration**************************\n\n";
		ppemin = 1.0;
		rf.dimimportance();
		for (int rep=0;rep<rf.nbootsamp;rep++){
			rf.nsel = rf.bootsamp[rep];
			cout<<"\nrep="<<rep+1<<"    nsel="<<rf.nsel<<"\n";
			cout<<rf.nstat<<" statistiques résumées en tout. Tirage de "<<rf.nvar<<" stat à chaque noeud\n";
			cout<<rf.ntrees<<" arbres à construire avec "<<rf.nsel<<" enregistrements par arbre\n";
			//fout<<rf.nstat<<" statistiques résumées en tout. Tirage de "<<rf.nvar<<" stat à chaque noeud\n";
			//fout<<rf.ntrees<<" arbres à construire avec "<<rf.nsel<<" enregistrements par arbre\n";
			fout<<"-----------Trees with "<<rf.bootsamp[rep]<<" data sets-------------\n";
			rf.growtrees(seed,rep);cout<<"......c'est fait.\n";
			rf.infermodel();
			ppe=rf.training_accuracy();
			if (ppe<=ppemin-0.001) {
				ppemin=ppe;repmin=rep;rf.ecrifich(nomrfmin);
			}
			for (int i=0;i<rf.ntrees;i++) {rf.tree[i].deletree();}
			rf.tree.clear();
			fout<<"-----------------------------------------------------------\n";
		}
		cout<<"\nRésultat de la calibration avec "<<rf.nsets<<" : nsel="<<rf.bootsamp[repmin]<<"  (prior predictive error = "<<ppemin<<")\n";
		fout<<"\nResults of calibration with  "<<rf.nsets<<" overall data sets : \n";
		fout<<"Minimum prior predictive error = "<<ppemin<<" obtained with trees of "<<rf.bootsamp[repmin]<<" data sets\n";
		lifich(nomrfmin);cout<<"apres lecture de la foret sélectionnée\n";
		calposterior(seed);
		var_importance3(repmin);
		
		duree=walltime(debut);
		cout<<"\n\ndurée totale ="<<TimeToStr(duree)<<"\n";
		fout<<"\n\nTotal duration ="<<TimeToStr(duree)<<"\n";
		fout.close();
	}
