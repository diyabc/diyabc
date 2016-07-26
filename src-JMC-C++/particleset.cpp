/*
 *
 * particleset.cppF
 *
 *  Created on: 23 sept. 2010
 *      Author: cornuet
 */

#include <string>
#include <vector>
#include <cstdlib>

#ifdef _OPENMP
#include <omp.h>
#else
#define omp_get_thread_num() 0
#endif

#include "header.hpp"
#include "reftable.hpp"
#include "particuleC.hpp"
#include "particleset.hpp"
#include "mesutils.hpp"

using namespace std;

extern vector<enregC> enreg;
extern long double** phistar;
extern long double** phistarOK;
extern int debuglevel;
extern ReftableC rt;
extern string scurfile;
extern DataC dataobs, datasim;
extern HeaderC header;
extern vector<LocusGroupC> groupe;
extern vector<ScenarioC> scenario;

long int ntentes, naccept;


/**
 * Structure ParticleSet : remplissage de la partie groupe de la particule p
 */
void ParticleSetC::setgroup(int p) {
	int ngr;
	ngr = groupe.size();
	if (debuglevel == 5) cout << "Dans setgroup ngr=" << ngr << "\n";
	this->particule[p].ngr = ngr;
	if (not this->particule[p].grouplist.empty()) this->particule[p].grouplist.clear();
	this->particule[p].grouplist = vector<LocusGroupC>(ngr);
	for (int gr = 1; gr < ngr; gr++) this->particule[p].grouplist[gr] = groupe[gr];
	//cout<<"Dans Particleset setgroup("<<p<<") groupe[1].mutmoy="<<groupe[1].mutmoy<<"  "<<this->particule[p].grouplist[1].mutmoy <<"\n";
	//cout<<"Dans Particleset setgroup("<<p<<") groupe[1].Pmoy="<<groupe[1].Pmoy<<"  "<<this->particule[p].grouplist[1].Pmoy<<"\n";
	//cout<<"Dans Particleset setgroup("<<p<<") groupe[1].snimoy="<<groupe[1].snimoy<<"  "<<this->particule[p].grouplist[1].snimoy<<"\n";
	for (int gr = 1; gr < ngr; gr++) {
		//particuleobs.grouplist[gr] = groupe[gr];
		for (int i = 0; i < groupe[gr].nstatsnp; i++) {
			this->particule[p].grouplist[gr].sumstatsnp[i].x = vector<long double>(groupe[gr].nloc);
			this->particule[p].grouplist[gr].sumstatsnp[i].w = vector<long double>(groupe[gr].nloc);
		}
	}
}

/**
 * Structure ParticleSet : remplissage de la partie locus de la particule p
 */
void ParticleSetC::setloci(int p) {
	int kmoy;
	this->particule[p].nloc = dataobs.nloc;//cout<<"nloc="<<dataobs.nloc<<"\n";
	this->particule[p].nsample = dataobs.nsample;
	if (not particule[p].locuslist.empty()) particule[p].locuslist.clear();
	this->particule[p].locuslist = vector<LocusC>(dataobs.nloc);
	for (int kloc = 0; kloc < dataobs.nloc; kloc++) {
		this->particule[p].locuslist[kloc].nsample = this->particule[p].nsample;
		this->particule[p].locuslist[kloc].coeffcoal = dataobs.locus[kloc].coeffcoal;
		this->particule[p].locuslist[kloc].type = dataobs.locus[kloc].type;
		this->particule[p].locuslist[kloc].groupe = dataobs.locus[kloc].groupe;
		if (dataobs.locus[kloc].type < 5) {
			kmoy = (dataobs.locus[kloc].maxi + dataobs.locus[kloc].mini) / 2;
			this->particule[p].locuslist[kloc].kmin = kmoy - ((dataobs.locus[kloc].motif_range / 2) - 1) * dataobs.locus[kloc].motif_size;
			this->particule[p].locuslist[kloc].kmax = this->particule[p].locuslist[kloc].kmin + (dataobs.locus[kloc].motif_range - 1) * dataobs.locus[kloc].motif_size;
			this->particule[p].locuslist[kloc].motif_size = dataobs.locus[kloc].motif_size;
			this->particule[p].locuslist[kloc].motif_range = dataobs.locus[kloc].motif_range;
			//this->particule[p].locuslist[kloc].mut_rate = dataobs.locus[kloc].mut_rate;
			//this->particule[p].locuslist[kloc].Pgeom = dataobs.locus[kloc].Pgeom;
			//this->particule[p].locuslist[kloc].sni_rate = dataobs.locus[kloc].sni_rate;
		}
		else if (dataobs.locus[kloc].type < 10) {
			this->particule[p].locuslist[kloc].dnalength = dataobs.locus[kloc].dnalength;
			this->particule[p].locuslist[kloc].kmin = 0;
			this->particule[p].locuslist[kloc].kmax = 0;
			//this->particule[p].locuslist[kloc].k2 =  dataobs.locus[kloc].k2;
			this->particule[p].locuslist[kloc].pi_A = dataobs.locus[kloc].pi_A;
			this->particule[p].locuslist[kloc].pi_C = dataobs.locus[kloc].pi_C;
			this->particule[p].locuslist[kloc].pi_G = dataobs.locus[kloc].pi_G;
			this->particule[p].locuslist[kloc].pi_T = dataobs.locus[kloc].pi_T;
			this->particule[p].locuslist[kloc].mutsit.resize(dataobs.locus[kloc].dnalength);
			for (int i = 0; i < this->particule[p].locuslist[kloc].dnalength; i++) this->particule[p].locuslist[kloc].mutsit[i] = dataobs.locus[kloc].mutsit[i];
			//std::cout <<"\n";
		}
		else if (dataobs.locus[kloc].type < 15) {
			//cout << "dans setloci locus "<<kloc<<"\n";
			this->particule[p].locuslist[kloc].samplesize = vector<int>(this->particule[p].nsample);
			for (int ech = 0; ech < this->particule[p].nsample; ech++) {
				this->particule[p].locuslist[kloc].samplesize[ech] = dataobs.locus[kloc].samplesize[ech];
			}
			//cout<<"apres samplesize\n";	
			this->particule[p].locuslist[kloc].ploidie = vector<vector<short int>>(this->particule[p].nsample);
			for (int ech = 0; ech < this->particule[p].nsample; ech++) {
				this->particule[p].locuslist[kloc].ploidie[ech] = vector<short int>(dataobs.nind[ech]);
				for (int i = 0; i < dataobs.nind[ech]; i++)
					this->particule[p].locuslist[kloc].ploidie[ech][i] = dataobs.locus[kloc].ploidie[ech][i];
			}
			//cout<<"apres ploidie\n";
			if (dataobs.locus[kloc].type == 13) {
				if (this->particule[p].samplesizeY.empty()) this->particule[p].samplesizeY = vector<int>(this->particule[p].nsample);
				for (int ech = 0; ech < this->particule[p].nsample; ech++) {
					if (this->particule[p].samplesizeY[ech] < this->particule[p].locuslist[kloc].samplesize[ech]) {
						this->particule[p].samplesizeY[ech] = this->particule[p].locuslist[kloc].samplesize[ech];
					}
				}
			}
			if (dataobs.locus[kloc].type == 14) {
				if (this->particule[p].samplesizeM.empty()) this->particule[p].samplesizeM = vector<int>(this->particule[p].nsample);
				for (int ech = 0; ech < this->particule[p].nsample; ech++) {
					if (this->particule[p].samplesizeM[ech] < this->particule[p].locuslist[kloc].samplesize[ech]) {
						this->particule[p].samplesizeM[ech] = this->particule[p].locuslist[kloc].samplesize[ech];
					}
				}
			}
		}
		else if (dataobs.locus[kloc].type < 20) {
			this->particule[p].locuslist[kloc].samplesize = vector<int>(this->particule[p].nsample);
			this->particule[p].locuslist[kloc].nreads = vector<int>(this->particule[p].nsample);
			this->particule[p].locuslist[kloc].nreads1 = vector<int>(this->particule[p].nsample);
			//if (p<10) cout<<"locus "<<kloc<<" dans setloci nsample ="<<this->particule[p].nsample<<"\n";
			for (int ech = 0; ech < this->particule[p].nsample; ech++) {
				this->particule[p].locuslist[kloc].samplesize[ech] = dataobs.locus[kloc].samplesize[ech];
				this->particule[p].locuslist[kloc].nreads[ech] = dataobs.locus[kloc].nreads[ech];
				//if (p<10) cout<<"samplesize["<<ech<<"]="<<this->particule[p].locuslist[kloc].samplesize[ech]<<"   nreads="<<this->particule[p].locuslist[kloc].nreads[ech]<<"\n";
			}
			//cout<<"locus "<<kloc<<"  apres samplesize\n";		
		}
	}
}

/**
 * Structure ParticleSet : remplissage de la partie scenarios de la particule p
 */
void ParticleSetC::setscenarios(int p, bool simulfile) {
	//cout<<"particuleset::setscenarios  nscenarios = "<<scenario.size()<<"   particule "<<p<<"\n";
	this->particule[p].nscenarios = scenario.size();
	//cout<<"particule["<<p<<"].nscenarios = "<<this->particule[p].nscenarios<<"\n";
	this->particule[p].drawuntil = header.drawuntil;
	if (not particule[p].scenario.empty()) particule[p].scenario.clear();
	this->particule[p].scenario = vector<ScenarioC>(scenario.size());
	int imax = (int)scenario.size();
	for (int i = 0; i < imax; i++) this->particule[p].scenario[i] = scenario[i];
	this->particule[p].scen = header.scen;
	if (simulfile) {
		int ip;
		for (int i = 0; i < imax; i++) {
			for (int ievent = 0; ievent < this->particule[p].scenario[i].nevent; ievent++) {
				this->particule[p].scenario[i].event[ievent].ecris();
				//cout<<"ievent="<<ievent<<"   <"<<this->particule[p].scenario[i].event[ievent].time<<">   {"<<this->particule[p].scenario[i].event[ievent].stime <<"}\n";
				if (this->particule[p].scenario[i].event[ievent].stime.length() > 0) {
					if (atoi(this->particule[p].scenario[i].event[ievent].stime.c_str()) == 0) {
						this->particule[p].scenario[i].event[ievent].time = this->particule[p].getvalue(this->particule[p].scenario[i].event[ievent].stime);
						this->particule[p].scen.event[ievent].time = this->particule[p].scenario[i].event[ievent].time;
						this->particule[p].scen.event[ievent].stime = this->particule[p].scenario[i].event[ievent].stime;
						/*
												for (ip=0;ip<this->particule[p].scenario[i].nparam;ip++) {
													if (this->particule[p].scenario[i].event[ievent].stime==this->particule[p].scenario[i].histparam[ip].name) break;
												}
												this->particule[p].scenario[i].event[ievent].time=this->particule[p].scenario[i].histparam[ip].value;
												this->particule[p].scenario[i].event[ievent].stime=IntToString(this->particule[p].scenario[i].event[ievent].time);
												this->particule[p].scen.event[ievent].time=this->particule[p].scenario[i].event[ievent].time;
												this->particule[p].scen.event[ievent].stime=this->particule[p].scenario[i].event[ievent].stime;
						*/
					}
				}
				cout << "ievent=" << ievent << "   <" << this->particule[p].scenario[i].event[ievent].time << ">   {" << this->particule[p].scenario[i].event[ievent].stime << "}\n";
				//cout<<"event "<<ievent<<"  action="<<this->particule[p].scenario[i].event[ievent].action<<"  time="<<this->particule[p].scenario[i].event[ievent].time<<"\n";
			}
		}
	}
	//this->particule[p].scen.ecris(simulfile);
	//cout<<"setscenarios nloc="<<this->particule[p].nloc<<"\n";
}

/**
 * Structure ParticleSet : réinitialisation de la particule p
 */
void ParticleSetC::resetparticle(int p) {
	for (int gr = 0; gr < this->particule[p].ngr; gr++) {
		if (groupe[gr].type == 0) { //MICROSAT
			//cout <<"MICROSAT "<<p<<"\n";
			this->particule[p].grouplist[gr].mutmoy = groupe[gr].mutmoy;
			this->particule[p].grouplist[gr].Pmoy = groupe[gr].Pmoy;
			this->particule[p].grouplist[gr].snimoy = groupe[gr].snimoy;
		}
		else if (groupe[gr].type == 1) { //SEQUENCES
			//cout<<"SEQUENCE\n";
			int kloc;
			for (int iloc = 0; iloc < this->particule[p].grouplist[gr].nloc; iloc++) {
				kloc = this->particule[p].grouplist[gr].loc[iloc];
				if (not this->particule[p].locuslist[kloc].mutsit.empty())this->particule[p].locuslist[kloc].mutsit.clear();
			}
			this->particule[p].grouplist[gr].musmoy = groupe[gr].musmoy; //musmoy
			if (groupe[gr].mutmod > 0) {
				this->particule[p].grouplist[gr].k1moy = groupe[gr].k1moy; //k1moy
			}
			if (groupe[gr].mutmod > 2) {
				this->particule[p].grouplist[gr].k2moy = groupe[gr].k2moy; //k2moy
			}
		}
		else if (groupe[gr].type == 2) { //SNP
			for (int i = 0; i < this->particule[p].grouplist[gr].nstatsnp; i++) {
				this->particule[p].grouplist[gr].sumstatsnp[i].defined = false;
			}
		}
		//for (int loc=0;loc<this->nloc;loc++) this->particule[p].libere(false);
	}
}

void ParticleSetC::init_particule(int npart, bool dnatrue, bool simulfile, int seed) {
	if (debuglevel == 5) cout << "debut de init_particule npart=" << npart << "\n";
	this->particule.resize(npart);
	//this->particule = vector<ParticleC>(npart);
	if (debuglevel == 5) cout << "après dimensionnement du vecteur de particules\n";
	for (int p = 0; p < npart; p++) {
		if (debuglevel == 5) cout << "avant set particule " << p << "\n";
		this->particule[p].dnatrue = dnatrue;
		this->particule[p].maf = dataobs.maf;
		this->particule[p].mrc = dataobs.mrc;
		if (debuglevel == 5) cout << "dnatrue\n";
		this->setgroup(p);
		if (debuglevel == 5) cout << "setgroup\n";
		this->setloci(p);
		if (debuglevel == 5) cout << "setloci   nloc=" << this->particule[p].nloc << "\n";
		this->setscenarios(p, simulfile);
		if (debuglevel == 5) cout << "setscenarios\n";
		this->particule[p].mw.randinit(p, seed);
		if (debuglevel == 5) cout << "                    apres init_particule 					particule " << p << "\n";
	}

}

void ParticleSetC::free_particule(int npart) {
	cout << "debut de free_particule\n";
	int n, imax;
	for (int p = 0; p < npart; p++) {
		//Partie grouplist
		int grmax = (int)groupe.size();
		for (int gr = 1; gr < grmax; gr++) {
			//particuleobs.grouplist[gr] = groupe[gr];
			for (int i = 0; i < groupe[gr].nstatsnp; i++) {
				this->particule[p].grouplist[gr].sumstatsnp[i].x.clear();
				this->particule[p].grouplist[gr].sumstatsnp[i].w.clear();
			}
			if (not this->particule[p].grouplist[gr].sumstatsnp.empty()) this->particule[p].grouplist[gr].sumstatsnp.clear();
			if (not this->particule[p].grouplist[gr].sumstat.empty()) this->particule[p].grouplist[gr].sumstat.clear();
			if (not this->particule[p].grouplist[gr].loc.empty()) this->particule[p].grouplist[gr].loc.clear();
		}
		if (not this->particule[p].grouplist.empty()) this->particule[p].grouplist.clear();
		//Partie locuslist
		for (int kloc = 0; kloc < dataobs.nloc; kloc++) {
			if (not this->particule[p].locuslist[kloc].mutsit.empty()) this->particule[p].locuslist[kloc].mutsit.clear();
			if (not this->particule[p].locuslist[kloc].sitmut.empty()) this->particule[p].locuslist[kloc].sitmut.clear();
			if (not this->particule[p].locuslist[kloc].sitmut2.empty()) this->particule[p].locuslist[kloc].sitmut2.clear();
			if (not this->particule[p].locuslist[kloc].tabsit.empty()) this->particule[p].locuslist[kloc].tabsit.clear();
			if (not this->particule[p].locuslist[kloc].freq.empty()) {
				imax = (int)this->particule[p].locuslist[kloc].freq.size();
				for (int i = 0; i < imax; i++) {
					if (not this->particule[p].locuslist[kloc].freq[i].empty()) this->particule[p].locuslist[kloc].freq[i].clear();
				}
			}
			if (not this->particule[p].locuslist[kloc].haplomic.empty()) {
				imax = (int)this->particule[p].locuslist[kloc].haplomic.size();
				for (int i = 0; i < imax; i++) {
					if (not this->particule[p].locuslist[kloc].haplomic[i].empty()) this->particule[p].locuslist[kloc].haplomic[i].clear();
				}
			}
			if (not this->particule[p].locuslist[kloc].haplodna.empty()) {
				imax = (int)this->particule[p].locuslist[kloc].haplodna.size();
				for (int i = 0; i < imax; i++) {
					if (not this->particule[p].locuslist[kloc].haplodna[i].empty()) this->particule[p].locuslist[kloc].haplodna[i].clear();
				}
			}
			if (not this->particule[p].locuslist[kloc].haplodnavar.empty()) {
				imax = (int)this->particule[p].locuslist[kloc].haplodnavar.size();
				for (int i = 0; i < imax; i++) {
					if (not this->particule[p].locuslist[kloc].haplodnavar[i].empty()) this->particule[p].locuslist[kloc].haplodnavar[i].clear();
				}
			}
			if (not this->particule[p].locuslist[kloc].haplosnp.empty()) {
				imax = (int)this->particule[p].locuslist[kloc].haplosnp.size();
				for (int i = 0; i < imax; i++) {
					if (not this->particule[p].locuslist[kloc].haplosnp[i].empty()) this->particule[p].locuslist[kloc].haplosnp[i].clear();
				}
			}
			if (not this->particule[p].locuslist[kloc].ploidie.empty()) {
				imax = (int)this->particule[p].locuslist[kloc].ploidie.size();
				for (int i = 0; i < imax; i++) {
					if (not this->particule[p].locuslist[kloc].ploidie[i].empty()) this->particule[p].locuslist[kloc].ploidie[i].clear();
				}
			}
		}
		if (not this->particule[p].locuslist.empty()) this->particule[p].locuslist.clear();
		//Partie scenario
		cout << "partie scenario  particule " << p << "\n";
		imax = (int)this->particule[p].scenario.size();
		for (int i = 0; i < imax; i++) {
			cout << "        scenario " << i << "\n";
			if (not this->particule[p].scenario[i].paramvar.empty()) this->particule[p].scenario[i].paramvar.clear();
			cout << "            paramvar\n";
			if (not this->particule[p].scenario[i].time_sample.empty()) this->particule[p].scenario[i].time_sample.clear();
			cout << "            time_sample\n";
			if (not this->particule[p].scenario[i].stime_sample.empty()) this->particule[p].scenario[i].stime_sample.clear();
			cout << "            stime_sample\n";
			if (not this->particule[p].scenario[i].event.empty()) this->particule[p].scenario[i].event.clear();
			cout << "            event\n";
			if (not this->particule[p].scenario[i].ne0.empty()) this->particule[p].scenario[i].ne0.clear();
			cout << "            ne0\n";
			if (not this->particule[p].scenario[i].histparam.empty()) this->particule[p].scenario[i].histparam.clear();
			cout << "            histparam\n";
			n = this->particule[p].scenario[i].condition.size();
			cout << "                    nconditions=" << n << "\n";
			if (not this->particule[p].scenario[i].condition.empty()) this->particule[p].scenario[i].condition.resize(0);
			cout << "            condition\n";
		}
		cout << "avant scenario.clear\n";
		this->particule[p].scenario.resize(3);
		cout << "apres scenario.resize(3)\n";
		//this->particule[p].scenario.clear();
		cout << "apres scenario.clear()\n";
	}
	cout << "avant this->particule.clear()\n";
	this->particule.resize(0);
	cout << "fin de free_particule\n";

}

/**
 * Structure ParticleSet : simulation de particules de paramètres fixés
 */
void ParticleSetC::dosimulstat(int debut, int npartic, bool dnatrue, bool multithread, int numscen, int seed, double** stat) {
	int scen = numscen - 1;
	this->npart = npartic;
	int k, ii, gr, ipart, nstat, nparam = scenario[scen].nparam;
	//string *ss;
	vector<int> sOK;
	sOK = vector<int>(npartic);
	bool simulfile = false;
	if (debut < 1) init_particule(npartic, dnatrue, simulfile, seed);
	else {
		for (int p = 0; p < npartic; p++) {
			this->resetparticle(p);
		}
	}
	for (int p = 0; p < npartic; p++) {
		k = debut + p;
		if (debuglevel == 5) cout << "p=" << p << "    k=" << k << "\n";
		ii = 0;
		for (int i = 0; i < nparam; i++) {
			if (not scenario[scen].histparam[i].prior.constant) {
				this->particule[p].scenario[scen].histparam[i].value = phistarOK[k][ii];
				this->particule[p].scenario[scen].histparam[i].prior.fixed = true;
				ii++;
			}
			else this->particule[p].scenario[scen].histparam[i].value = scenario[scen].histparam[i].prior.mini;
			if (debuglevel == 5) {
				cout << scenario[scen].histparam[i].name << "=" << this->particule[p].scenario[scen].histparam[i].value << "   ";
				cout << "phistar[" << k << "][" << ii - 1 << "]=" << phistarOK[k][ii - 1] << "\n";
			}
		}

		if (debuglevel == 5) cout << "\n";
		//cout<<"apres la copie des paramètres historiques\n";
		int grmax = (int)groupe.size();
		for (gr = 1; gr < grmax; gr++)
			if (groupe[gr].type == 0) { //MICROSAT
				//cout <<"MICROSAT particule "<<p<<"   gr="<<gr<<"   k="<<k<<"   ii="<<ii<<"\n";
				//cout<<"phistarOK["<<k<<"]["<<ii<<"]="<<phistarOK[k][ii]<<"\n";
				if (not groupe[gr].priormutmoy.constant) {
					//cout<<" priormutmoy non constant "<<phistarOK[k][ii]<<" \n";
					this->particule[p].grouplist[gr].mutmoy = phistarOK[k][ii];
					ii++;
					this->particule[p].grouplist[gr].priormutmoy.fixed = true;
				}
				//cout<<this->particule[p].grouplist[gr].mutmoy<<"\n";
				//cout<<"phistarOK["<<k<<"]["<<ii<<"]="<<phistarOK[k][ii]<<"    ";
				if (not groupe[gr].priorPmoy.constant) {
					//cout<<" priorPmoy non constant "<<phistarOK[k][ii]<<" \n";
					this->particule[p].grouplist[gr].Pmoy = phistarOK[k][ii];
					ii++;
					this->particule[p].grouplist[gr].priorPmoy.fixed = true;
				}
				//cout<<this->particule[p].grouplist[gr].Pmoy<<"\n";
				//cout<<"phistarOK["<<k<<"]["<<ii<<"]="<<phistarOK[k][ii]<<"    ";
				if (not groupe[gr].priorsnimoy.constant) {
					//cout<<" priorsnimoy non constant "<<phistarOK[k][ii]<<" \n";
					this->particule[p].grouplist[gr].snimoy = phistarOK[k][ii];
					ii++;
					this->particule[p].grouplist[gr].priorsnimoy.fixed = true;
				}
				//cout<<this->particule[p].grouplist[gr].snimoy<<"\n";
			}
			else if (groupe[gr].type == 1) { //SEQUENCES
				//cout<<"SEQUENCE\n";
				if (not groupe[gr].priormusmoy.constant) { //musmoy
					this->particule[p].grouplist[gr].musmoy = phistarOK[k][ii];
					ii++;
					this->particule[p].grouplist[gr].priormusmoy.fixed = true;
				}
				if (groupe[gr].mutmod > 0) {
					if (not groupe[gr].priork1moy.constant) { //k1moy
						this->particule[p].grouplist[gr].k1moy = phistarOK[k][ii];
						ii++;
						this->particule[p].grouplist[gr].priork1moy.fixed = true;
					}
				}
				if (groupe[gr].mutmod > 2) {
					if (not groupe[gr].priork2moy.constant) { //k2moy
						this->particule[p].grouplist[gr].k2moy = phistarOK[k][ii];
						ii++;
						this->particule[p].grouplist[gr].priork2moy.fixed = true;
					}
				}
			}
		//cout<<"\n";
	}
	//if (multithread) cout<<"calcul multithread\n"; else cout<<"pas de calcul multithread\n";
#pragma omp parallel for shared(sOK) private(gr) if(multithread)
	for (ipart = 0; ipart < npartic; ipart++) {
		if (debuglevel == 5) cout << "avant dosimulpart de la particule " << ipart << "   avec scenario " << numscen << "\n";
		sOK[ipart] = this->particule[ipart].dosimulpart(numscen);
		if (debuglevel == 5) cout << sOK[ipart] << "apres dosimulpart de la particule " << ipart << "\n";
		if (sOK[ipart] == 0) {
			for (gr = 1; gr < this->particule[ipart].ngr; gr++) {
				this->particule[ipart].docalstat(gr);
			}
		}
	}
	cout << "after simulation of " << npart << " particles (scenario " << numscen << ")\n";
	///////////////fin du pragma
	for (int ipart = 0; ipart < npartic; ipart++) {
		if (sOK[ipart] == 0) {
			nstat = 0;
			for (int gr = 1; gr < this->particule[ipart].ngr; gr++) {
				//cout<<"gr="<<gr<<"\n";
				for (int st = 0; st < this->particule[ipart].grouplist[gr].nstat; st++) {
					//cout<<"ipart="<<ipart<<"   gr="<<gr<<"   st="<<st<<"\n";
					/////////////// filtrage des Aml et aml début
					if (this->particule[ipart].grouplist[gr].sumstat[st].val == -9999.0) {
						if ((this->particule[ipart].grouplist[gr].sumstat[st].cat == 12)or (this->particule[ipart].grouplist[gr].sumstat[st].cat == -14))
							this->particule[ipart].grouplist[gr].sumstat[st].val = 0.5;
					}
					//cout<<"avant la fin du filtrage debut+ipart="<<debut+ipart<<"   nstat="<<nstat<<"\n";
					/////////////// filtrage des Aml et aml fin
					stat[debut + ipart][nstat] = this->particule[ipart].grouplist[gr].sumstat[st].val;
					//cout<<"coucou\n";
					//cout<<this->particule[ipart].grouplist[gr].sumstat[st].val<<" ("<<stat[debut+ipart][nstat]<<")  ";
					nstat++;
				}
				//cout<<"\n";
			}
		}
	}
	//cout<<"avant delete sOK \n";
	//cout<<"sOK size = "<<sOK.size()<<"\n";
	if (not sOK.empty()) sOK.clear();
	if (debuglevel == 5) cout << "fin de dosimulstat\n";
}

/**
 * Structure ParticleSet : simulation de particules de paramètres fixés deuxième version
 */
void ParticleSetC::dosimulstat2(int debut, int npartic, vector<vector<long double>>& phistarOKK, bool dnatrue, bool multithread, int numscen, int seed, double** stat) {
	//cout<<"debut de dosimulstat2\n";
	int scen = numscen - 1;
	this->npart = npartic;
	int k, ii, gr, ipart, nstat, nparam = scenario[scen].nparam;
	//string *ss;
	vector<int> sOK;
	sOK = vector<int>(npartic);
	bool simulfile = false;
	if (debut < 100) {
		//cout<<"avant init_particule\n";
		init_particule(npartic, dnatrue, simulfile, seed);
	}
	else {
		for (int p = 0; p < npartic; p++) {
			this->resetparticle(p);
		}
	}
	for (int p = 0; p < npartic; p++) {
		//k=debut+p;
		k = p;
		if (debuglevel == 5) cout << "p=" << p << "    k=" << k << "\n";
		ii = 0;
		for (int i = 0; i < nparam; i++) {
			if (not scenario[scen].histparam[i].prior.constant) {
				this->particule[p].scenario[scen].histparam[i].value = phistarOKK[k][ii];
				this->particule[p].scenario[scen].histparam[i].prior.fixed = true;
				ii++;
			}
			else this->particule[p].scenario[scen].histparam[i].value = scenario[scen].histparam[i].prior.mini;
			if (debuglevel == 5) {
				cout << scenario[scen].histparam[i].name << "=" << this->particule[p].scenario[scen].histparam[i].value << "   ";
				cout << "phistar[" << k << "][" << ii - 1 << "]=" << phistarOKK[k][ii - 1] << "\n";
			}
		}

		if (debuglevel == 5) cout << "\n";
		//cout<<"apres la copie des paramètres historiques\n";
		int grmax = (int)groupe.size();
		for (gr = 1; gr < grmax; gr++)
			if (groupe[gr].type == 0) { //MICROSAT
				//cout <<"MICROSAT particule "<<p<<"   gr="<<gr<<"   k="<<k<<"   ii="<<ii<<"\n";
				//cout<<"phistarOKK["<<k<<"]["<<ii<<"]="<<phistarOKK[k][ii]<<"\n";
				if (not groupe[gr].priormutmoy.constant) {
					//cout<<" priormutmoy non constant "<<phistarOKK[k][ii]<<" \n";
					this->particule[p].grouplist[gr].mutmoy = phistarOKK[k][ii];
					ii++;
					this->particule[p].grouplist[gr].priormutmoy.fixed = true;
				}
				//cout<<this->particule[p].grouplist[gr].mutmoy<<"\n";
				//cout<<"phistarOKK["<<k<<"]["<<ii<<"]="<<phistarOKK[k][ii]<<"    ";
				if (not groupe[gr].priorPmoy.constant) {
					//cout<<" priorPmoy non constant "<<phistarOKK[k][ii]<<" \n";
					this->particule[p].grouplist[gr].Pmoy = phistarOKK[k][ii];
					ii++;
					this->particule[p].grouplist[gr].priorPmoy.fixed = true;
				}
				//cout<<this->particule[p].grouplist[gr].Pmoy<<"\n";
				//cout<<"phistarOKK["<<k<<"]["<<ii<<"]="<<phistarOKK[k][ii]<<"    ";
				if (not groupe[gr].priorsnimoy.constant) {
					//cout<<" priorsnimoy non constant "<<phistarOKK[k][ii]<<" \n";
					this->particule[p].grouplist[gr].snimoy = phistarOKK[k][ii];
					ii++;
					this->particule[p].grouplist[gr].priorsnimoy.fixed = true;
				}
				//cout<<this->particule[p].grouplist[gr].snimoy<<"\n";
			}
			else if (groupe[gr].type == 1) { //SEQUENCES
				//cout<<"SEQUENCE\n";
				if (not groupe[gr].priormusmoy.constant) { //musmoy
					this->particule[p].grouplist[gr].musmoy = phistarOKK[k][ii];
					ii++;
					this->particule[p].grouplist[gr].priormusmoy.fixed = true;
				}
				if (groupe[gr].mutmod > 0) {
					if (not groupe[gr].priork1moy.constant) { //k1moy
						this->particule[p].grouplist[gr].k1moy = phistarOKK[k][ii];
						ii++;
						this->particule[p].grouplist[gr].priork1moy.fixed = true;
					}
				}
				if (groupe[gr].mutmod > 2) {
					if (not groupe[gr].priork2moy.constant) { //k2moy
						this->particule[p].grouplist[gr].k2moy = phistarOKK[k][ii];
						ii++;
						this->particule[p].grouplist[gr].priork2moy.fixed = true;
					}
				}
			}
		//cout<<"\n";
	}
	//if (multithread) cout<<"calcul multithread\n"; else cout<<"pas de calcul multithread\n";
#pragma omp parallel for shared(sOK) private(gr) if(multithread)
	for (ipart = 0; ipart < npartic; ipart++) {
		if (debuglevel == 5) cout << "avant dosimulpart de la particule " << ipart << "   avec scenario " << numscen << "\n";
		sOK[ipart] = this->particule[ipart].dosimulpart(numscen);
		if (debuglevel == 5) cout << sOK[ipart] << " apres dosimulpart de la particule " << ipart << "\n";
		if (sOK[ipart] == 0) {
			for (gr = 1; gr < this->particule[ipart].ngr; gr++) {
				for (int st = 0; st < this->particule[ipart].grouplist[gr].nstat; st++)this->particule[ipart].grouplist[gr].sumstat[st].val = 0.0;
				this->particule[ipart].docalstat(gr);
			}
		}
	}
	if (debuglevel == 5) cout << "after simulation of " << npart << " particles (scenario " << numscen << ")    ";
	///////////////fin du pragma
	for (int ipart = 0; ipart < npartic; ipart++) {
		if (sOK[ipart] == 0) {
			nstat = 0;
			for (int gr = 1; gr < this->particule[ipart].ngr; gr++) {
				//cout<<"gr="<<gr<<"    stat.size = "<<stat[ipart].size()<<"      nstat = "<<this->particule[ipart].grouplist[gr].nstat<<"\n";
				for (int st = 0; st < this->particule[ipart].grouplist[gr].nstat; st++) {
					//if ((ipart<5)and(nstat<5)) cout<<"ipart="<<ipart<<"   gr="<<gr<<"   st="<<st<<"      ";
					/////////////// filtrage des Aml et aml début
					if (this->particule[ipart].grouplist[gr].sumstat[st].val == -9999.0) {
						if ((this->particule[ipart].grouplist[gr].sumstat[st].cat == 12)or (this->particule[ipart].grouplist[gr].sumstat[st].cat == -14))
							this->particule[ipart].grouplist[gr].sumstat[st].val = 0.5;
					}
					//cout<<"avant la fin du filtrage ipart="<<ipart<<"   nstat="<<nstat<<"\n";
					/////////////// filtrage des Aml et aml fin
					stat[ipart][nstat] = this->particule[ipart].grouplist[gr].sumstat[st].val;
					//cout<<"coucou\n";
					//if ((ipart<5)and(nstat<5)) cout<<this->particule[ipart].grouplist[gr].sumstat[st].val<<" ("<<stat[ipart][nstat]<<")  \n";
					nstat++;
				}
				if ((ipart < 5)and (nstat < 5)) cout << "\n";
			}
		}
	}
	//cout<<"avant delete sOK \n";
	//cout<<"sOK size = "<<sOK.size()<<"\n";
	if (not sOK.empty()) sOK.clear();
	//free_particule(npartic);
	if (debuglevel == 5) cout << "fin de dosimulstat2\n";
}

/**
 * Structure ParticleSet : simulation des particules utilisées pour le model checking
 */
void ParticleSetC::dosimulphistar(int npart, bool dnatrue, bool multithread, bool firsttime, int numscen, int seed, int nsel) {
	int scen = numscen - 1;
	int ii, ipart, gr, nstat, k, nparam = scenario[scen].nparam, np, iscen, categ;
	bool simulfile = false;
	this->npart = npart;
	int ns;
	vector<int> sOK;
	vector<string> ss;
	sOK = vector<int>(npart);
	if (debuglevel == 5) cout << "avant firsttime=" << firsttime << "\n";
	if (firsttime) init_particule(npart, dnatrue, simulfile, seed);
	else {
		for (int p = 0; p < this->npart; p++) {
			this->resetparticle(p);
		}
	}
	if (debuglevel == 5) cout << "apres firsttime\n";
	for (int p = 0; p < this->npart; p++) {
		k = this->particule[0].mw.rand0(nsel);
		if (debuglevel == 55) cout << "p=" << p << "    k=" << k << "\n";
		ii = 0;
		for (int i = 0; i < nparam; i++) {
			if (not scenario[scen].histparam[i].prior.constant) {
				if (debuglevel == 55) cout << "phistar[" << k << "][" << ii << "]=" << phistarOK[k][ii] << "\n";
				this->particule[p].scenario[scen].histparam[i].value = phistarOK[k][ii];
				this->particule[p].scenario[scen].histparam[i].prior.fixed = true;
				ii++;
			}
			else this->particule[p].scenario[scen].histparam[i].value = scenario[scen].histparam[i].prior.mini;
			if (debuglevel == 5) cout << scenario[scen].histparam[i].name << "=" << this->particule[p].scenario[scen].histparam[i].value << "   ";
		}
		if (debuglevel == 5) cout << "\n";
		//cout<<"apres la copie des paramètres historiques\n";
		int grmax = (int)groupe.size();
		for (gr = 1; gr < grmax; gr++)
			if (groupe[gr].type == 0) { //MICROSAT
				//cout <<"MICROSAT "<<p<<"\n";
				//cout<<"phistar["<<k<<"]["<<ii<<"]="<<phistar[k][ii]<<"    \n";
				if (not groupe[gr].priormutmoy.constant) {
					this->particule[p].grouplist[gr].mutmoy = phistarOK[k][ii];
					ii++;
					this->particule[p].grouplist[gr].priormutmoy.fixed = true;
				}
				//cout<<this->particule[p].grouplist[gr].mutmoy<<"\n";
				//cout<<"phistar["<<k<<"]["<<ii<<"]="<<phistar[k][ii]<<"    ";
				if (not groupe[gr].priorPmoy.constant) {
					this->particule[p].grouplist[gr].Pmoy = phistarOK[k][ii];
					ii++;
					this->particule[p].grouplist[gr].priorPmoy.fixed = true;
				}
				//cout<<this->particule[p].grouplist[gr].Pmoy<<"\n";
				//cout<<"phistar["<<k<<"]["<<ii<<"]="<<phistarOK[k][ii]<<"    ";
				if (not groupe[gr].priorsnimoy.constant) {
					this->particule[p].grouplist[gr].snimoy = phistarOK[k][ii];
					ii++;
					this->particule[p].grouplist[gr].priorsnimoy.fixed = true;
				}
				//cout<<this->particule[p].grouplist[gr].snimoy<<"\n";
			}
			else if (groupe[gr].type == 1) { //SEQUENCES
				//cout<<"SEQUENCE\n";
				if (not groupe[gr].priormusmoy.constant) { //musmoy
					this->particule[p].grouplist[gr].musmoy = phistarOK[k][ii];
					ii++;
					this->particule[p].grouplist[gr].priormusmoy.fixed = true;
				}
				if (groupe[gr].mutmod > 0) {
					if (not groupe[gr].priork1moy.constant) { //k1moy
						this->particule[p].grouplist[gr].k1moy = phistarOK[k][ii];
						ii++;
						this->particule[p].grouplist[gr].priork1moy.fixed = true;
					}
				}
				if (groupe[gr].mutmod > 2) {
					if (not groupe[gr].priork2moy.constant) { //k2moy
						this->particule[p].grouplist[gr].k2moy = phistarOK[k][ii];
						ii++;
						this->particule[p].grouplist[gr].priork2moy.fixed = true;
					}
				}
			}
		//cout<<"\n";
	}
	//cout<<"apres le remplissage de "<<npart<<" particules\n";
	//if (num_threads>0) omp_set_num_threads(num_threads);_
#pragma omp parallel for shared(sOK) private(gr) if(multithread)
	for (ipart = 0; ipart < this->npart; ipart++) {
		cout << "avant dosimulpart de la particule " << ipart << "\n";
		sOK[ipart] = this->particule[ipart].dosimulpart(numscen);
		cout << sOK[ipart] << "apres dosimulpart de la particule " << ipart << "\n";
		if (sOK[ipart] == 0) {
			for (gr = 1; gr < this->particule[ipart].ngr; gr++) {
				this->particule[ipart].docalstat(gr);
				cout << "statistiques de la particules : \n";
				for (int st = 0; st < this->particule[ipart].grouplist[gr].nstat; st++) cout << this->particule[ipart].grouplist[gr].sumstat[st].val << "\n";
			}
		}
	}
	cout << "apres la simulation de " << npart << " particules\n";
	///////////////fin du pragma
	for (int ipart = 0; ipart < this->npart; ipart++) {
		if (sOK[ipart] == 0) {
			enreg[ipart].numscen = 1;
			if (this->particule[ipart].nscenarios > 1) {
				enreg[ipart].numscen = this->particule[ipart].scen.number;
			}

			//cout<<"dans particleset\n";
			for (int j = 0; j < this->particule[ipart].scen.nparamvar; j++) {
				enreg[ipart].param[j] = this->particule[ipart].scen.paramvar[j];
				//cout<<this->particule[ipart].scen.paramvar[j]<<"  ("<<enreg[ipart].param[j]<<")     ";
			}
			//cout <<"\n";
			nstat = 0;
			for (int gr = 1; gr < this->particule[ipart].ngr; gr++) {
				for (int st = 0; st < this->particule[ipart].grouplist[gr].nstat; st++) {
					/////////////// filtrage des Aml et aml début
					if (this->particule[ipart].grouplist[gr].sumstat[st].val == -9999.0) {
						if ((this->particule[ipart].grouplist[gr].sumstat[st].cat == 12)or (this->particule[ipart].grouplist[gr].sumstat[st].cat == -14))
							this->particule[ipart].grouplist[gr].sumstat[st].val = 0.5;
					}
					/////////////// filtrage des Aml et aml fin
					enreg[ipart].stat[nstat] = this->particule[ipart].grouplist[gr].sumstat[st].val;
					cout << this->particule[ipart].grouplist[gr].sumstat[st].val << " (" << enreg[ipart].stat[nstat] << ")  ";
					nstat++;
				}
				cout << "\n";
			}
			enreg[ipart].message = "OK";
		}
		else {
			enreg[ipart].message = "A gene genealogy failed in scenario ";
			enreg[ipart].message += IntToString(this->particule[ipart].scen.number);
			enreg[ipart].message += ". Check consistency of this scenario over possible prior ranges of time event parameters.";
		}
	}

	int nph, npm, nn = 0, ccc, pa, ip, iq;
	splitwords(header.entetehist, " ", ss);
	nph = ss.size();
	for (int m = 0; m < nph; m++) {
		ccc = ss[m][0];
		if (not isalnum(ccc)) nn++;
	}
	if (header.entetemut.length() > 10) {
		splitwords(header.entetemut, " ", ss);
		npm = ss.size();
	}
	else npm = 0;
	nph -= nn;
	FILE* pFile;
	bool trouve, trouve2;
	if (firsttime) pFile = fopen(scurfile.c_str(), "w");
	else pFile = fopen(scurfile.c_str(), "a");
	if (firsttime) fprintf(pFile, "%s\n", header.entete.c_str());
	splitwords(header.entete, " ", ss);
	ns = ss.size();
	if (debuglevel == 5) cout << "nph=" << nph << "   npm=" << npm << "   ns=" << ns << "\n";
	np = ns - header.nstat - 1;
	for (int ipart = 0; ipart < this->npart; ipart++) {
		//cout<<"ipart="<<ipart<<"    sOK[ipart]="<<sOK[ipart]<<"   message="<<enreg[ipart].message<<"\n";
		if (sOK[ipart] == 0) {
			//cout<<enreg[ipart].numscen<<"\n";
			fprintf(pFile, "%3d  ", enreg[ipart].numscen);
			iscen = enreg[ipart].numscen - 1;
			//cout<<"scenario "<<enreg[ipart].numscen<<"\n";
			if (debuglevel == 5) cout << header.nparamtot << "   " << rt.nhistparam[iscen] << "   " << np << "\n";
			pa = 0;
			for (int j = 1; j < nph; j++) {
				trouve = false;
				ip = -1;
				while ((not trouve)and (ip < rt.nhistparam[iscen] - 1)) {
					ip++;
					trouve = (ss[j] == rt.histparam[iscen][ip].name);
					//cout<<"->"<<ss[j]<<"<-   ?   ->"<< header.scenario[iscen].histparam[ip].name<<"<-"<<trouve<<"\n";
				}
				if (trouve) {
					trouve2 = false;
					iq = -1;
					while ((not trouve2)and (iq < header.nparamtot - 1)) {
						iq++;
						trouve2 = (ss[j] == header.histparam[iq].name);
					}
					if (trouve2) categ = header.histparam[iq].category;
					else categ = 0;
					if (categ < 2) fprintf(pFile, "  %12.0f", enreg[ipart].param[ip]);
					else fprintf(pFile, "  %12.3f", enreg[ipart].param[ip]);
					pa++;
				}
				else fprintf(pFile, "              ");
			}
			if (debuglevel == 5) cout << "pa=" << pa << "   rt.nparam[iscen]=" << rt.nparam[iscen] << "\n";
			for (int j = 0; j < npm; j++) fprintf(pFile, "  %12.3e", enreg[ipart].param[pa + j]);
			for (int st = 0; st < header.nstat; st++) fprintf(pFile, "  %12.6f", enreg[ipart].stat[st]);
			fprintf(pFile, "\n");
			if (debuglevel == 55) {
				cout << "parametres de la particule" << ipart << "\n";
				for (int j = 0; j < npm + pa; j++) cout << enreg[ipart].param[j] << "   ";
				cout << "\n";
				cout << "stat de la particule " << ipart << "\n";
				for (int st = 0; st < header.nstat; st++) cout << enreg[ipart].stat[st] << "   ";
				cout << "\n\n";
			}
		}
	}
	fclose(pFile);
	sOK.clear();
	//cout<<"fin de dosimulphistar\n";
}

/*void ParticleSetC::libere(int npart) {

}*/


/**
 * Structure ParticleSet : simulation des particules utilisées pour la table de référence, le biais et la confiance
 */
void ParticleSetC::dosimultabref(int npart, bool dnatrue, bool multithread, bool firsttime, int numscen, int seed, int depuis) {
	int ipart, gr, nstat, pa, ip, iscen, np, ns, ccc;
	bool trouve;
	this->npart = npart;
	int* sOK;
	vector<string> ss;
	//cout <<"debut de dosimultabref  nscenarios="<<scenario.size()<<"  ngroupes="<<groupe.size()<<"\n";
	if (debuglevel == 5) cout << "debut de dosimultabref\n";
	fflush(stdin);
	sOK = new int[npart];
	if (debuglevel == 5) cout << "dosimultabref npart = " << npart << "\n";
	if (firsttime) init_particule(npart, dnatrue, false, seed);
	else {
		for (int p = 0; p < this->npart; p++) {
			this->resetparticle(p);
		}
	}
	if (debuglevel == 5) cout << "avant pragma npart = " << this->npart << "\n";
	//if (num_threads>0) omp_set_num_threads(num_threads);_
	//ntentes=0;naccept=0;
#pragma omp parallel for shared(sOK) private(gr) if(multithread)
	for (ipart = 0; ipart < this->npart; ipart++) {
		//cout<<"this->particule["<<ipart<<"].grouplist[1].mutmoy = "<<this->particule[ipart].grouplist[1].mutmoy<<"\n";
		//cout<<"this->particule["<<ipart<<"].grouplist[1].priormutmoy.constant = "<<this->particule[ipart].grouplist[1].priormutmoy.constant<<"\n";
		if (debuglevel == 5) cout << "avant dosimulpart de la particule " << ipart << "\n";
		//cout <<"   nloc="<<this->particule[ipart].nloc<<"\n";
		sOK[ipart] = this->particule[ipart].dosimulpart(numscen);
		if (debuglevel == 5) cout << "apres dosimulpart de la particule " << ipart << "\n";

		if (sOK[ipart] == 0) {
			if (debuglevel == 5) cout << "\navant docalstat de la particule " << ipart << "   scenario " << this->particule[ipart].scen.number << "\n";
			for (gr = 1; gr < this->particule[ipart].ngr; gr++) this->particule[ipart].docalstat(gr);
			//this->particule[ipart].libere();
		}
		if (debuglevel == 5) cout << "apres docalstat de la particule " << ipart << "\n";
	}
	///fin du pragma
	if (debuglevel == 5) cout << "apres pragma\n";
	//cout<<"acceptation ="<<(double)(100*naccept)/(double)ntentes<<"   ("<<naccept<<" sur "<<ntentes<<")\n";
	//exit(1);
	int nph, npm, nn = 0;
	splitwords(header.entetehist, " ", ss);
	nph = ss.size();
	for (int m = 0; m < nph; m++) {
		ccc = ss[m][0];
		if (not isalnum(ccc)) nn++;
	}
	nph -= nn;
	//nph--;
	if (debuglevel == 5) cout << header.entetehist << "\n";
	if (header.entetemut.length() > 10) {
		splitwords(header.entetemut, " ", ss);
		npm = ss.size();
	}
	else npm = 0;
	//cout<<"nph="<<nph<<"   npm="<<npm<<"\n";
	for (int ipart = 0; ipart < this->npart; ipart++) {
		enreg[ipart].numscen = 1;
		if (this->particule[ipart].nscenarios > 1) {
			enreg[ipart].numscen = this->particule[ipart].scen.number;
		}
		if (sOK[ipart] == 0) {
			if (debuglevel == 5) cout << "dans particleset ipart=" << ipart << "     nparamvar=" << this->particule[ipart].scen.nparamvar << "\n";
			for (int j = 0; j < this->particule[ipart].scen.nparamvar; j++) {
				//for (int j=0;j<nph+npm;j++) {	
				if (debuglevel == 5) cout << "this->particule[" << ipart << "].scen.paramvar[" << j << "] = " << this->particule[ipart].scen.paramvar[j] << "     ";
				enreg[ipart].param[j] = this->particule[ipart].scen.paramvar[j];
				//cout<<"enreg["<<ipart<<"].param["<<j<<"] = "<<enreg[ipart].param[j]<<"\n";
				if (debuglevel == 5) cout << this->particule[ipart].scen.paramvar[j] << "  (" << enreg[ipart].param[j] << ")\n";
			}
			if (debuglevel == 5) cout << "\n";
			nstat = 0;
			for (int gr = 1; gr < this->particule[ipart].ngr; gr++) {
				/*cout<<"\ngroupe "<<gr<<"\n";
				for (int st=0;st<this->particule[ipart].grouplist[gr].nstat;st++) cout<<this->particule[ipart].grouplist[gr].sumstat[st].val<<"   \n";
				cout<<"\n";*/
				for (int st = 0; st < this->particule[ipart].grouplist[gr].nstat; st++) {
					/////////////// filtrage des Aml et aml début
					if (this->particule[ipart].grouplist[gr].sumstat[st].val == -9999.0) {
						if ((this->particule[ipart].grouplist[gr].sumstat[st].cat == 12)or (this->particule[ipart].grouplist[gr].sumstat[st].cat == -14))
							this->particule[ipart].grouplist[gr].sumstat[st].val = 0.5;
					}
					/////////////// filtrage des Aml et aml fin
					enreg[ipart].stat[nstat] = this->particule[ipart].grouplist[gr].sumstat[st].val;
					nstat++;/*cout<<nstat<<"\n";*/
				}
			}
			if (debuglevel == 5) {
				for (int st = 0; st < nstat; st++) cout << enreg[ipart].stat[st] << "   ";
				cout << "\n\n";
			}
			enreg[ipart].message = "OK";
		}
		else {
			if (this->particule[ipart].drawuntil) {
				//cout<<"A gene genealogy failed in scenario "<<this->particule[ipart].scen.number<<"\n";
				enreg[ipart].message = "A gene genealogy failed in scenario ";
				enreg[ipart].message += IntToString(this->particule[ipart].scen.number);
				enreg[ipart].message += ". Check consistency of this scenario over possible prior ranges of time event parameters.";
				enreg[ipart].message = "sOK[ipart] = " + IntToString(sOK[ipart]);
			}
			else enreg[ipart].message = "NOTOK";
		}
	}
	if (debuglevel == 5) cout << "apres remplissage des enreg\n";
	fflush(stdin);
	if (firsttime) {
		if (debuglevel == 5) cout << "FIRSTTIME\n";
		FILE* pFile;
		int categ, iq, icur = 0;
		bool trouve2;
		if (depuis == 0) {
			scurfile = header.pathbase + "first_records_of_the_reference_table_" + IntToString(icur) + ".txt";
			while (ifstream(scurfile.c_str())) {
				icur++;
				//scurfile=this->header.pathbase+"courant_"+IntToString(icur)+".log";
				scurfile = header.pathbase + "first_records_of_the_reference_table_" + IntToString(icur) + ".txt";
			}
		}
		if (debuglevel == 5) cout << "ecriture du fichier courant = " << scurfile << "\n";
		pFile = fopen(scurfile.c_str(), "w");
		fprintf(pFile, "%s\n", header.entete.c_str());
		if (debuglevel == 5) {
			cout << "FISTTIME\n";
			cout << header.entete << "\n\n";
			cout << header.entetehist << "\n\n";
			cout << header.entetemut << "\n\n";
			cout << header.entetestat << "\n\n";
		}
		/*int nph,npm;
		ss=splitwords(header.entetehist," ",&nph);
		if (debuglevel==5) cout<<header.entetehist<<"\n";
		if (header.entetemut.length()>10) ss=splitwords(header.entetemut," ",&npm); else npm=0;*/
		splitwords(header.entete, " ", ss);
		ns = ss.size();
		if (debuglevel == 5) cout << "nph=" << nph << "   npm=" << npm << "   ns=" << ns << "\n";
		np = ns - header.nstat - 1;
		//cout<<"ns="<<ns<<"  nparam="<<np<<"   nparamut="<<rt.nparamut<<"   nstat="<<header.nstat<<"\n";
		//cout<<"nph="<<nph<<"    npm="<<npm<<"\n";
		//cout<<"npart="<<this->npart<<"\n";
		for (int ipart = 0; ipart < this->npart; ipart++) {
			//cout<<"ipart="<<ipart<<"    sOK[ipart]="<<sOK[ipart]<<"   message="<<enreg[ipart].message<<"\n";
			if (sOK[ipart] == 0) {
				//cout<<enreg[ipart].numscen<<"\n";
				fprintf(pFile, "%3d  ", enreg[ipart].numscen);
				iscen = enreg[ipart].numscen - 1;
				//cout<<"scenario "<<enreg[ipart].numscen<<"\n";
				if (debuglevel == 5) cout << header.nparamtot << "   " << rt.nhistparam[iscen] << "   " << np << "\n";
				pa = 0;
				for (int j = 1; j < nph; j++) {
					trouve = false;
					ip = -1;
					while ((not trouve)and (ip < rt.nhistparam[iscen] - 1)) {
						ip++;
						trouve = (ss[j] == rt.histparam[iscen][ip].name);
						//cout<<"->"<<ss[j]<<"<-   ?   ->"<< header.scenario[iscen].histparam[ip].name<<"<-"<<trouve<<"\n";
					}
					if (trouve) {
						trouve2 = false;
						iq = -1;
						while ((not trouve2)and (iq < header.nparamtot - 1)) {
							iq++;
							trouve2 = (ss[j] == header.histparam[iq].name);
						}
						if (trouve2) categ = header.histparam[iq].category;
						else categ = 0;
						if (categ < 2) fprintf(pFile, "  %12.0f", enreg[ipart].param[ip]);
						else fprintf(pFile, "  %12.3f", enreg[ipart].param[ip]);
						pa++;
					}
					else fprintf(pFile, "              ");
				}
				if (debuglevel == 5) cout << "pa=" << pa << "   rt.nparam[iscen]=" << rt.nparam[iscen] << "\n";
				for (int j = 0; j < npm; j++) fprintf(pFile, "  %12.3e", enreg[ipart].param[pa + j]);
				for (int st = 0; st < header.nstat; st++) fprintf(pFile, "  %12.6f", enreg[ipart].stat[st]);
				fprintf(pFile, "\n");
				if (debuglevel == 5) {
					cout << "parametres de la particule" << ipart << "\n";
					for (int j = 0; j < npm + pa; j++) cout << enreg[ipart].param[j] << "   ";
					cout << "\n";
					cout << "\nstat de la particule " << ipart << "\n";
					for (int st = 0; st < header.nstat; st++) cout << enreg[ipart].stat[st] << "   ";
					cout << "\n";
				}
			}
		}
		fclose(pFile);
		if (debuglevel == 5) cout << "apres fclose\n";
	}
	else {
		if (depuis != 0) {
			FILE* pFile;
			int categ, iq;
			bool trouve2;
			pFile = fopen(scurfile.c_str(), "a");
			splitwords(header.entete, " ", ss);
			ns = ss.size();
			if (debuglevel == 5) cout << "nph=" << nph << "   npm=" << npm << "   ns=" << ns << "\n";
			np = ns - header.nstat - 1;
			for (int ipart = 0; ipart < this->npart; ipart++) {
				if (sOK[ipart] == 0) {
					//cout<<enreg[ipart].numscen<<"\n";
					fprintf(pFile, "%3d  ", enreg[ipart].numscen);
					iscen = enreg[ipart].numscen - 1;
					//cout<<"scenario "<<enreg[ipart].numscen<<"\n";
					if (debuglevel == 5) cout << header.nparamtot << "   " << rt.nhistparam[iscen] << "   " << np << "\n";
					pa = 0;
					for (int j = 1; j < nph; j++) {
						trouve = false;
						ip = -1;
						while ((not trouve)and (ip < rt.nhistparam[iscen] - 1)) {
							ip++;
							trouve = (ss[j] == rt.histparam[iscen][ip].name);
							//cout<<"->"<<ss[j]<<"<-   ?   ->"<< header.scenario[iscen].histparam[ip].name<<"<-"<<trouve<<"\n";
						}
						if (trouve) {
							trouve2 = false;
							iq = -1;
							while ((not trouve2)and (iq < header.nparamtot - 1)) {
								iq++;
								trouve2 = (ss[j] == header.histparam[iq].name);
							}
							if (trouve2) categ = header.histparam[iq].category;
							else categ = 0;
							if (categ < 2) fprintf(pFile, "  %12.0f", enreg[ipart].param[ip]);
							else fprintf(pFile, "  %12.3f", enreg[ipart].param[ip]);
							pa++;
						}
						else fprintf(pFile, "              ");
					}
					if (debuglevel == 5) cout << "pa=" << pa << "   rt.nparam[iscen]=" << rt.nparam[iscen] << "\n";
					for (int j = 0; j < npm; j++) fprintf(pFile, "  %12.3e", enreg[ipart].param[pa + j]);
					for (int st = 0; st < header.nstat; st++) fprintf(pFile, "  %12.6f", enreg[ipart].stat[st]);
					fprintf(pFile, "\n");
					if (debuglevel == 5) {
						cout << "parametres de la particule" << ipart << "\n";
						for (int j = 0; j < npm + pa; j++) cout << enreg[ipart].param[j] << "   ";
						cout << "\n";
						cout << "\nstat de la particule " << ipart << "\n";
						for (int st = 0; st < header.nstat; st++) cout << enreg[ipart].stat[st] << "   ";
						cout << "\n";
					}
				}
			}
			fclose(pFile);
		}
	}
	//cout<<"apres remise en ordre des enreg.param\n";
	if (debuglevel == 5) cout << "avant delete sOK\n";
	delete [] sOK;
	sOK = NULL;
	if (debuglevel == 5) cout << "apres delete sOK\n";
	//exit(1);
}

/**
 * Structure ParticleSet : simulation des particules utilisées pour la création de fichiers genepop
 */
string* ParticleSetC::simulgenepop(int npart, bool multithread, int seed) {
	bool dnatrue = true, simulfile = true;
	int numscen = 1;
	this->npart = npart;
	int sOK;
	string* ss;
	this->particule = vector<ParticleC>(this->npart);
	ss = new string[npart];
	//cout<<"avant le remplissage des "<<this->npart<<" particules\n";
	dataobs = datasim;
	init_particule(npart, dnatrue, simulfile, seed);
	/*for (int p=0;p<this->npart;p++) {
		this->particule[p].dnatrue = dnatrue;
		//this->setdata(p);//cout<<"apres setdata\n";
		this->setgroup(p);//cout<<"apres setgroup\n";
		this->setloci(p);//cout<<"apres setloci\n";
		this->setscenarios(p,true);//cout<<"apres setscenario\n";
		this->particule[p].mw.randinit(p,seed);
	}*/
	//cout<<"avant omp\n";
	//cout<<"Dans Particleset avant le dosimulpart snimoy="<<this->particule[0].grouplist[1].snimoy<<"\n";
#pragma omp parallel for shared(sOK) if(multithread)
	for (int ipart = 0; ipart < this->npart; ipart++) {
		sOK = this->particule[ipart].dosimulpart(numscen);
		//cout<<"sOK["<<ipart<<"]="<<sOK[ipart]<<"\n";
		if (sOK == 0) ss[ipart] = this->particule[ipart].dogenepop();
		else ss[ipart] = "";
		//cout<<ss[ipart]<<"\n";
	}
	return ss;
}

/**
 * Structure ParticleSet : simulation des particules utilisées pour la création de fichiers SNP
 */
string* ParticleSetC::simuldataSNP(int npart, bool multithread, int seed) {
	int numscen = 1;
	this->npart = npart;
	int* sOK;
	string* ss;
	sOK = new int[npart];
	this->particule = vector<ParticleC>(this->npart);
	ss = new string[npart];
	cout << "avant le remplissage des " << this->npart << " particules\n";
	dataobs = datasim;
	//cout<<"apres dataobs=datasim\n";
	for (int p = 0; p < this->npart; p++) {
		this->particule[p].maf = datasim.maf;
		//this->setdata(p);//cout<<"apres setdata\n";
		//this->setgroup(p);cout<<"apres setgroup\n"; // Je supprime cette ligne car pas de groupe de locus dans les parties simul de jeux de donnees. Pierre
		this->particule[p].ngr = groupe.size();
		this->particule[p].grouplist = vector<LocusGroupC>(this->particule[p].ngr + 1);
		for (int gr = 0; gr < this->particule[p].ngr + 1; gr++) {
			this->particule[p].grouplist[gr].type = 2;
			this->particule[p].grouplist[gr].nloc = groupe[gr].nloc;

		} //cout<<"avant setloci\n";
		this->setloci(p);//cout<<"apres setloci\n";
		this->setscenarios(p, true);//cout<<"apres setscenario\n";
		this->particule[p].mw.randinit(p, seed);
	}
	//cout<<"avant omp\n";
#pragma omp parallel for shared(sOK) if(multithread)
	for (int ipart = 0; ipart < this->npart; ipart++) {
		//cout<<"avant dosimulpart\n";
		sOK[ipart] = this->particule[ipart].dosimulpart(numscen);
		//cout<<"sOK["<<ipart<<"]="<<sOK[ipart]<<"\n";
		if (sOK[ipart] == 0) ss[ipart] = this->particule[ipart].dodataSNP();
		else ss[ipart] = "";
		//cout<<ss[ipart]<<"\n";
	}
	return ss;
}


/**
 * Structure ParticleSet : simulation des particules utilisées pour la création de fichiers SNP poolseq
 */
string* ParticleSetC::simuldataSNPool(int npart, bool multithread, int seed) {
	int numscen = 1;
	this->npart = npart;
	int* sOK;
	string* ss;
	sOK = new int[npart];
	this->particule = vector<ParticleC>(this->npart);
	ss = new string[npart];
	cout << "SNP Poolseq avant le remplissage des " << this->npart << " particules\n";
	dataobs = datasim;
	cout << "apres dataobs=datasim\n";
	for (int p = 0; p < this->npart; p++) {
		cout << "debut de la particule " << p + 1 << "\n";
		this->particule[p].mrc = datasim.mrc;
		this->particule[p].lambdapoolseq = datasim.lambdapoolseq;
		cout << "apres mrc\n";
		//this->setdata(p);//cout<<"apres setdata\n";
		//this->setgroup(p);cout<<"apres setgroup\n"; // Je supprime cette ligne car pas de groupe de locus dans les parties simul de jeux de donnees. Pierre
		this->particule[p].ngr = groupe.size();
		this->particule[p].grouplist = vector<LocusGroupC>(this->particule[p].ngr);
		cout << "avant la boucle sur les groupes\n";
		for (int gr = 0; gr < this->particule[p].ngr; gr++) {
			this->particule[p].grouplist[gr].type = 2;
			cout << "apres type=2\n";
			cout << "groupe[" << gr << "].nloc=" << groupe[gr].nloc << "\n";
			this->particule[p].grouplist[gr].nloc = groupe[gr].nloc;
			cout << "après nloc=" << groupe[gr].nloc << "\n";
		}
		cout << "avant setloci\n";
		this->setloci(p);
		cout << "apres setloci\n";
		this->setscenarios(p, true);
		cout << "apres setscenario\n";
		this->particule[p].mw.randinit(p, seed);
	}
	//cout<<"avant omp\n";
#pragma omp parallel for shared(sOK) if(multithread)
	for (int ipart = 0; ipart < this->npart; ipart++) {
		//cout<<"avant dosimulpart\n";
		sOK[ipart] = this->particule[ipart].dosimulpart(numscen);
		//cout<<"sOK["<<ipart<<"]="<<sOK[ipart]<<"\n";
		if (sOK[ipart] == 0) ss[ipart] = this->particule[ipart].dodataSNPool();
		else ss[ipart] = "";
		//cout<<ss[ipart]<<"\n";
	}
	return ss;
}
