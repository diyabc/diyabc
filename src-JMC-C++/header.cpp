/*



*/
#include <cstdio>
#include <cmath>
#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>
//#include <complex>
#include <iomanip>
#include <utility>
#include <string>
#include <cstring>
#include <algorithm>
#include <cstdlib>
#include <stdexcept>
#include <time.h>

#include "randomgenerator.h"
#include "mesutils.h"
#include "data.h"
#include "history.h"
#include "particuleC.h"
#include "history.h"


extern int debuglevel;
/*
#ifndef PARTICULEC
#include "particuleC.cpp"
#define PARTICULEC
#endif

#ifndef DATA
#include "data.cpp"
#define DATA
#endif
*/

#include "header.h"

extern string * stat_type;
extern int* stat_num;
extern string path;
extern bool randomforest,datafile;
extern DataC dataobs,datasim;
extern vector <LocusGroupC> groupe;
extern vector <ScenarioC> scenario;

ParticleC particuleobs;

/**
 * Copie du contenu d'une classe MutParameterC
 */
/*  MutParameterC::MutParameterC (MutParameterC const & source) {
	this->name = source.name;
	this->groupe = source.groupe;
	this->category = source.category;
	this->value = source.value;
	this->prior = source.prior;
  }*/
  
/**
 * Definition de l'operateur = pour une instance de la classe MutParameterC
 */
  MutParameterC & MutParameterC::operator= (MutParameterC const & source) {
	if (this== &source) return *this;
	this->name = source.name;
	this->groupe = source.groupe;
	this->category = source.category;
	this->value = source.value;
	this->prior = source.prior;
	return *this;
  }
  
  
/**
 * Copie du contenu d'une classe HeaderC
 */
 /* HeaderC::HeaderC (HeaderC const & source) {
    this->message = source.message;
    this->datafilename = source.datafilename;
    this->entete = source.entete;
    this->entetehist = source.entetehist;
    this->entetemut = source.entetemut;
    this->entetemut0 = source.entetemut0;
    this->entetestat = source.entetestat;
    this->pathbase = source.pathbase;
    this->nparamtot = source.nparamtot;
    this->nstat = source.nstat;
    this->nscenarios = source.nscenarios;
    this->nconditions = source.nconditions;
    this->ngroupes = source.ngroupes;
    this->nparamut = source.nparamut;
    this->nsimfile = source.nsimfile;
    this->drawuntil = source.drawuntil;
    this->reference = source.reference;
    this->message = source.message;
	this->threshold = source.threshold;
	this->reffreqmin = source.reffreqmin;
    this->drawuntil = source.drawuntil;
	this->reference = source.reference;
//
	this->scen = source.scen;
	this->dataobs = source.dataobs;
///// Tableaux à dimmensionner	
	this->stat_obs = new float[this->nstat];
	for (int i=0;i<this->nstat;i++) this->stat_obs[i]=source.stat_obs[i];
	this->mutparam = new MutParameterC[this->nparamut];
	for (int i=0;i<this->nparamut;i++) this->mutparam[i]=source.mutparam[i];
	this->statname = new string[this->nstat];
	for (int i=0;i<this->nstat;i++) this->statname[i]=source.statname[i];
	this->histparam = new HistParameterC[this->nparamtot];
	for (int i=0;i<this->nparamtot;i++) this->histparam[i] = source.histparam[i];
	this->condition = new ConditionC[this->nconditions];
	for (int i=0;i<this->nconditions;i++) this->condition[i] = source.condition[i];
	scenario = new ScenarioC[this->nscenarios];
	for (int i=0;i<this->nscenarios;i++) scenario[i] = source.scenario[i];
	this->groupe = new LocusGroupC[this->ngroupes];
	for (int i=0;i<this->ngroupes;i++) this->groupe[i] = source.groupe[i];
}
*/
/**
 * Definition de l'operateur = pour une instance de la classe HeaderC
 */
  /*HeaderC & HeaderC::operator= (HeaderC const & source) {
	if (this == &source)  return *this;
// Suppression des tableaux non nuls
	if (this->stat_obs != NULL) delete [] this->stat_obs;
	if (this->mutparam != NULL) delete [] this->mutparam;
	if (this->statname != NULL) delete [] this->statname;
	if (this->histparam != NULL) delete [] this->histparam;
	if (this->condition != NULL) delete [] this->condition;
	if (scenario != NULL) delete [] scenario;
	if (this->groupe != NULL) delete [] this->groupe;
// recopie des constituants	
    this->message = source.message;
    this->datafilename = source.datafilename;
    this->entete = source.entete;
    this->entetehist = source.entetehist;
    this->entetemut = source.entetemut;
    this->entetemut0 = source.entetemut0;
    this->entetestat = source.entetestat;
    this->pathbase = source.pathbase;
    this->nparamtot = source.nparamtot;
    this->nstat = source.nstat;
    this->nstatsnp = source.nstatsnp;
	this->nscenarios = source.nscenarios;
    this->nconditions = source.nconditions;
    this->ngroupes = source.ngroupes;
    this->nparamut = source.nparamut;
    this->nsimfile = source.nsimfile;
    this->drawuntil = source.drawuntil;
    this->reference = source.reference;
    this->message = source.message;
	this->threshold = source.threshold;
	this->reffreqmin = source.reffreqmin;
    this->drawuntil = source.drawuntil;
	this->reference = source.reference;

	this->scen = source.scen;
	this->dataobs = source.dataobs;
	if (source.stat_obs != NULL){
		this->stat_obs = new float[this->nstat];
		for (int i=0;i<this->nstat;i++) this->stat_obs[i]=source.stat_obs[i];
	}
	if (source.stat_obs != NULL){
		this->mutparam = new MutParameterC[this->nparamut];
		for (int i=0;i<this->nparamut;i++) this->mutparam[i]=source.mutparam[i];
	}
	if (source.statname != NULL){
		this->statname = new string[this->nstat];
		for (int i=0;i<this->nstat;i++) this->statname[i]=source.statname[i];
	}
	if (source.histparam != NULL){
		this->histparam = new HistParameterC[this->nparamtot];
		for (int i=0;i<this->nparamtot;i++) this->histparam[i] = source.histparam[i];
	}
	if (source.condition != NULL){
		this->condition = new ConditionC[this->nconditions];
		for (int i=0;i<this->nconditions;i++) this->condition[i] = source.condition[i];
	}
	if (source.scenario != NULL){
		scenario = new ScenarioC[this->nscenarios];
		for (int i=0;i<this->nscenarios;i++) scenario[i] = source.scenario[i];
	}
	if (source.groupe != NULL){
		this->groupe = new LocusGroupC[this->ngroupes+1];
		for (int i=0;i<this->ngroupes+1;i++) this->groupe[i] = source.groupe[i];
	}
	return *this;
  }	
	*/
	
	
void MutParameterC::ecris() {
	cout<<"    groupe="<<this->groupe<<"   category="<<this->category<<"\n";
	//prior.ecris();
}

void HeaderC::assignloc(int gr){
	groupe[gr].nloc = 0;
	for (int loc=0;loc<dataobs.nloc;loc++) {
		if (dataobs.locus[loc].groupe==gr) groupe[gr].nloc++;
	}
	if (debuglevel==2) cout<<"assignloc nloc="<<groupe[gr].nloc<<"\n";
	groupe[gr].loc = vector<int>(groupe[gr].nloc);
	int iloc=-1;
	for (int i=0;i<dataobs.nloc;i++) {
		if (dataobs.locus[i].groupe==gr) {iloc++;groupe[gr].loc[iloc] = i;}
	}
}

/*ScenarioC superscen() {
                int neventm,nn0m,nsampm,nparamm,nparamvarm,nconditionm;
                neventm=0;    for (int i=0;i<this->nscenarios;i++) {if (neventm<scenario[i].nevent)          neventm=scenario[i].nevent;}
                nn0m=0;       for (int i=0;i<this->nscenarios;i++) {if (nn0m<scenario[i].nn0)                nn0m=scenario[i].nn0;}
                nsampm=0;     for (int i=0;i<this->nscenarios;i++) {if (nsampm<scenario[i].nsamp)            nsampm=scenario[i].nsamp;}
                nparamm=0;    for (int i=0;i<this->nscenarios;i++) {if (nparamm<scenario[i].nparam)         nparamm=scenario[i].nparam;}
                nparamvarm=0; for (int i=0;i<this->nscenarios;i++) {if (nparamvarm<scenario[i].nparamvar)    nparamvarm=scenario[i].nparamvar;}
                nconditionm=0;for (int i=0;i<this->nscenarios;i++) {if (nconditionm<scenario[i].nconditions) nconditionm=scenario[i].nconditions;}
                cout<<"neventm="<<neventm<<"   nn0m="<<nn0m<<"   nsampm="<<nsampm<<"   nparamm="<<nparamm<<"nparamvarm="<<nparamvarm<<"   ncondm="<<nconditionm<<"\n";
                this->scen.event = new EventC[neventm];
                this->scen.ne0   = new Ne0C[nn0m];

                for (int i=0;i<nn0m;i++) this->scen.ne0[i].name="blabla";
                this->scen.time_sample = new int[nsampm];
                this->scen.histparam = new HistParameterC[nparamm];
                this->scen.paramvar = new double[nparamvarm+3];
                if (nconditionm>0) this->scen.condition = new ConditionC[nconditionm];

        }*/

int HeaderC::readHeaderDebut(ifstream & file){
	string s1;
	int nl=0;
	getline(file,this->datafilename);nl++;
	this->pathbase=path;
	//cout <<"avant loadfromfile du fichier "<<this->datafilename<<"\n";
	int error = dataobs.loadfromfile(path+this->datafilename);
	//cout<<"apres loadfromfile  error="<<error<<"\n";
	if (error != 0) {
		this->message = dataobs.message;
		throw std::runtime_error(message);
		return error;
	}
	if (debuglevel==1) cout<<"lecture du fichier data terminée\n";
	getline(file,s1);nl++;
	this->nparamtot=getwordint(s1,1);		cout<<"nparamtot="<<this->nparamtot<<"\n";
	this->nstat=getwordint(s1,4);			cout<<"nstat="<<this->nstat<<"\n";
	//cout<<"avant scenarios\n";fflush(stdin);
	return 0;
}

int HeaderC::readHeaderScenarios(ifstream & file){
		//Partie Scenarios
	string s1;
	int nl = 0;
	getline(file,s1);nl++;   //cout<<s1<<"\n";     //ligne vide
	getline(file,s1);nl++;    //cout<<s1<<"\n";    // nombre de scenarios
	this->nscenarios=getwordint(s1,1); cout<<this->nscenarios<<" scenario(s)\n";
	string** sl; sl = new string*[this->nscenarios];
	int * nlscen = new int[this->nscenarios];

	scenario = vector<ScenarioC>(this->nscenarios);
	for (int i=0;i<this->nscenarios;i++) {nlscen[i]=getwordint(s1,3+i);}
	for (int i=0;i<this->nscenarios;i++) {
		sl[i] = new string[nlscen[i]];
		getline(file,s1); nl++;   //cout<<s1<<"\n";
		if (s1.find("scenario")!=0) {
			this->message = "Error when reading header.txt file :keyword <scenario> expected at line "+IntToString(nl)+". Check the number of lines of each scenario at line 4.";
			throw std::runtime_error(message);
			return 1;
		}
		scenario[i].number = getwordint(s1,2);
		scenario[i].prior_proba = getwordfloat(s1,3);
		scenario[i].nparam = 0;
		scenario[i].nparamvar=0;
		for (int j=0;j<nlscen[i];j++) {getline(file,sl[i][j]);nl++;/*cout<<sl[i][j]<<"\n";*/}
		this->message = scenario[i].read_events(nlscen[i],sl[i]);
		if (this->message!="") {
			this->message = "Error when reading  header.txt file : "+this->message+" in scenario "+IntToString(i+1);
			return 1;
		}
		//scenario[i].ecris();
		//cout<<"apres read_events\n";
		this->message = scenario[i].checklogic();
		if (this->message!="") {
			this->message = "Error when reading  header.txt file : "+this->message+" in scenario "+IntToString(i+1);
			return 1;
		}
		//cout<<"apres checklogic\n";
	}
	this->reference=false;

	for (int i=0;i<this->nscenarios;i++) delete []sl[i];
	delete [] sl;
	delete [] nlscen;
	if (debuglevel==2) cout<<"header.txt : fin de la lecture de la partie scénarios\n";
	return 0;
}

int HeaderC::readHeaderHistParam(ifstream & file){
	string s1, s2;
	vector<string> ss, ss2;
	int j, k;

	//cout <<"avant histparam\n";fflush(stdin);
	getline(file,s1);       //ligne vide
	getline(file,s1);
	//cout<<s1<<"\n";
	splitwords(s1," ",ss);
	s2=ss[3].substr(1,ss[3].length()-2);
	splitwords(s2,",",ss2);
	this->nparamtot = atoi(ss2[0].c_str());
	this->nconditions = atoi(ss2[1].c_str());
	//cout<<"dans header nconditions="<<this->nconditions<<"\n";
	//cout<<"avant  new HistParameterC[this->nparamtot]\n";

	this->histparam = vector<HistParameterC>(this->nparamtot);
	this->condition.clear();
	for (int i=0;i<this->nparamtot;i++) {
		getline(file,s1);
		//cout<<s1<<"\n";
		splitwords(s1," ",ss); 
		this->histparam[i].name=ss[0];
		//cout<<this->histparam[i].name<<"\n";
		if (ss[1]=="N") this->histparam[i].category = 0;
		else if  (ss[1]=="T") this->histparam[i].category = 1;
		else if  (ss[1]=="A") this->histparam[i].category = 2;
		this->histparam[i].prior.readprior(ss[2]);
		//    cout<<"readHeader "<<this->histparam[i].name;
		//    if (this->histparam[i].prior.constant) cout<<"   constant\n"; else cout<<"   variable\n";
	}
	//cout<<"apres les readprior\n";
	this->drawuntil=true;
	if (this->nconditions>0) {
		this->condition = vector<ConditionC>(this->nconditions);
		for (int i=0;i<this->nconditions;i++) {
			getline(file,s1);
			//cout<<s1<<"\n";
			this->condition[i].readcondition(s1);
		}
		getline(file,s1);
		if (s1 != "DRAW UNTIL") this->drawuntil=false; else  this->drawuntil=true;
	}
	if (debuglevel==2) cout<<"header.txt : fin de la lecture de la partie priors des paramètres démographiques\n";
	//retour sur les parametres spécifiques à chaque scenario;
	//cout <<"avant retour sur histparam/scenario\n";fflush(stdin);
	for (int i=0;i<this->nscenarios;i++) {
		scenario[i].nparamvar=0;
		for (int j=0;j<scenario[i].nparam;j++) {
			int k=0;
			while (scenario[i].histparam[j].name!=this->histparam[k].name) {k++;}
			//scenario[i].histparam[j].prior = copyprior(this->histparam[k].prior);
			scenario[i].histparam[j].prior = this->histparam[k].prior;
			if (not this->histparam[k].prior.constant) {
				//scenario[i].paramvar[scenario[i].nparamvar].prior = copyprior(this->histparam[k].prior);
				scenario[i].nparamvar++;
			}
			//scenario[i].histparam[j].ecris(true);
		}
		//cout<<"scenario "<<i+1<<"   "<<scenario[i].nparam<<" param et "<<scenario[i].nparamvar<<" paramvar\n "<<flush;
	}
	//retour sur les conditions spécifiques à chaque scenario
	//cout <<"avant retour sur conditions\n";fflush(stdin);
	if (this->nconditions>0) {
		for (int i=0;i<this->nscenarios;i++) {
			int nc=0;
			for (j=0;j<this->nconditions;j++) {
				int np=0;
				for (k=0;k<scenario[i].nparam;k++) {
					//cout<<this->condition[j].param1<<"   "<<this->condition[j].param2<<"   "<<scenario[i].histparam[k].name<<"  np="<<np<<"\n";
					if (this->condition[j].param1 == scenario[i].histparam[k].name) np++;
					if (this->condition[j].param2 == scenario[i].histparam[k].name) np++;
				}
				if (np==2) nc++;
			}
			scenario[i].nconditions=nc;
			//cout <<"header.scenario["<<i<<"].nconditions="<<scenario[i].nconditions<<"\n";
			scenario[i].condition = vector<ConditionC>(nc);
			nc=0;
			while (nc<scenario[i].nconditions) {
				for (j=0;j<this->nconditions;j++) {
					int np=0;
					for (k=0;k<scenario[i].nparam;k++) {
						if (this->condition[j].param1 == scenario[i].histparam[k].name) np++;
						if (this->condition[j].param2 == scenario[i].histparam[k].name) np++;
					}
					if (np==2) {
						scenario[i].condition[nc] = this->condition[j];
						//scenario[i].condition[nc].ecris();
						nc++;
					}
				}
			}
			//cout <<"dans readheader  \n";

		}
	} else for (int i=0;i<this->nscenarios;i++) scenario[i].nconditions=0;
	return 0;
}

int HeaderC::readHeaderLoci(ifstream & file){
	string s1;
	vector<string> ss;
	int k, k1, k2, nss, gr, grm, nsg, nl,typ;
		//Partie loci description
	//cout <<"avant partie loci\n";fflush(stdin);
	getline(file,s1);       //ligne vide
	getline(file,s1);      //cout<<"readHeaderLoci  s1="<<s1<<"\n"; //ligne "loci description"
	if (dataobs.filetype==0){ //fichier GENEPOP
		grm=1;                  //nombre maximal de groupes
		for (int loc=0;loc<dataobs.nloc;loc++){
			getline(file,s1);
			splitwords(s1," ",ss); nss = ss.size();
			k=0;while (ss[k].find("[")==string::npos) k++;
			if (ss[k]=="[M]") {
				s1=ss[k+1].substr(1,ss[k+1].length());gr=atoi(s1.c_str());dataobs.locus[loc].groupe=gr;if (gr>grm) grm=gr;
				dataobs.locus[loc].motif_size=atoi(ss[k+2].c_str());dataobs.locus[loc].motif_range=atoi(ss[k+3].c_str());
			}
			else if (ss[k]=="[S]") {
				s1=ss[k+1].substr(1,ss[k+1].length());gr=atoi(s1.c_str());dataobs.locus[loc].groupe=gr;if (gr>grm) grm=gr;
				dataobs.locus[loc].mutsit.resize(dataobs.locus[loc].dnalength);
				//cout<<dataobs.locus[loc].dnalength<<"\n";
				//dataobs.locus[loc].dnalength=atoi(ss[k+2].c_str());  //inutile variable déjà renseignée
			}
		}
	} else {			//fichier SNP
		//cout<<"fichier SNP dans redheaderLoci\n";
		this->ngroupes=getwordint(s1,3); cout<<s1<<"\n";cout<<"this->ngroupes="<<this->ngroupes<<"\n";
		groupe = vector<LocusGroupC>(this->ngroupes+1);
		groupe[0].nloc = dataobs.nloc;//cout<<"groupe[0].nloc="<<groupe[0].nloc<<"\n";
		groupe[0].type = 2;
		int prem;
		for (int loc=0;loc<dataobs.nloc;loc++) dataobs.locus[loc].groupe=0;
		for (gr=1;gr<=this->ngroupes;gr++) {
			groupe[gr].type = 2;
			getline(file,s1);
			splitwords(s1," ",ss); nss = ss.size();
			nsg=(nss-3)/2;
			//cout<<"nss = "<<nss<<"   nsg = "<<nsg<<"\n";
			groupe[gr].nloc=0;
			groupe[gr].type=2;
			prem = getwordint(s1,nss)-1;
			cout<<"prem = "<<prem<<"\n";
			for (k=0;k<nsg;k++) groupe[gr].nloc +=getwordint(s1,2*k+1); cout<<"groupe["<<gr<<"].nloc="<<groupe[gr].nloc<<"\n";
			groupe[0].nloc -= groupe[gr].nloc;
			if (debuglevel==2) for (int kk=0;kk<=this->ngroupes;kk++) cout<<"groupe["<<kk<<"].nloc = "<<groupe[kk].nloc<<"\n";
			groupe[gr].loc = vector<int>(groupe[gr].nloc);
			k1=0;
			for (k=0;k<nsg;k++) {
				nl=getwordint(s1,2*k+1);
				k2=0;
				if (ss[2*k+1].find("<A>")!=string::npos) typ=10;
				else if (ss[2*k+1].find("<H>")!=string::npos) typ=11;
				else if (ss[2*k+1].find("<X>")!=string::npos) typ=12;
				else if (ss[2*k+1].find("<Y>")!=string::npos) typ=13;
				else if (ss[2*k+1].find("<M>")!=string::npos) typ=14;
				for (int loc=prem;loc<dataobs.nloc;loc++){
					if (dataobs.locus[loc].type==typ) {
						groupe[gr].loc[k1]=loc;
						dataobs.locus[loc].groupe=gr;
						k1++;
						k2++;
					}
					if (k2 == nl) break;
					if ((loc+1==dataobs.nloc)and(k2<nl)) {
						this->message = "Not enough loci of type ";
						if (typ==10) this->message += "<A>";
						if (typ==11) this->message += "<H>";
						if (typ==12) this->message += "<X>";
						if (typ==13) this->message += "<Y>";
						if (typ==14) this->message += "<M>";
						this->message += " found in the datafile from position "+IntToString(prem)+". This is perhaps due to the elimination of "+IntToString(-groupe[0].nloc)+" monomorphic loci from the original dataset.";
						return 1;
					}
				}
			}
		}
		//cout<<"avant recopie des locus\n";
		vector <LocusC> pivlocus;
		
		for (int i=0;i<groupe[1].nloc;i++) pivlocus.push_back(dataobs.locus[groupe[1].loc[i]]);
		cout<<"dataobs.locus[0].type="<<dataobs.locus[groupe[1].loc[0]].type<<"\n";
		cout<<"pivlocus[0].type="<<pivlocus[0].type<<"\n";
		for (int loc=0;loc<groupe[1].nloc;loc++) dataobs.locus[loc] = pivlocus[loc];
		cout<<"dataobs.locus[0].type="<<dataobs.locus[0].type<<"\n";
		pivlocus.clear();
		for (int loc=dataobs.nloc-1;loc>=groupe[1].nloc;loc--) dataobs.locus.pop_back();
		cout<<"dataobs.locus.size()="<<dataobs.locus.size()<<"\n";		
		dataobs.nloc = (int)dataobs.locus.size();
		for (int loc=0;loc<dataobs.nloc;loc++) groupe[1].loc[loc]=loc;
		

		
		/*
		//cout<<"groupe[0].nloc = "<<groupe[0].nloc<<"\n";
		if (groupe[0].nloc>0) {
			groupe[0].loc = vector<int>(groupe[0].nloc);
			k1=0;
			for (int loc=0;loc<dataobs.nloc;loc++) {
				//cout<<"loc="<<loc<<"\n";
				if (dataobs.locus[loc].groupe==0) {groupe[0].loc[k1]=loc;k1++;}
			}
		}*/
		//for (int i=0;i<groupe[1].nloc;i++) cout<<"Locus numéro "<<i+1<<" = locus initial "<<groupe[1].loc[i]<<"   (type "<<dataobs.locus[groupe[1].loc[i]].type<<")\n";
	}
	
	if (debuglevel==2) cout<<"header.txt : fin de la lecture de la partie description des locus\n";
	return 0;
}



int HeaderC::readHeaderGroupPrior(ifstream & file){
	string s1;
	vector<string> ss, ss1;
	int gr, j;
	if (dataobs.filetype==0){
		//cout <<"avant partie group priors\n";fflush(stdin);
		getline(file,s1);       //ligne vide
		getline(file,s1);     // cout<<"readHeaderGroupPrior s1="<<s1<<"\n"; //ligne "group prior"
		this->ngroupes=getwordint(s1,3);
		//cout<<"header.ngroupes="<<this->ngroupes;
		groupe = vector<LocusGroupC>(this->ngroupes+1);
		this->assignloc(0);
		//cout<<"on attaque les groupes : analyse des priors nombre de groupes="<<this->ngroupes <<"\n";
		for (gr=1;gr<=this->ngroupes;gr++){
			getline(file,s1);
			//cout<<s1<<"\n";
			splitwords(s1," ",ss);
			this->assignloc(gr);
			if (ss[2]=="[M]") {
				groupe[gr].type=0;
				getline(file,s1); splitwords(s1," ",ss1);
				groupe[gr].priormutmoy.readprior(ss1[1]);
				//cout<<"readHeaderGroupPrior priormutmoy.constant="<<groupe[gr].priormutmoy.constant<<"\n";
				groupe[gr].priormutmoy.fixed=false;
				if (groupe[gr].priormutmoy.constant) groupe[gr].mutmoy=groupe[gr].priormutmoy.mini;
				else {groupe[gr].mutmoy=-1.0;for (int i=0;i<this->nscenarios;i++) {scenario[i].nparamvar++;}}
				getline(file,s1);splitwords(s1," ",ss1);
				groupe[gr].priormutloc.readprior(ss1[1]);
				getline(file,s1);splitwords(s1," ",ss1);
				groupe[gr].priorPmoy.readprior(ss1[1]);
				groupe[gr].priorPmoy.fixed=false;
				if (groupe[gr].priorPmoy.constant) groupe[gr].Pmoy=groupe[gr].priorPmoy.mini;
				else {groupe[gr].Pmoy=-1.0;for (int i=0;i<this->nscenarios;i++) {scenario[i].nparamvar++;}}

				getline(file,s1);splitwords(s1," ",ss1);
				groupe[gr].priorPloc.readprior(ss1[1]);
				//cout<<"Ploc    ";groupe[gr].priorPloc.ecris();

				getline(file,s1);splitwords(s1," ",ss1);
				groupe[gr].priorsnimoy.readprior(ss1[1]);
				groupe[gr].priorsnimoy.fixed=false;
				if (groupe[gr].priorsnimoy.constant) groupe[gr].snimoy=groupe[gr].priorsnimoy.mini;
				else {groupe[gr].snimoy=-1.0;for (int i=0;i<this->nscenarios;i++) {scenario[i].nparamvar++;}}
				getline(file,s1);splitwords(s1," ",ss1);
				groupe[gr].priorsniloc.readprior(ss1[1]);
				//cout<<"sniloc  ";groupe[gr].priorsniloc.ecris();

			} else if (ss[2]=="[S]") {
				groupe[gr].type=1;
				getline(file,s1);splitwords(s1," ",ss1);
				groupe[gr].priormusmoy.readprior(ss1[1]);
				groupe[gr].priormusmoy.fixed=false;
				if (groupe[gr].priormusmoy.constant) groupe[gr].musmoy=groupe[gr].priormusmoy.mini;
				else {groupe[gr].musmoy=-1.0;for (int i=0;i<this->nscenarios;i++) {scenario[i].nparamvar++;}}

				getline(file,s1);splitwords(s1," ",ss1); 
				groupe[gr].priormusloc.readprior(ss1[1]);

				getline(file,s1);splitwords(s1," ",ss1);
				groupe[gr].priork1moy.readprior(ss1[1]);
				groupe[gr].priork1moy.fixed=false;
				if (groupe[gr].priork1moy.constant) groupe[gr].k1moy=groupe[gr].priork1moy.mini;
				else groupe[gr].k1moy=-1.0;

				getline(file,s1);splitwords(s1," ",ss1);
				groupe[gr].priork1loc.readprior(ss1[1]);

				getline(file,s1);splitwords(s1," ",ss1); 
				groupe[gr].priork2moy.readprior(ss1[1]);
				groupe[gr].priork2moy.fixed=false;
				if (groupe[gr].priork2moy.constant) groupe[gr].k2moy=groupe[gr].priork2moy.mini;
				else groupe[gr].k2moy=-1.0;

				getline(file,s1);splitwords(s1," ",ss1);
				groupe[gr].priork2loc.readprior(ss1[1]);

				getline(file,s1);splitwords(s1," ",ss1);
				groupe[gr].p_fixe=atof(ss1[2].c_str());groupe[gr].gams=atof(ss1[3].c_str());
				if (ss1[1]=="JK") groupe[gr].mutmod=0;
				else if (ss1[1]=="K2P") groupe[gr].mutmod=1;
				else if (ss1[1]=="HKY") groupe[gr].mutmod=2;
				else if (ss1[1]=="TN") groupe[gr].mutmod=3;
				//cout<<"mutmod = "<<groupe[gr].mutmod <<"\n";
				if (groupe[gr].mutmod>0) {
					if (not groupe[gr].priork1moy.constant)for (int i=0;i<this->nscenarios;i++) {scenario[i].nparamvar++;}
				}
				if (groupe[gr].mutmod==3) {
					if (not groupe[gr].priork2moy.constant)for (int i=0;i<this->nscenarios;i++) {scenario[i].nparamvar++;}
				}
			} // fin du if portant sur ss[2]
		}
		//cout<<"debuglevel = "<<debuglevel<<"\n";
		if (debuglevel==2) cout<<"header.txt : fin de la lecture de la partie définition des groupes\n";
		//Mise à jour des locus séquences
		MwcGen mwc;
		mwc.randinit(999,time(NULL));
		int nsv;
		bool nouveau;
		if (debuglevel==2) cout<<"avant la boucle sur les locus\n";
		for (int loc=0;loc<dataobs.nloc;loc++){
			gr=dataobs.locus[loc].groupe;
			if ((dataobs.locus[loc].type>4)and(dataobs.locus[loc].type<10)and(gr>0)){
				nsv = floor(dataobs.locus[loc].dnalength*(1.0-0.01*groupe[gr].p_fixe)+0.5);
				if (debuglevel==2) cout<<"nsv = "<<nsv<<"\n";
				if (debuglevel==2) cout<<"mutsit.size = "<<dataobs.locus[loc].mutsit.size()<<"\n";
				if (debuglevel==2) cout<<"groupe[gr].gams = "<<groupe[gr].gams<<"\n";
				if (debuglevel==2) cout<<"dataobs.locus[loc].dnalength = "<<dataobs.locus[loc].dnalength<<"\n";
				for (int i=0;i<dataobs.locus[loc].dnalength;i++) {dataobs.locus[loc].mutsit[i] = mwc.ggamma3(1.0,groupe[gr].gams);}
				if (debuglevel==2) cout<<"apres tirage dans gamma3\n";
				int *sitefix;
				sitefix=new int[dataobs.locus[loc].dnalength-nsv];
				for (int i=0;i<dataobs.locus[loc].dnalength-nsv;i++) {
					if (i==0) sitefix[i]=mwc.rand0(dataobs.locus[loc].dnalength);
					else {
						do {
							sitefix[i]=mwc.rand0(dataobs.locus[loc].dnalength);
							nouveau=true;j=0;
							while((nouveau)and(j<i)) {nouveau=(sitefix[i]!=sitefix[j]);j++;}
						} while (not nouveau);
					}
					dataobs.locus[loc].mutsit[i] = 0.0;
				}
				delete [] sitefix;
				double s=0.0;
				for (int i=0;i<dataobs.locus[loc].dnalength;i++) s += dataobs.locus[loc].mutsit[i];
				for (int i=0;i<dataobs.locus[loc].dnalength;i++) dataobs.locus[loc].mutsit[i] /=s;
			}
		}
		if (debuglevel==2) cout<<"header.txt : fin de la mise à jour des locus séquences\n";
		//cout<<"avant la mise à jour des paramvar\n";fflush(stdin);
	} //fin de la condition fichier=GENEPOP
	else 
	//Mise à jour des paramvar
	for (int i=0;i<this->nscenarios;i++) {
		scenario[i].paramvar = vector<double>(scenario[i].nparamvar);
		for (int j=0;j<scenario[i].nparamvar;j++)scenario[i].paramvar[j]=-1.0;
		//scenario[i].ecris();
	}

	return 0;
}

int HeaderC::readHeaderAllStat(ifstream & file) {
	string s1;
	int j, k, gr,nstatgr,nsamp=dataobs.nsample;
	//cout<<"debut de readHeaderAllStat\n";
	//cout<<"nsamp="<<nsamp<<"    "<<"nstat="<<this->nstat<<"\n";
    //cout<<"Lecture bidon du header\n";
	getline(file,s1);//cout<<s1<<"\n";       //ligne vide
	do {getline(file,s1);/*cout<<s1<<"\n";*/} while (s1!="");
	getline(file,this->entete);     //ligne entete
	this->entetehist=this->entete.substr(0,this->entete.length()-14*(nparamut+nstat));
	if (nparamut>0) this->entetemut=this->entete.substr(this->entetehist.length(),14*nparamut);else this->entetemut="";
	//cout<<this->entete<<"\n";
	
	this->nstat=0;
    //cout<<"Fin de la lecture bidon du header\n";
	vector <StatsnpC> statsnp;
	StatsnpC stsnp;
	bool trouve;
	statsnp.resize(0);
	int catsnp;
	k=0;
	for (gr=1;gr<=this->ngroupes;gr++) {
// COMPTAGE DES STAT
		nstatgr=0;
		if (groupe[gr].type==0) {   //MICROSAT
			nstatgr +=4*nsamp;   //"NAL","HET","VAR","MGW"
			if (nsamp>1) nstatgr +=8*nsamp*(nsamp-1)/2;  //"N2P","H2P","V2P","FST","LIK","DAS","DM2"
			if (nsamp>2) nstatgr +=nsamp*(nsamp-1)*(nsamp-2)/2;   //AML
		} else if (groupe[gr].type==1) {   //DNA SEQUENCE
			nstatgr +=8*nsamp;   //"NHA","NSS","MPD","VPD","DTA","PSS","MNS","VNS"
			if (nsamp>1) nstatgr +=5*nsamp*(nsamp-1)/2;  //"NH2","NS2","MP2","MPB","HST"
			if (nsamp>2) nstatgr +=nsamp*(nsamp-1)*(nsamp-2)/2;   //"SML"
		} else if (groupe[gr].type==2) {   //SNP
			nstatgr +=4*nsamp;   //"HP0","HM1","HV1","HMO"
			if (nsamp>1) nstatgr +=8*nsamp*(nsamp-1)/2;  
			//"NP0","NM1","NV1","NMO","FP0","FM1","FV1","FMO"
			if (nsamp>2) nstatgr +=4*nsamp*(nsamp-1)*(nsamp-2)/2;   //"AP0","AM1","AV1","AMO"
		}
		this->nstat += nstatgr;
		groupe[gr].nstat=nstatgr;//cout<<"gr="<<gr<<"   nstat="<<nstatgr<<"\n";
//DEFINITION DES STAT
		groupe[gr].sumstat = vector<StatC>(nstatgr);//cout<<"dimensionnement\n";
		//cout<<"type du groupe : "<<groupe[gr].type<<"\n";
		if (groupe[gr].type==0) { //MICROSAT 
			for (int i=1;i<=4;i++) {
				for (int sa=1;sa<=nsamp;sa++){
					groupe[gr].sumstat[k].cat=i;
					groupe[gr].sumstat[k].samp=sa;
					s1=stat_type[i]+"_"+IntToString(gr)+"_"+IntToString(sa);
					this->entetestat +=centre(s1,14);
					//cout<<"k="<<k<<"  "<<this->entetestat<<"\n";
					k++;
				}
			}
			if (nsamp>1) {
				for (int i=5;i<=11;i++) {
					for (int sa=1;sa<=nsamp;sa++) {
						for (int sa1=sa+1;sa1<=nsamp;sa1++) {
							groupe[gr].sumstat[k].cat=i;
							groupe[gr].sumstat[k].samp=sa;
							groupe[gr].sumstat[k].samp1=sa1;
							s1=stat_type[i]+"_"+IntToString(gr)+"_"+IntToString(sa)+"&"+IntToString(sa1);
							this->entetestat +=centre(s1,14);
					//cout<<"k="<<k<<"  "<<this->entetestat<<"\n";
							k++;
							if (i==9){
								groupe[gr].sumstat[k].cat=i;
								groupe[gr].sumstat[k].samp=sa1;
								groupe[gr].sumstat[k].samp1=sa;
								s1=stat_type[i]+"_"+IntToString(gr)+"_"+IntToString(sa1)+"&"+IntToString(sa);
								this->entetestat +=centre(s1,14);
								k++;
							}
						}
					}
				}
			}
			if (nsamp>2) {
				for (int i=12;i<=12;i++) {
					for (int sa=1;sa<=nsamp;sa++) {
						for (int sa1=1;sa1<=nsamp;sa1++) {
						    for (int sa2=1;sa2<=nsamp;sa2++) {
						        if ((sa1!=sa)and(sa2!=sa)and(sa2>sa1)) {
									groupe[gr].sumstat[k].cat=i;
									groupe[gr].sumstat[k].samp=sa;
									groupe[gr].sumstat[k].samp1=sa1;
									groupe[gr].sumstat[k].samp2=sa2;
									s1=stat_type[i]+"_"+IntToString(gr)+"_"+IntToString(sa)+"&"+IntToString(sa1)+"&"+IntToString(sa2);
									this->entetestat +=centre(s1,14);
									k++;
								}
							}
						}
					}
				}
			}
		}
		if (groupe[gr].type==1) {   //DNA SEQUENCE
			for (int i=-1;i>=-8;i--) {
				for (int sa=1;sa<=nsamp;sa++){
					groupe[gr].sumstat[k].cat=i;
					groupe[gr].sumstat[k].samp=sa;
					j=12;do {j++;}while(i!=stat_num[j]);
					s1=stat_type[j]+"_"+IntToString(gr)+"_"+IntToString(sa);
					this->entetestat +=centre(s1,14);
					
					k++;
				}
			}
			if (nsamp>1) {
				for (int i=-9;i>=-13;i--) {
					for (int sa=1;sa<=nsamp;sa++) {
						for (int sa1=sa+1;sa1<=nsamp;sa1++) {
							groupe[gr].sumstat[k].cat=i;
							groupe[gr].sumstat[k].samp=sa;
							groupe[gr].sumstat[k].samp1=sa1;
							j=12;do {j++;}while(i!=stat_num[j]);
							s1=stat_type[j]+"_"+IntToString(gr)+"_"+IntToString(sa)+"&"+IntToString(sa1);
							this->entetestat +=centre(s1,14);
							k++;
						}
					}
				}
			}
			if (nsamp>2) {
				for (int i=-14;i>=-14;i--) {
					for (int sa=1;sa<=nsamp;sa++) {
						for (int sa1=1;sa1<=nsamp;sa1++) {
						    for (int sa2=1;sa2<=nsamp;sa2++) {
						        if ((sa1!=sa)and(sa2!=sa)and(sa2>sa1)) {
									groupe[gr].sumstat[k].cat=i;
									groupe[gr].sumstat[k].samp=sa;
									groupe[gr].sumstat[k].samp1=sa1;
									groupe[gr].sumstat[k].samp2=sa2;
									j=12;do {j++;}while(i!=stat_num[j]);
									s1=stat_type[j]+"_"+IntToString(gr)+"_"+IntToString(sa)+"&"+IntToString(sa1)+"&"+IntToString(sa2);
									this->entetestat +=centre(s1,14);
									k++;
								}
							}
						}
					}
				}
			}
		}
		if (groupe[gr].type==2) {   //SNP
			for (int i=21;i<=24;i++) {
				catsnp = (i-21)/4;
				for (int sa=1;sa<=nsamp;sa++){
					groupe[gr].sumstat[k].cat=i;
					groupe[gr].sumstat[k].samp=sa;
					j=25;do {j++;}while(i!=stat_num[j]);
					s1=stat_type[j]+"_"+IntToString(gr)+"_"+IntToString(sa);
					this->entetestat +=centre(s1,14);
					trouve=false;
					if (statsnp.size()>0){
						for (int jj=0; jj < (int)statsnp.size();jj++) {
							trouve=((statsnp[jj].cat==catsnp)and(statsnp[jj].samp==groupe[gr].sumstat[k].samp));
							if (trouve) {groupe[gr].sumstat[k].numsnp=jj;break;}
						}
					}
					if (not trouve) {
						stsnp.cat=catsnp;
						stsnp.samp=groupe[gr].sumstat[k].samp;
						stsnp.defined=false;
						groupe[gr].sumstat[k].numsnp=statsnp.size();
						statsnp.push_back(stsnp);
					}
					k++;
				}
			} 
			//cout<<"fin des sumstat 21 à 24 statsnp.size ="<<statsnp.size()<<"\n";
			if (nsamp>1) {
				for (int i=25;i<=32;i++) {
				    catsnp = (i-21)/4;
					for (int sa=1;sa<=nsamp;sa++) {
						for (int sa1=sa+1;sa1<=nsamp;sa1++) {
							groupe[gr].sumstat[k].cat=i;
							groupe[gr].sumstat[k].samp=sa;
							groupe[gr].sumstat[k].samp1=sa1;
							j=25;do {j++;}while(i!=stat_num[j]);
							s1=stat_type[j]+"_"+IntToString(gr)+"_"+IntToString(sa)+"&"+IntToString(sa1);
							this->entetestat +=centre(s1,14);
							trouve=false;
							if (statsnp.size()>0){
								for (int jj=0; jj < (int)statsnp.size(); jj++) {
									trouve=((statsnp[jj].cat==catsnp)and(statsnp[jj].samp==groupe[gr].sumstat[k].samp)and(statsnp[jj].samp1==groupe[gr].sumstat[k].samp1));
									if (trouve) {groupe[gr].sumstat[k].numsnp=jj;break;}
								}
							}
							if (not trouve) {
								stsnp.cat=catsnp;
								stsnp.samp=groupe[gr].sumstat[k].samp;
								stsnp.samp1=groupe[gr].sumstat[k].samp1;
								stsnp.defined=false;
								groupe[gr].sumstat[k].numsnp=statsnp.size();
								statsnp.push_back(stsnp);
							}
							k++;
						}
					}
				}
				//cout<<"fin des sumstat 25 à 32 statsnp.size ="<<statsnp.size()<<"\n";
			}
			if (nsamp>2) {
				for (int i=33;i<=36;i++) {
					catsnp = (i-21)/4;
					for (int sa=1;sa<=nsamp;sa++) {
						for (int sa1=1;sa1<=nsamp;sa1++) {
						    for (int sa2=1;sa2<=nsamp;sa2++) {
						        if ((sa1!=sa)and(sa2!=sa)and(sa2>sa1)) {
								    groupe[gr].sumstat[k].cat=i;
								    groupe[gr].sumstat[k].samp=sa;
								    groupe[gr].sumstat[k].samp1=sa1;
								    groupe[gr].sumstat[k].samp2=sa2;
								    j=25;do {j++;}while(i!=stat_num[j]);
								    s1=stat_type[j]+"_"+IntToString(gr)+"_"+IntToString(sa)+"&"+IntToString(sa1)+"&"+IntToString(sa2);
								    this->entetestat +=centre(s1,14);
								    //cout<<"samples "<<sa<<", "<<sa1<<" et "<<sa2<<"\n";
								    trouve=false;
								    if (statsnp.size()>0){
									    for (int jj = 0; jj < (int)statsnp.size(); jj++) {
										    trouve=((statsnp[jj].cat==catsnp)and(statsnp[jj].samp==groupe[gr].sumstat[k].samp)and(statsnp[jj].samp1==groupe[gr].sumstat[k].samp1)and(statsnp[jj].samp2==groupe[gr].sumstat[k].samp2));
										    if (trouve) {groupe[gr].sumstat[k].numsnp=jj;break;}
									    }
								    }
								    if (not trouve) {
									    stsnp.cat=catsnp;
									    stsnp.samp=groupe[gr].sumstat[k].samp;
									    stsnp.samp1=groupe[gr].sumstat[k].samp1;
									    stsnp.samp2=groupe[gr].sumstat[k].samp2;
									    stsnp.defined=false;
								    	groupe[gr].sumstat[k].numsnp=statsnp.size();
									    statsnp.push_back(stsnp);
								    }
								    k++;
								}
							}
						}
					}
				}
			} //cout<<"fin des sumstat 33 à 36 statsnp.size ="<<statsnp.size()<<"\n";
		}
		groupe[gr].nstatsnp=statsnp.size();
		//cout<<"groupe[gr].nstatsnp="<<groupe[gr].nstatsnp<<"\n";
		if (groupe[gr].nstatsnp>0){
			groupe[gr].sumstatsnp = vector<StatsnpC>(groupe[gr].nstatsnp);
			for (int i=0;i<groupe[gr].nstatsnp;i++){
				groupe[gr].sumstatsnp[i].cat=statsnp[i].cat;
				groupe[gr].sumstatsnp[i].samp=statsnp[i].samp;
				groupe[gr].sumstatsnp[i].samp1=statsnp[i].samp1;
				groupe[gr].sumstatsnp[i].samp2=statsnp[i].samp2;
				groupe[gr].sumstatsnp[i].defined=statsnp[i].defined;
				//cout<<"recopie de sumstatsnp["<<i<<"]\n";
			}
			//cout<<"fin de la boucle\n";
			statsnp.clear();
			//cout<<"apres statsnp.clear()\n";
		}
		//cout<<"fin du traitement du groupe "<<gr<<"\n";
	}
	splitwords(this->entetestat, " ", statname);
	//for(int i=0;i<(int)statname.size();i++) cout<<this->statname[i]<<"  ";cout<<"\n";
	//cout<<"nstat="<<this->nstat<<"      nss="<<nss<<"\n";
	return 0;
}

int HeaderC::readHeaderGroupStat(ifstream & file) {
	string s1;
	vector<string> ss, ss1;
	int j, k, nss, gr;
	//Partie group statistics
	if (debuglevel==2) cout<<"debut des group stat\n";
	this->nstat=0;
	getline(file,s1);       //ligne vide
	getline(file,s1);       //ligne "group group statistics"
	//cout <<"s1="<<s1<<"\n";
	for (gr=1;gr<=this->ngroupes;gr++) {
		getline(file,s1);
		//cout <<"s1="<<s1<<"\n";
		splitwords(s1," ",ss); nss = ss.size();

		//cout <<"s1="<<s1<<"   ss[3]="<< ss[3] <<"   atoi = "<< atoi(ss[3].c_str()) <<"\n";
		groupe[gr].nstat = getwordint(ss[nss-1],0);
		this->nstat +=groupe[gr].nstat;
		groupe[gr].sumstat = vector<StatC>(groupe[gr].nstat);
		k=0;
		vector <StatsnpC> statsnp;
		StatsnpC stsnp;
		bool trouve;
		statsnp.resize(0);
		int catsnp;
		if (debuglevel==2) cout<<"nstat="<<groupe[gr].nstat<<"\n";
		while (k<groupe[gr].nstat) {
			if (debuglevel==2) cout <<"stat "<<k<<"    groupe "<<gr<<"\n";
			getline(file,s1);
			if (debuglevel==2) cout <<"s1="<<s1<<"\n";
			splitwords(s1," ",ss); nss = ss.size();
			if (debuglevel==2) cout << ss[0]<<"\n";
			j=0;while (ss[0]!=stat_type[j]) {
				//if (debuglevel==2)	cout << ss[0] << " " << j << " "<< stat_type[j] << endl;
				j++;
			}
			if (debuglevel==2) cout<<"j="<<j<<"\n";
			if (groupe[gr].type==0) {   //MICROSAT
				if (stat_num[j]<5) {
					for (int i=1;i<nss;i++) {
						groupe[gr].sumstat[k].cat=stat_num[j];
						groupe[gr].sumstat[k].samp=atoi(ss[i].c_str());
						k++;
					}
				} else if (stat_num[j]<12) {
					for (int i=1;i<nss;i++) {
						groupe[gr].sumstat[k].cat=stat_num[j];
						splitwords(ss[i],"&",ss1);
						groupe[gr].sumstat[k].samp=atoi(ss1[0].c_str());
						groupe[gr].sumstat[k].samp1=atoi(ss1[1].c_str());
						k++;
					}
				} else if (stat_num[j]==12) {
					for (int i=1;i<nss;i++) {
						groupe[gr].sumstat[k].cat=stat_num[j];
						splitwords(ss[i],"&",ss1);
						groupe[gr].sumstat[k].samp=atoi(ss1[0].c_str());
						groupe[gr].sumstat[k].samp1=atoi(ss1[1].c_str());
						groupe[gr].sumstat[k].samp2=atoi(ss1[2].c_str());
						k++;
					}
				}
			} else if (groupe[gr].type==1) {   //DNA SEQUENCE
				if (stat_num[j]>-9) {
					for (int i=1;i<nss;i++) {
						groupe[gr].sumstat[k].cat=stat_num[j];
						groupe[gr].sumstat[k].samp=atoi(ss[i].c_str());
						k++;
					}
				} else if (stat_num[j]>-14) {
					for (int i=1;i<nss;i++) {
						groupe[gr].sumstat[k].cat=stat_num[j];
						splitwords(ss[i],"&",ss1);
						groupe[gr].sumstat[k].samp=atoi(ss1[0].c_str());
						groupe[gr].sumstat[k].samp1=atoi(ss1[1].c_str());
						k++;
					}
				} else if (stat_num[j]==-14) {
					for (int i=1;i<nss;i++) {
						groupe[gr].sumstat[k].cat=stat_num[j];
						splitwords(ss[i],"&",ss1); 
						groupe[gr].sumstat[k].samp=atoi(ss1[0].c_str());
						groupe[gr].sumstat[k].samp1=atoi(ss1[1].c_str());
						groupe[gr].sumstat[k].samp2=atoi(ss1[2].c_str());
						k++;
					}
				}
			} else if (groupe[gr].type==2) {   //SNP
				//cout<<"statnum="<<stat_num[j]<<"\n";
				if (stat_num[j]<50) {
					catsnp = (stat_num[j]-21)/4;
					if (debuglevel==2) cout<<"stat_num["<<j<<"]="<<stat_num[j]<<"   catsnp="<<catsnp<<"\n";
				}
				if (stat_num[j]<25) {
					if (debuglevel==2) cout<<"nss="<<nss<<"\n";
					for (int i=1;i<nss;i++) {
						if (debuglevel==2) cout<<"i="<<i<<"\n";
						groupe[gr].sumstat[k].cat=stat_num[j];
						groupe[gr].sumstat[k].samp=atoi(ss[i].c_str());
						trouve=false;
						if (debuglevel==2) cout<<"statsnp.size="<<statsnp.size()<<"\n";
						if (statsnp.size()>0){
							for (int jj=0; jj < (int)statsnp.size();jj++) {
								trouve=((statsnp[jj].cat==catsnp)and(statsnp[jj].samp==groupe[gr].sumstat[k].samp));
								if (trouve) {groupe[gr].sumstat[k].numsnp=jj;break;}
							}
						}
						if (debuglevel==2) cout<<"trouve="<<trouve<<"\n";
						if (not trouve) {
							stsnp.cat=catsnp;
							stsnp.samp=groupe[gr].sumstat[k].samp;
							stsnp.defined=false;
							groupe[gr].sumstat[k].numsnp=statsnp.size();
							//cout<<"avant le push_back\n";
							statsnp.push_back(stsnp);
							//cout<<"apres le push_back\n";
						}
						//cout<<"                          numsnp = "<<groupe[gr].sumstat[k].numsnp<<"\n";
						k++;
					}
				} else if ((stat_num[j]>24)and(stat_num[j]<33)) {
					for (int i=1;i<nss;i++) {
						groupe[gr].sumstat[k].cat=stat_num[j];
						splitwords(ss[i],"&",ss1);
						groupe[gr].sumstat[k].samp=atoi(ss1[0].c_str());
						groupe[gr].sumstat[k].samp1=atoi(ss1[1].c_str());
						trouve=false;
						if (statsnp.size()>0){
							for (int jj=0; jj < (int)statsnp.size(); jj++) {
								trouve=((statsnp[jj].cat==catsnp)and(statsnp[jj].samp==groupe[gr].sumstat[k].samp)and(statsnp[jj].samp1==groupe[gr].sumstat[k].samp1));
								if (trouve) {groupe[gr].sumstat[k].numsnp=jj;break;}
							}
						}
						//cout<<"statsnp.size = "<<statsnp.size()<<"   trouve = "<<trouve<<"\n";
						if (not trouve) {
							stsnp.cat=catsnp;
							stsnp.samp=groupe[gr].sumstat[k].samp;
							stsnp.samp1=groupe[gr].sumstat[k].samp1;
							stsnp.defined=false;
							groupe[gr].sumstat[k].numsnp=statsnp.size();
							statsnp.push_back(stsnp);
						}
						//cout<<"                      numsnp = "<<groupe[gr].sumstat[k].numsnp<<"\n";
						k++;
					}
				} else if ((stat_num[j]>32)and (stat_num[j]<50)) {
					for (int i=1;i<nss;i++) {
						groupe[gr].sumstat[k].cat=stat_num[j];
						splitwords(ss[i],"&",ss1);
						groupe[gr].sumstat[k].samp=atoi(ss1[0].c_str());
						groupe[gr].sumstat[k].samp1=atoi(ss1[1].c_str());
						groupe[gr].sumstat[k].samp2=atoi(ss1[2].c_str());
						trouve=false;
						if (statsnp.size()>0){
							for (int jj = 0; jj < (int)statsnp.size(); jj++) {
								trouve=((statsnp[jj].cat==catsnp)and(statsnp[jj].samp==groupe[gr].sumstat[k].samp)and(statsnp[jj].samp1==groupe[gr].sumstat[k].samp1)and(statsnp[jj].samp2==groupe[gr].sumstat[k].samp2));
								if (trouve) {groupe[gr].sumstat[k].numsnp=jj;break;}
							}
						}
						if (not trouve) {
							stsnp.cat=catsnp;
							stsnp.samp=groupe[gr].sumstat[k].samp;
							stsnp.samp1=groupe[gr].sumstat[k].samp1;
							stsnp.samp2=groupe[gr].sumstat[k].samp2;
							stsnp.defined=false;
							groupe[gr].sumstat[k].numsnp=statsnp.size();
							statsnp.push_back(stsnp);
						}
						//cout<<"                      numsnp = "<<groupe[gr].sumstat[k].numsnp<<"\n";
						k++;
					} 
				} else if (stat_num[j]==50) {
					groupe[gr].sumstat[k].cat=stat_num[j];
					k++;
					cout<<"k="<<k<<"\n";
				}
			}
			if (debuglevel==2) cout<<"fin de la stat "<<k<<"\n";
		}
		groupe[gr].nstatsnp=statsnp.size();
		//cout<<"groupe[gr].nstatsnp="<<groupe[gr].nstatsnp<<"\n";
		if (groupe[gr].nstatsnp>0){
			groupe[gr].sumstatsnp = vector<StatsnpC>(groupe[gr].nstatsnp);
			for (int i=0;i<groupe[gr].nstatsnp;i++){
				groupe[gr].sumstatsnp[i].cat=statsnp[i].cat;
				groupe[gr].sumstatsnp[i].samp=statsnp[i].samp;
				groupe[gr].sumstatsnp[i].samp1=statsnp[i].samp1;
				groupe[gr].sumstatsnp[i].samp2=statsnp[i].samp2;
				groupe[gr].sumstatsnp[i].defined=statsnp[i].defined;
			}
			statsnp.clear();
		}
		//for (int i=0;i<groupe[gr].nstat;i++) cout<<groupe[gr].sumstat[i].cat<<"   "<<groupe[gr].sumstat[i].numsnp<<"\n";
	}
	if (debuglevel==2) cout<<"header.txt : fin de la lecture des summary stats\n";
	return 0;
}



int HeaderC::buildSuperScen(){
	//cout<<"avant superscen\n";
	this->scen.nevent=0;
	for (int i=0;i<this->nscenarios;i++) {if (this->scen.nevent<scenario[i].nevent) this->scen.nevent=scenario[i].nevent;}
	this->scen.nn0=0;
	for (int i=0;i<this->nscenarios;i++) {if (this->scen.nn0<scenario[i].nn0)       this->scen.nn0 =scenario[i].nn0;}
	this->scen.nsamp=0;
	for (int i=0;i<this->nscenarios;i++) {if (this->scen.nsamp<scenario[i].nsamp)   this->scen.nsamp=scenario[i].nsamp;}
	this->scen.nparam=0;
	for (int i=0;i<this->nscenarios;i++) {if (this->scen.nparam<scenario[i].nparam)  this->scen.nparam=scenario[i].nparam;}
	this->scen.nparamvar=0;
	for (int i=0;i<this->nscenarios;i++) {if (this->scen.nparamvar<scenario[i].nparamvar) this->scen.nparamvar=scenario[i].nparamvar;}
	this->scen.nconditions=0;
	for (int i=0;i<this->nscenarios;i++) {if (this->scen.nconditions<scenario[i].nconditions) this->scen.nconditions=scenario[i].nconditions;}
	this->scen.event = vector<EventC>(this->scen.nevent);
	this->scen.ne0   = vector<Ne0C>(this->scen.nn0);

	this->scen.time_sample = vector<int>(this->scen.nsamp);
	this->scen.stime_sample = vector<string>(this->scen.nsamp);
	this->scen.histparam = vector<HistParameterC>(this->scen.nparam);
	this->scen.paramvar = vector<double>(this->scen.nparamvar+3);
	if (this->scen.nconditions>0) this->scen.condition = vector<ConditionC>(this->scen.nconditions);
	for (int i=0;i<this->scen.nsamp;i++) this->scen.time_sample[i]=0;
	for (int i=0;i<this->scen.nsamp;i++) this->scen.stime_sample[i]=" ";
	for (int i=0;i<this->scen.nparamvar+3;i++) this->scen.paramvar[i]=-1;
	PriorC pr;pr.loi="uniforme";pr.mini=0.0;pr.maxi=1.0;pr.ndec=3;pr.constant=false;
	HistParameterC pp;pp.name="bidon";pp.category=0;pp.value=-1; pp.prior = pr; //pp.prior=copyprior(pr);
	for (int i=0;i<this->scen.nparam;i++)  this->scen.histparam[i] = pp;
	//this->scen.ecris();
	//cout<<"apres superscen\n";
	if (debuglevel==2) cout<<"header.txt : fin de l'établissement du superscen\n";
	//for(int i=0;i<this->nscenarios;i++) scenario[i].ecris();
	return 0;
}

int HeaderC::buildMutParam(){
	int gr;
	this->nparamut=0;
	for (gr=1;gr<=this->ngroupes;gr++) {
		if (groupe[gr].type==0) {
			if (not groupe[gr].priormutmoy.constant) this->nparamut++;
			if (not groupe[gr].priorPmoy.constant) this->nparamut++;
			if (not groupe[gr].priorsnimoy.constant) this->nparamut++;
		} else if (groupe[gr].type==1) {
			if (not groupe[gr].priormusmoy.constant) this->nparamut++;
			//cout<<"calcul de nparamut\n";
			if (groupe[gr].mutmod>0) {if (not groupe[gr].priork1moy.constant) this->nparamut++;}
			if (groupe[gr].mutmod==3){if (not groupe[gr].priork2moy.constant) this->nparamut++;}
			//cout<<"fin du calcul de nparamut = "<<this->nparamut<<"\n";
		}
	}
	this->mutparam = vector<MutParameterC>(this->nparamut);
	this->nparamut=0;
	for (gr=1;gr<=this->ngroupes;gr++) {
		if (groupe[gr].type==0) {
			if (not groupe[gr].priormutmoy.constant) {
				this->mutparam[nparamut].name="µmic_"+IntToString(gr);
				this->mutparam[nparamut].groupe=gr;
				this->mutparam[nparamut].category=0;
				//this->mutparam[nparamut].prior = copyprior(groupe[gr].priormutmoy);
				this->mutparam[nparamut].prior = groupe[gr].priormutmoy;
				this->nparamut++;
			}
			if (not groupe[gr].priorPmoy.constant) {
				this->mutparam[nparamut].name="pmic_"+IntToString(gr);
				this->mutparam[nparamut].groupe=gr;
				this->mutparam[nparamut].category=1;
				//this->mutparam[nparamut].prior = copyprior(groupe[gr].priorPmoy);
				this->mutparam[nparamut].prior = groupe[gr].priorPmoy;
				this->nparamut++;
			}
			if (not groupe[gr].priorsnimoy.constant) {
				this->mutparam[nparamut].name="snimic_"+IntToString(gr);
				this->mutparam[nparamut].groupe=gr;
				this->mutparam[nparamut].category=2;
				// this->mutparam[nparamut].prior = copyprior(groupe[gr].priorsnimoy);
				this->mutparam[nparamut].prior = groupe[gr].priorsnimoy;
				this->nparamut++;
			}
		} else if (groupe[gr].type==1) {
			if (not groupe[gr].priormusmoy.constant){
				this->mutparam[nparamut].name="useq_"+IntToString(gr);
				this->mutparam[nparamut].groupe=gr;
				this->mutparam[nparamut].category=3;
				// this->mutparam[nparamut].prior = copyprior(groupe[gr].priormusmoy);
				this->mutparam[nparamut].prior = groupe[gr].priormusmoy;
				this->nparamut++;
			}
			if ((groupe[gr].mutmod>0)and(not groupe[gr].priork1moy.constant)) {
				this->mutparam[nparamut].name="k1seq_"+IntToString(gr);
				this->mutparam[nparamut].groupe=gr;
				this->mutparam[nparamut].category=4;
				// this->mutparam[nparamut].prior = copyprior(groupe[gr].priork1moy);
				this->mutparam[nparamut].prior = groupe[gr].priork1moy;
				this->nparamut++;
			}
			if ((groupe[gr].mutmod==3)and(not groupe[gr].priork2moy.constant)) {
				this->mutparam[nparamut].name="k2seq_"+IntToString(gr);
				this->mutparam[nparamut].groupe=gr;
				this->mutparam[nparamut].category=5;
				// this->mutparam[nparamut].prior = copyprior(groupe[gr].priork2moy);
				this->mutparam[nparamut].prior = groupe[gr].priork2moy;
				this->nparamut++;
			}
		}
	}
	if (debuglevel==2) cout<<"header.txt : fin de l'établissement des paramètres mutationnels   nparamut="<<this->nparamut<<"\n";
	return 0;
}

int HeaderC::readHeaderEntete(ifstream & file){
	string s1;
	vector<string> ss;
	int nss;
	//Entete du fichier reftable
	getline(file,s1);       //ligne vide
	getline(file,this->entete);     //ligne entete
	this->statname = vector<string>(this->nstat);
	splitwords(this->entete," ",ss); nss = ss.size();
	
	if (debuglevel==2) cout<<"nss="<<nss<<"   nparamtot="<<nparamtot<<"   nparamut="<<nparamut<<"   nstat="<<nstat<<"\n";
	for (int i=nstat-1;i>=0;i--) this->statname[i]=ss[--nss];
	this->entetehist=this->entete.substr(0,this->entete.length()-14*(nparamut+nstat));
	if (nparamut>0) this->entetemut=this->entete.substr(this->entetehist.length(),14*nparamut);else this->entetemut="";
	this->entetestat=this->entete.substr(this->entetehist.length()+this->entetemut.length(),14*nstat);
	if (debuglevel==2) cout<<"les trois entetes\n";
	if (debuglevel==2) cout<<">>>"<<this->entetehist<<"<<<\n";
	if (debuglevel==2) cout<<">>>"<<this->entetemut<<"<<<\n";
	if (debuglevel==2) cout<<">>>"<<this->entetestat<<"<<<\n";
	//for (int i=0;i<this->nstat;i++) cout<<this->statname[i]<<"   ";cout<<"\n";
	//cout<<"scenarios à la fin de readheader\n";
	return 0;
}

int HeaderC::readHeader(string headerfilename){
	string s1; // string s2,**sl,*ss,*ss1,*ss2;
	int error = 0; // int *nlscen,nss,nss1,j,k,gr,grm,k1,cat,nl=0;
	//cout<<"debut de readheader\n";
	//cout<<"readHeader headerfilename = "<<headerfilename<<"\n";
	ifstream file(headerfilename.c_str(), ios::in);
	if (!file.is_open()) {
		this->message = "Header file "+headerfilename+" not found";
		cout<<this->message<<"\n";
		throw std::runtime_error(this->message);
		return 1;
	} else this->message="";

	error = readHeaderDebut(file); cout<<"----------------------------------------apres readHeaderDebut\n";
	if(error != 0) return error;

	error = readHeaderScenarios(file);cout<<"----------------------------------------apres readHeaderScenarios\n";
	if(error != 0) return error;

	error = readHeaderHistParam(file);cout<<"----------------------------------------apres readHeaderHistParam\n";
	if(error != 0) return error;

	error = readHeaderLoci(file);cout<<"----------------------------------------apres readHeaderLoci\n";
	if(error != 0) return error;

	error = readHeaderGroupPrior(file);cout<<"----------------------------------------apres readHeaderGroupPrior\n";
	if(error != 0) return error;

	error = buildSuperScen();cout<<"----------------------------------------apres buildSuperScen\n";
	if(error != 0) return error;
	
	if (not randomforest) {error = readHeaderGroupStat(file);cout<<"----------------------------------------apres readHeaderGroupStat\n";}
	else {error =readHeaderAllStat(file);cout<<"----------------------------------------apres readHeaderAllStat\n";}
	if(error != 0) return error;

	error = buildMutParam();cout<<"----------------------------------------apres buildMutParam\n";
	if(error != 0) return error;

	if (not randomforest){
	error = readHeaderEntete(file);cout<<"----------------------------------------apres readHeaderEntete\n";
	if(error != 0) return error;}
    /*this->threshold = 0.;
	if (not file.eof()){
		getline(file,s1);
		getline(file,s1);
		this->reffreqmin = getwordfloat(s1,0);
		getline(file, s1);
		particuleobs.locpol = getwordfloat(s1,0);
		getline(file, s1);
		this->threshold = getwordfloat(s1, 0);
		cout <<"refreqmin="<<this->reffreqmin<<"   locpol="<<particuleobs.locpol<<"   threshold="<<this->threshold<<"\n";
	}*/
	cout<<"fin de readheader  nloc retenus="<<dataobs.nloc<<"\n";
	return 0;
}

int HeaderC::readHeadersimDebut(ifstream & file){
	string s1;
	vector<string> ss;
	int nss,*nf,*nm;
	if (debuglevel==2) cout<<"debut de readheadersim\n";
	getline(file, s1);
	splitwords(s1," ",ss); nss = ss.size();
	this->datafilename=ss[0];
	this->nsimfile=atoi(ss[1].c_str());
	if (nss>2) datasim.sexratio=atof(ss[2].c_str()); else datasim.sexratio=0.5;
	cout<<"sexratio="<<datasim.sexratio<<"\n";
	getline(file,s1);
	splitwords(s1," ",ss); nss = ss.size();
	datasim.nsample=atoi(ss[0].c_str());
	datasim.nind.resize(datasim.nsample);
	datasim.indivsexe.resize(datasim.nsample);
	nf=new int[datasim.nsample];
	nm=new int[datasim.nsample];
	for (int i=0;i<datasim.nsample;i++) {
		getline(file,s1);
		splitwords(s1," ",ss); nss = ss.size();
		nf[i]=atoi(ss[0].c_str());
		nm[i]=atoi(ss[1].c_str());
		datasim.nind[i]=nf[i]+nm[i];
		datasim.indivsexe[i].resize(datasim.nind[i]);
		for (int j=0;j<nf[i];j++) datasim.indivsexe[i][j]=2;
		for (int j=nf[i];j<datasim.nind[i];j++) datasim.indivsexe[i][j]=1;
	}
	cout<<"nom générique : "<<this->datafilename<<"\n";
	cout<<"nombre de fichiers à simuler = "<<this->nsimfile<<"\n";
	cout<<"nombre d'échantillons = "<<datasim.nsample<<"\n";
	for (int i=0;i<datasim.nsample;i++) {
		cout<<"Echantillon "<<i+1<<"  :  "<<nf[i]<<" femelles et "<<datasim.nind[i]-nf[i]<<" males\n";
	}
	datasim.nmisshap=0;
	datasim.nmissnuc=0;
	return 0;
}

int HeaderC::readHeadersimScenario(ifstream & file){
	string s1;
	vector<string> ss;
	int nlscen;
	getline(file,s1);   //cout<<s1<<"\n";     //ligne vide
	getline(file,s1);    //cout<<s1<<"\n";    // nombre de scenarios
	this->nscenarios=1; //cout<<this->nscenarios<<"\n";
	splitwords(s1," ",ss);
	nlscen=getwordint(ss[1],1);
	scenario = vector<ScenarioC>(this->nscenarios);
	vector<string>sl (nlscen);
	scenario[0].number = 1;
	scenario[0].prior_proba = 1.0;
	scenario[0].nparam = 0;
	scenario[0].nparamvar=0;
	for (int j=0;j<nlscen;j++) {getline(file,sl[j]);/*cout<<sl[j]<<"\n";*/}
	scenario[0].read_events(nlscen, &(sl[0]));
	//cout<<"apres read_events\n";
	//cout<<"dans readHeadersimScenario\n";
	//scenario[0].ecris();
	if (debuglevel==2) cout<<"header.txt : fin de la lecture de la partie scénario(s)\n";
	return 0;
}

int HeaderC::readHeadersimHistParam(std::ifstream & file){
	string s1;
	vector<string> ss;
	
	getline(file,s1);       //ligne vide
	getline(file,s1);
	//cout<<s1<<"\n";
	splitwords(s1," ",ss); 
	this->nparamtot=getwordint(ss[2],1);
	scenario[0].nconditions=0;
	this->histparam = vector<HistParameterC>(this->nparamtot);
	scenario[0].histparam = vector<HistParameterC>(this->nparamtot);
	scenario[0].nparam = this->nparamtot;
	for (int i=0;i<this->nparamtot;i++) {
		getline(file,s1);  //cout<<s1<<"\n";
		splitwords(s1," ",ss);
		this->histparam[i].name=ss[0];
		if (ss[1]=="N") this->histparam[i].category = 0;
		else if  (ss[1]=="T") this->histparam[i].category = 1;
		else if  (ss[1]=="A") this->histparam[i].category = 2;
		this->histparam[i].value=atof(ss[2].c_str());
		scenario[0].histparam[i] = this->histparam[i];
		scenario[0].histparam[i].prior.constant=true;
		scenario[0].histparam[i].prior.fixed=true;
	}
	if (debuglevel==2) cout<<"header.txt : fin de la lecture de la partie priors des paramètres démographiques\n";
	return 0;	
}

int HeaderC::readHeadersimLoci(std::ifstream & file){
	string s1;
	vector<string> ss;
	int grm,gr,kloc,loctyp;
	double som;
	getline(file,s1); //cout<<s1<<"    ligne vide ? \n";       //ligne vide
	getline(file,s1);  //cout<<s1<<"     loci description ?\n";
	splitwords(s1," ",ss);
	datasim.nloc=getwordint(s1,3);
	cout<<"nloc="<<datasim.nloc<<"\n";
	datasim.locus = vector<LocusC>(datasim.nloc);
	grm=1;
	datasim.filetype=0;
	for (int loc=0;loc<datasim.nloc;loc++){
		getline(file,s1);
		splitwords(s1," ",ss);
		datasim.locus[loc].name = strdup(ss[0].c_str());
		if (ss[1]=="<A>") datasim.locus[loc].type =0;
		else if (ss[1]=="<H>")   datasim.locus[loc].type =1;
		else if (ss[1]=="<X>")   datasim.locus[loc].type =2;
		else if (ss[1]=="<Y>")   datasim.locus[loc].type =3;
		else if (ss[1]=="<M>")   datasim.locus[loc].type =4;
		datasim.cal_coeffcoal(loc);
		if (ss[2]=="[M]") {
			s1=ss[3].substr(1);gr=atoi(s1.c_str());datasim.locus[loc].groupe=gr;if (gr>grm) grm=gr;
			datasim.locus[loc].motif_size=atoi(ss[4].c_str());
			datasim.locus[loc].motif_range=atoi(ss[5].c_str());
			datasim.locus[loc].mini=100;datasim.locus[loc].maxi=300;
		}
		else if (ss[2]=="[S]") {
			datasim.locus[loc].type +=5;
			s1=ss[3].substr(1);gr=atoi(s1.c_str());datasim.locus[loc].groupe=gr;if (gr>grm) grm=gr;
			datasim.locus[loc].dnalength=atoi(ss[4].c_str());
			datasim.locus[loc].mutsit.resize(datasim.locus[loc].dnalength);
			som=0.0;
			datasim.locus[loc].pi_A=atof(ss[5].c_str());som +=datasim.locus[loc].pi_A;
			datasim.locus[loc].pi_C=atof(ss[6].c_str());som +=datasim.locus[loc].pi_C;
			datasim.locus[loc].pi_G=atof(ss[7].c_str());som +=datasim.locus[loc].pi_G;
			datasim.locus[loc].pi_T=atof(ss[8].c_str());som +=datasim.locus[loc].pi_T;
			datasim.locus[loc].pi_A /= som;
			datasim.locus[loc].pi_C /= som;
			datasim.locus[loc].pi_G /= som;
			datasim.locus[loc].pi_T /= som;
			//cout<<datasim.locus[loc].dnalength<<"\n";
			//datasim.locus[loc].dnalength=atoi(ss[k+2].c_str());  //inutile variable déjà renseignée
		}
		else if (ss[2]=="[P]") {
			datasim.filetype=1;
			kloc=getwordint(s1,1);
			loctyp=datasim.locus[loc].type+10;
			s1=ss[3].substr(1);gr=atoi(s1.c_str());if (gr>grm) grm=gr;
			this->ngroupes=grm;
			for (int k=0;k<kloc;k++){
				datasim.locus[loc+k].type =loctyp;
				datasim.locus[loc+k].groupe=gr;
				//cout<<"locus "<<loc+k<<"  type="<<datasim.locus[loc+k].type<<"  groupe="<<datasim.locus[loc+k].groupe<<"\n";
				//cout<<"locus "<<loc+k<<"   sexratio="<<datasim.sexratio<<"\n";
				datasim.cal_coeffcoal(loc+k);
				//cout<<"apres cal_coeffcoal\n";
			}
			loc +=kloc-1;
		}
	}
	for (int loc=0;loc<datasim.nloc;loc++) datasim.locus[loc].nsample = datasim.nsample;
	datasim.catexist = vector<bool>(5);
	datasim.ssize.resize(5);
	for (int i=0;i<5;i++) datasim.catexist[i]=false;
	for (int locustype=0;locustype<5;locustype++) {
		if (not datasim.catexist[locustype]) {
			for (int loc=0;loc<datasim.nloc;loc++) {
				if ((datasim.locus[loc].type % 5) == locustype) {
					datasim.ssize[locustype].resize(datasim.nsample);
					datasim.locus[loc].samplesize = vector<int>(datasim.nsample);
					datasim.locus[loc].ploidie = vector < vector <short int > >(datasim.nsample);
					for (int sa=0;sa<datasim.nsample;sa++) {
						datasim.locus[loc].ploidie[sa] = vector <short int >(datasim.nind[sa]);
						datasim.ssize[locustype][sa]=0;
						for (int ind=0;ind<datasim.nind[sa];ind++) {
							//if ((datasim.locus[loc].type==10)or((datasim.locus[loc].type==12)and(datasim.indivsexe[sa][ind]==2))) datasim.ssize[locustype][sa] +=2;
							switch(locustype) {
								case 0 : datasim.ssize[locustype][sa] +=2;
										 datasim.locus[loc].samplesize[sa] +=2;
										 datasim.locus[loc].ploidie[sa][ind] = 2;
								         break;
								case 1 : datasim.ssize[locustype][sa] +=1;
										 datasim.locus[loc].samplesize[sa] +=1;
										 datasim.locus[loc].ploidie[sa][ind] = 1;
								         break;
								case 2 : datasim.ssize[locustype][sa] +=datasim.indivsexe[sa][ind];
										 datasim.locus[loc].samplesize[sa] +=datasim.indivsexe[sa][ind];
										 datasim.locus[loc].ploidie[sa][ind] = datasim.indivsexe[sa][ind];
								         break;
								case 3 : if (datasim.indivsexe[sa][ind]==1) {
											datasim.ssize[locustype][sa] +=1;
											datasim.locus[loc].samplesize[sa] +=1;
											datasim.locus[loc].ploidie[sa][ind] = 1;
										 } else datasim.locus[loc].ploidie[sa][ind] = 0;
								         break;
								case 4 : datasim.ssize[locustype][sa] +=1;
										 datasim.locus[loc].samplesize[sa] +=1;
										 datasim.locus[loc].ploidie[sa][ind] = 1;
								         break;
							}
/*							if( (locustype == 0) or ( (locustype == 2)and(datasim.indivsexe[sa][ind]==2) ) ) {
								datasim.ssize[locustype][sa] +=2;
							} else {
								datasim.ssize[locustype][sa] +=1;
							}*/
						}
					}
					datasim.catexist[locustype] = true;
				}
			}
		}
	}
	if (debuglevel==2) cout<<"header.txt : fin de la lecture de la partie description des locus\n";
	return 0;
}

int HeaderC::readHeadersimGroupPrior(std::ifstream & file){
	string s1;
	vector<string> ss, ss1;
	//int nss,nss1;
	getline(file,s1); //cout<<s1<<"    ligne vide ? \n";      //ligne vide
	getline(file,s1);       //ligne "group prior"
	this->ngroupes=getwordint(s1,2);
	cout<<"header.ngroupes="<<this->ngroupes<<"\n";
	groupe = vector<LocusGroupC>(this->ngroupes+1);
	this->assignloc(0);
	//cout<<"on attaque les groupes : analyse des priors nombre de groupes="<<this->ngroupes <<"\n";
	for (int gr=1;gr<=this->ngroupes;gr++){
		getline(file,s1);
		splitwords(s1," ",ss); 
		this->assignloc(gr);
		if (ss[2]=="[M]") {
			groupe[gr].type=0;
			getline(file,s1);splitwords(s1," ",ss1); groupe[gr].mutmoy = atof(ss1[1].c_str());
			groupe[gr].priormutmoy.fixed=true;groupe[gr].priormutmoy.constant=true;
			getline(file,s1);splitwords(s1," ",ss1); groupe[gr].priormutloc.readprior(ss1[1]);

			getline(file,s1);splitwords(s1," ",ss1); groupe[gr].Pmoy   = atof(ss1[1].c_str());
			groupe[gr].priorPmoy.fixed=true;groupe[gr].priorPmoy.constant=true;
			getline(file,s1);splitwords(s1," ",ss1); groupe[gr].priorPloc.readprior(ss1[1]);

			getline(file,s1);splitwords(s1," ",ss1); groupe[gr].snimoy = atof(ss1[1].c_str());
			groupe[gr].priorsnimoy.fixed=true;groupe[gr].priorsnimoy.constant=true;
			getline(file,s1);splitwords(s1," ",ss1); groupe[gr].priorsniloc.readprior(ss1[1]);
			//cout<<"sniloc  ";groupe[gr].priorsniloc.ecris();

		} else if (ss[2]=="[S]") {
			groupe[gr].type=1;
			getline(file,s1);splitwords(s1," ",ss1); groupe[gr].musmoy = atof(ss1[1].c_str());
			groupe[gr].priormusmoy.fixed=true;groupe[gr].priormusmoy.constant=true;
			getline(file,s1);splitwords(s1," ",ss1); groupe[gr].priormusloc.readprior(ss1[1]);

			getline(file,s1);splitwords(s1," ",ss1);
			if (ss1[0]!="MODEL") {
				groupe[gr].k1moy  = atof(ss1[1].c_str());
				groupe[gr].priork1moy.fixed=true;groupe[gr].priork1moy.constant=true;
				getline(file,s1);splitwords(s1," ",ss1); groupe[gr].priork1loc.readprior(ss1[1]);
				getline(file,s1);splitwords(s1," ",ss1);			
			}
			if (ss1[0]!="MODEL") {
				groupe[gr].k2moy  = atof(ss1[1].c_str());
				groupe[gr].priork2moy.fixed=true;groupe[gr].priork2moy.constant=true;
				getline(file,s1);splitwords(s1," ",ss1); groupe[gr].priork2loc.readprior(ss1[1]);
				getline(file,s1);splitwords(s1," ",ss1); 
			}

			groupe[gr].p_fixe=atof(ss1[2].c_str());groupe[gr].gams=atof(ss1[3].c_str());
			if (ss1[1]=="JK") groupe[gr].mutmod=0;
			else if (ss1[1]=="K2P") groupe[gr].mutmod=1;
			else if (ss1[1]=="HKY") groupe[gr].mutmod=2;
			else if (ss1[1]=="TN") groupe[gr].mutmod=3;
		}
		groupe[gr].nstat=0;
	}
	if (debuglevel==2) cout<<"header.txt : fin de la lecture de la partie définition des groupes\n";
	return 0;
}

int HeaderC::readHeadersimGroupSNP(){
	groupe = vector<LocusGroupC>(this->ngroupes+1);
	this->assignloc(0);
	for (int gr=1;gr<=this->ngroupes;gr++) {
		this->assignloc(gr);
		groupe[gr].type = 2;
	}
	return 0;
}

int HeaderC::readHeadersimFin(){
	//Mise à jour des locus séquences
	MwcGen mwc;
	mwc.randinit(999,time(NULL));
	int nsv,gr,j;
	bool nouveau;
	if (datasim.filetype==0) {
		for (int loc=0;loc<datasim.nloc;loc++){
			gr=datasim.locus[loc].groupe;
			if ((gr>0)and(datasim.locus[loc].type>4)){
				nsv = floor(datasim.locus[loc].dnalength*(1.0-0.01*groupe[gr].p_fixe)+0.5);
				for (int i=0;i<datasim.locus[loc].dnalength;i++) datasim.locus[loc].mutsit[i] = mwc.ggamma3(1.0,groupe[gr].gams);
				int *sitefix;
				sitefix=new int[datasim.locus[loc].dnalength-nsv];
				for (int i=0;i<datasim.locus[loc].dnalength-nsv;i++) {
					if (i==0) sitefix[i]=mwc.rand0(datasim.locus[loc].dnalength);
					else {
						do {
							sitefix[i]=mwc.rand0(datasim.locus[loc].dnalength);
							nouveau=true;j=0;
							while((nouveau)and(j<i)) {nouveau=(sitefix[i]!=sitefix[j]);j++;}
						} while (not nouveau);
					}
					datasim.locus[loc].mutsit[i] = 0.0;
				}
				delete [] sitefix;
				double s=0.0;
				for (int i=0;i<datasim.locus[loc].dnalength;i++) s += datasim.locus[loc].mutsit[i];
				for (int i=0;i<datasim.locus[loc].dnalength;i++) datasim.locus[loc].mutsit[i] /=s;
			}
		}
	}
	else {
		this->nparamut=0;
	}
	if (debuglevel==2) cout<<"header.txt : fin de la mise à jour des locus séquences\n";
	this->scen = scenario[0];
	//this->scen.ecris();
	if (debuglevel==2) cout<<"header.txt : apres la copie du scénario\n";
	return 0;
}

int HeaderC::readHeadersim(string headersimfilename){
	int error;
//	if (debuglevel==2) cout<<"debut de readheadersim\n";
	//cout<<"readHeader headerfilename = "<<headerfilename<<"\n";
	ifstream file(headersimfilename.c_str(), ios::in);
	if (!file.is_open()) {
		this->message = "HeaderSim  File "+string(headersimfilename)+" not found";
		cout<<this->message<<"\n";
		return 1;
	} else this->message="";
	

	error = readHeadersimDebut(file);
	if(error != 0) return error;

	//Partie Scenarios
	error = readHeadersimScenario(file);
	if(error != 0) return error;

	//Partie historical parameters
	error = readHeadersimHistParam(file);
	if(error != 0) return error;
	
	//Partie loci description
	error = readHeadersimLoci(file);
	if(error != 0) return error;

	//Partie group priors
	if (datasim.filetype==0){
		error = readHeadersimGroupPrior(file);
		if(error != 0) return error;
	}
	if (datasim.filetype==1){
		error = readHeadersimGroupSNP();
		if(error != 0) return error;
	}
	
	error = readHeadersimFin();
	return error;
}

string HeaderC::calstatobs(string statobsfilename) {
	stringstream erreur;
	int jstat;
	vector <vector <int> > ast;
	cout<<"calstatobs ("<<statobsfilename<<")\n";
	particuleobs.dnatrue = true;
	particuleobs.nsample = dataobs.nsample0; //cout<<"particuleobs.nsample = "<<particuleobs.nsample<<"\n";
	particuleobs.dnatrue = true;
	//partie GROUPES
	int ngr = groupe.size(); 
	if (debuglevel==2) cout<<"ngr="<<ngr<<"\n";
	particuleobs.ngr = ngr;
	particuleobs.grouplist = vector <LocusGroupC>(ngr);
	for (int i=0;i<ngr;i++) particuleobs.grouplist[i] = groupe[i];
	if (debuglevel==2) cout<<"Dans calstatobs nloc="<<groupe[1].nloc<<"\n";
	/*if (groupe[0].nloc>0){
		particuleobs.grouplist[0].loc  = vector<int>(groupe[0].nloc);
		for (int i=0;i<groupe[0].nloc;i++) particuleobs.grouplist[0].loc[i] = groupe[0].loc[i];
	}*/
	for (int gr=1;gr<ngr;gr++) {
		//particuleobs.grouplist[gr] = groupe[gr];
		for (int i=0;i<groupe[gr].nstatsnp;i++){
			particuleobs.grouplist[gr].sumstatsnp[i].x = vector <long double>(groupe[gr].nloc);
			particuleobs.grouplist[gr].sumstatsnp[i].w = vector <long double>(groupe[gr].nloc);
		}
	}
	//for (int i=0;i<particuleobs.grouplist[1].nstat;i++) cout<<"particuleobs.grouplist[1].sumstat["<<i<<"].numsnp = "<<particuleobs.grouplist[1].sumstat[i].numsnp<<"\n";
/*	//////////////////////////////////////////////
	for (int gr=1;gr<=ngr;gr++) {
		if (debuglevel==2) cout <<"groupe "<<gr<<"\n";
		particuleobs.grouplist[gr].type =groupe[gr].type;
		particuleobs.grouplist[gr].nloc = groupe[gr].nloc;
		particuleobs.grouplist[gr].loc  = new int[groupe[gr].nloc];
		for (int i=0;i<groupe[gr].nloc;i++) particuleobs.grouplist[gr].loc[i] = groupe[gr].loc[i];
		particuleobs.grouplist[gr].nstat=groupe[gr].nstat;
		particuleobs.grouplist[gr].sumstat = new StatC[groupe[gr].nstat];
		//cout<<"calstatobs nstat["<<gr<<"]="<<groupe[gr].nstat<<"\n";
		for (int i=0;i<groupe[gr].nstat;i++){
			particuleobs.grouplist[gr].sumstat[i].cat   = groupe[gr].sumstat[i].cat;
			particuleobs.grouplist[gr].sumstat[i].samp  = groupe[gr].sumstat[i].samp;
			particuleobs.grouplist[gr].sumstat[i].samp1 = groupe[gr].sumstat[i].samp1;
			particuleobs.grouplist[gr].sumstat[i].samp2 = groupe[gr].sumstat[i].samp2;
			particuleobs.grouplist[gr].sumstat[i].numsnp = groupe[gr].sumstat[i].numsnp;
		}
		particuleobs.grouplist[gr].nstatsnp=groupe[gr].nstatsnp;
		if (groupe[gr].nstatsnp>0){
			particuleobs.grouplist[gr].sumstatsnp = new StatsnpC[groupe[gr].nstatsnp];
			for (int i=0;i<groupe[gr].nstatsnp;i++){
				particuleobs.grouplist[gr].sumstatsnp[i].cat   = groupe[gr].sumstatsnp[i].cat;
				particuleobs.grouplist[gr].sumstatsnp[i].samp  = groupe[gr].sumstatsnp[i].samp;
				particuleobs.grouplist[gr].sumstatsnp[i].samp1 = groupe[gr].sumstatsnp[i].samp1;
				particuleobs.grouplist[gr].sumstatsnp[i].samp2 = groupe[gr].sumstatsnp[i].samp2;
				particuleobs.grouplist[gr].sumstatsnp[i].defined = groupe[gr].sumstatsnp[i].defined;
				//particuleobs.grouplist[gr].sumstatsnp[i].sorted = groupe[gr].sumstatsnp[i].sorted;
				particuleobs.grouplist[gr].sumstatsnp[i].x = new long double[groupe[gr].nloc];
				particuleobs.grouplist[gr].sumstatsnp[i].w = new long double[groupe[gr].nloc];
			}
		}
	}
*/	//////////////////////////////////////
	if (debuglevel==2) cout<<"apres GROUPS\n";
	//partie LOCUSLIST
	int kmoy;
	particuleobs.nloc = dataobs.nloc;
	if (debuglevel==2) cout<<"avant la partie locus  nloc="<<dataobs.nloc<<"\n";
	particuleobs.locuslist =vector <LocusC>(dataobs.nloc);
	if (debuglevel==2) cout<<"après l'allocation mémoire des locus\n";
	for (int kloc=0;kloc<dataobs.nloc;kloc++){
		//if (debuglevel==2) cout<<"Locus "<<kloc<<"\n";
		particuleobs.locuslist[kloc].nsample = dataobs.nsample;
		particuleobs.locuslist[kloc] = dataobs.locus[kloc];
//////////////////////////////////////		
		//cat = particuleobs.locuslist[kloc].type % 5;
		//if (debuglevel==2) cout<<"kloc="<<kloc<<"   cat="<<cat<<"\n";
		//particuleobs.locuslist[kloc].name =  new char[strlen(this->dataobs.locus[kloc].name)+1];
		//strcpy(particuleobs.locuslist[kloc].name,this->dataobs.locus[kloc].name);
		//if (debuglevel==2) cout<<"locus "<<kloc<<"   groupe "<<particuleobs.locuslist[kloc].groupe<<"   type="<<this->dataobs.locus[kloc].type<<"\n";
		if (dataobs.locus[kloc].type < 5) {
			kmoy=(dataobs.locus[kloc].maxi+dataobs.locus[kloc].mini)/2;
			particuleobs.locuslist[kloc].kmin=kmoy-((dataobs.locus[kloc].motif_range/2)-1)*dataobs.locus[kloc].motif_size;
			particuleobs.locuslist[kloc].kmax=particuleobs.locuslist[kloc].kmin+(dataobs.locus[kloc].motif_range-1)*dataobs.locus[kloc].motif_size;
			particuleobs.locuslist[kloc].motif_size = dataobs.locus[kloc].motif_size;
			particuleobs.locuslist[kloc].motif_range = dataobs.locus[kloc].motif_range;
			if ((particuleobs.locuslist[kloc].kmin>dataobs.locus[kloc].mini)or(particuleobs.locuslist[kloc].kmax<dataobs.locus[kloc].maxi)) {
				erreur <<"Job aborted : motif range at locus "<<kloc+1<<" is not large enough to include all observed alleles.\n";
				throw std::range_error( erreur.str() ); //exit(1);
			}
		}else if (dataobs.locus[kloc].type < 10) {
			/*cout<<"type du locus = "<<this->dataobs.locus[kloc].type<<"\n";
			particuleobs.locuslist[kloc].dnalength =  this->dataobs.locus[kloc].dnalength;
			particuleobs.locuslist[kloc].pi_A = this->dataobs.locus[kloc].pi_A ;
			particuleobs.locuslist[kloc].pi_C =  this->dataobs.locus[kloc].pi_C;
			particuleobs.locuslist[kloc].pi_G =  this->dataobs.locus[kloc].pi_G;
			particuleobs.locuslist[kloc].pi_T =  this->dataobs.locus[kloc].pi_T;
			cout<<"apres les pi\n";
			particuleobs.locuslist[kloc].haplodna = new string*[particuleobs.data.nsample];
			cout<<"nsample="<<particuleobs.data.nsample<<"\n";
			for (int sa=0;sa<particuleobs.data.nsample;sa++){
				particuleobs.locuslist[kloc].haplodna[sa] = new string[particuleobs.data.ssize[cat][sa]];
				cout<<"sa="<<sa<<"  ss="<<particuleobs.data.ssize[cat][sa]<<"\n";
				for (int i=0;i<particuleobs.data.ssize[cat][sa];i++){
					cout<<"i="<<i<<"   "<<this->dataobs.locus[kloc].haplodna[sa][i]<<"\n";
					particuleobs.locuslist[kloc].haplodna[sa][i] =this->dataobs.locus[kloc].haplodna[sa][i];
					
				}
			}
			cout<<"apres type<10\n";*/
		}else {
			//if (debuglevel==2) cout<<"coucou\n";
			particuleobs.locuslist[kloc].samplesize =  vector<int>(particuleobs.nsample);
			for (int ech=0;ech<particuleobs.nsample;ech++) 
				particuleobs.locuslist[kloc].samplesize[ech] = dataobs.locus[kloc].samplesize[ech];
			//if (debuglevel==2) cout<<"coucou-2\n";
				
			particuleobs.locuslist[kloc].ploidie = vector < vector <short int > >(particuleobs.nsample);
			//if (debuglevel==2) cout<<"coucou-3\n";
			for (int ech=0;ech<particuleobs.nsample;ech++) {
				particuleobs.locuslist[kloc].ploidie[ech] = vector <short int> (dataobs.nind[ech]);
			//if (debuglevel==2) cout<<"coucou-4\n";
				for (int i=0;i<dataobs.nind[ech];i++) {
					particuleobs.locuslist[kloc].ploidie[ech][i] = dataobs.locus[kloc].ploidie[ech][i];
				}
			}
			//if (debuglevel==2) cout<<"coucou-5\n";
			
			
			particuleobs.locuslist[kloc].weight = 1.0;
			particuleobs.locuslist[kloc].haplosnp = vector < vector <short int> >(particuleobs.nsample);
			for (int sa=0;sa<particuleobs.nsample;sa++){
				particuleobs.locuslist[kloc].haplosnp[sa] = vector < short int>(particuleobs.locuslist[kloc].samplesize[sa]);
				int ii=0,imax=(int)dataobs.locus[kloc].haplosnp[sa].size();
				for (int i=0;i<imax;i++){
					if (dataobs.locus[kloc].haplosnp[sa][i] != 9) { 
						particuleobs.locuslist[kloc].haplosnp[sa][ii] =dataobs.locus[kloc].haplosnp[sa][i];
						//if (kloc==0) cout<<"particuleobs.locuslist["<<kloc<<"].haplosnp["<<sa<<"]["<<ii<<"]="<<particuleobs.locuslist[kloc].haplosnp[sa][ii]<<"\n";
						ii++;
					}
				}
			}
		}
	}
	if (debuglevel==2) cout<<"apres la partie locus\n";
	//string *sb,
	string ent;
	if (debuglevel==2) cout<<"this->entetestat = "<<this->entetestat<<"\n";
	//sb = splitwords(this->entete," ",&j);
	//cout<<"this->entete = "<<this->entete<<"\nn="<<j<<"   this->nstat="<<this->nstat<<"\n";
	ent=entetestat+"\n";
	//for (int k=j-this->nstat;k<j;k++) ent=ent+centre(sb[k],14);
	//ent=ent+"\n";
	//if (debuglevel==2)cout<<"entete : "<<entete<<"\n";
	//delete []sb;
	jstat=0;
	if (debuglevel==2) cout<<"avant FILE *fobs\n";
	ofstream fobs;
	fobs.open(statobsfilename.c_str());
	fobs<<ent<<"\n";
	fobs<<setiosflags(ios::fixed)<<setprecision(8);
	ast.resize(particuleobs.ngr+1);
	for(int gr=1;gr<particuleobs.ngr;gr++) {
		if (debuglevel==2)cout<<"avant calcul des statobs du groupe "<<gr<<"\n";
		particuleobs.docalstat(gr);
		jstat +=particuleobs.grouplist[gr].nstat;
		if (debuglevel==2)cout<<"apres calcul des statobs du groupe "<<gr<<"\n";
		ast[gr].resize(0);
		for (int j=0;j<particuleobs.grouplist[gr].nstat;j++){ 
			//cout<<"particuleobs.grouplist["<<gr<<"].sumstat["<<j<<"].val="<<particuleobs.grouplist[gr].sumstat[j].val<<"\n";
			if (particuleobs.grouplist[gr].sumstat[j].val!=-9999.0){
				fobs<<setw(12)<<particuleobs.grouplist[gr].sumstat[j].val<<"  ";
			} else {
				ast[gr].push_back(j);
			};
		}
	}
	fobs<<"\n";
	fobs.close();
	int nast=0;for(int gr=1;gr<particuleobs.ngr;gr++) nast+=ast[gr].size();
	//cout<<"nast="<<nast<<"\n";
	if (nast>0) {
		string message;
		message = "The following introgression rate summary statistics : ";
		for(int gr=1;gr<particuleobs.ngr;gr++) {
			if (ast[gr].size()>0) {
				message = message + "\nGroup "+IntToString(gr)+" : ";
				int jmax=(int)ast[gr].size();
				for (int j=0;j<jmax;j++) {
					message += IntToString(particuleobs.grouplist[gr].sumstat[ast[gr][j]].samp)+"&";
					message += IntToString(particuleobs.grouplist[gr].sumstat[ast[gr][j]].samp1)+"&";
					message += IntToString(particuleobs.grouplist[gr].sumstat[ast[gr][j]].samp2)+"  ";
				}
			}
		}
		if (nast==1)  message += "\ncannot be estimated in the observed data set. Please don't use it.";
				else  message += "\ncannot be estimated in the observed data set. Please don't use them.";
		return message;
	}
	
	this->stat_obs = vector<float>(this->nstat);
	jstat=0;
	for(int gr=1;gr<particuleobs.ngr;gr++) {
		for (int j=0;j<particuleobs.grouplist[gr].nstat;j++) {
			this->stat_obs[jstat] = (float)particuleobs.grouplist[gr].sumstat[j].val;
			//cout<<"stat_obs["<<j<<"]="<<this->stat_obs[j]<<"\n";
			if (this->stat_obs[jstat]!=-9999.0) jstat++;
		}
	}
	//cout<<"header.nstat = "<<this->nstat<<"     jstat="<<jstat<<"\n";
	if (this->nstat<11) {
		for (int j=0;j<this->nstat;j++) cout<<"stat_obs["<<j<<"]="<<this->stat_obs[j]<<"\n";
	} else {
		for (int j=0;j<5;j++)  cout<<"stat_obs["<<j<<"]="<<this->stat_obs[j]<<"\n";
		cout<<"   .../...\n";
		for (int j=this->nstat-6;j<this->nstat;j++)  cout<<"stat_obs["<<j<<"]="<<this->stat_obs[j]<<"\n";
	}
	if (debuglevel==2) cout<<"avant libere \n";
	//particuleobs.libereobs(true);
	if (debuglevel==2) cout<<"fin de calstatobs\n";
	return "OK";
}

/**
* lit le fichier des statistiques observées (placées dans double *stat_obs)
*/
/*    double* read_statobs(char *statobsfilename) {
       string entete,ligne,*ss;
       int ns;
       ifstream file(statobsfilename, ios::in);
       getline(file,entete);
       getline(file,ligne);
       file.close();
       ss=splitwords(ligne," ",&ns);
       //cout<<"statobs ns="<<ns<<"\n";
       if (ns!=this->nstat) exit(1);
       this->stat_obs = new double[ns];
       for (int i=0;i<ns;i++) this->stat_obs[i]=atof(ss[i].c_str());
       delete []ss;
       return this->stat_obs;
       //for (int i=0;i<ns;i++) cout<<stat_obs[i]<<"   ";
       //cout<<"\n";
    }*/



