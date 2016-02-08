/*
 * history.h
 *
 *  Created on: 9 déc. 2011
 *      Author: ppudlo
 */

#ifndef HISTORY_H_
#define HISTORY_H_

#include <string>
#include "randomgenerator\.hpp"


/**
 * struct StatC :éléments de définition d'une summary statistic
 */
struct StatC
{
  int cat,samp,samp1,samp2,group,numsnp;
  long double val;
//  StatC(StatC const & source);
  StatC & operator= (StatC const & source);
};

/**
 * Classe StatsnpC :éléments de définition d'une summary statistic pour les snp
 */
class StatsnpC
{
public:
  int cat,samp,samp1,samp2,group,n;
  std::vector <long double> x,w;
  long double sw;
  bool defined;
  StatsnpC & operator= (StatsnpC const & source);
  ~StatsnpC(){
    if(not x.empty()) x.clear();
    if(not w.empty()) w.clear();
  };
};





/**
 * Classe PriorC :éléments de définition d'un prior de paramètre historique
 */
class PriorC
{
public:
  std::string loi;
  double mini,maxi,mean,sdshape;
  int ndec;
  bool constant,fixed;
  PriorC() {
	  loi = "";
}
//  PriorC(PriorC const & source);
  PriorC & operator= (PriorC const & source);
  void ecris();
  double drawfromprior(MwcGen & mw);
  void readprior(std::string ss);
};

/**
 * Classe ConditionC :éléments de définition d'une condition sur un couple de paramètres historiques
 */

class ConditionC
{
public:
  std::string param1,param2,operateur;
	ConditionC(){
		param1 = "";
		param2 = "";
		operateur = "";
	};
  ConditionC & operator= (ConditionC const & source);
  void ecris();
  void readcondition(std::string ss);

};

/**
 * Classe EventC :éléments de définition d'un événement populationnel
 */
class EventC
{
public:
  char action;   //"V"=VarNe "M"=Merge  "S"=Split  "E" = sample/echantillon
  int pop,pop1,pop2,sample,Ne,time;
  double admixrate;
  int numevent0,nindMref,nindFref;
  std::string stime, sNe, sadmixrate;
  // char *stime,*sNe,*sadmixrate;
  // int ltime,lNe,ladmixrate;


  /**
   * Déclarations des méthodes (PP)
   */

  // Définition de l'opérateur < pour le tri des événements populationnels (PP)
  friend bool operator< (const EventC & lhs, const EventC & rhs);
	
	EventC & operator= (EventC const & source);
  // Pour écrire l'événement dans la sortie standard
  void ecris();
};

/**
 * struct Ne0C : éléments de définition d'un effectif efficace à la ligne 1 d'un scénario
 */
class Ne0C
{
public:
  int val;
  std::string name;
  
  Ne0C & operator= (Ne0C const & source);
};

/**
 * Classe HistparameterC :éléments de définition d'un paramètre historique
 */
class HistParameterC
{
public:
  std::string name;
  int category;   //0 pour N, 1 pour T et 2 pour A
  double value;
  PriorC prior;
  HistParameterC(){
	  name = "";
}
//  HistParameterC(HistParameterC const & source);
  HistParameterC & operator= (HistParameterC const & source);
  void ecris(bool simulfile);
};


/**
 * struct LocusGroupC : éléments de définition d'un groupe de locus
 */
class LocusGroupC
{
  public:
    std::vector<int> loc;     // loc=numeros des locus du groupe
    int nloc,nstat,nstatsnp;
    int type;                      //O= microsat, 1=sequence
    double p_fixe,gams;
    double musmoy,mutmoy,Pmoy,snimoy;
    double k1moy,k2moy;
    int  mutmod;
    PriorC priormusmoy,priork1moy,priork2moy,priormusloc,priork1loc,priork2loc;
    PriorC priormutmoy,priorPmoy,priorsnimoy,priormutloc,priorPloc,priorsniloc;
    std::vector<StatC> sumstat;
    std::vector<StatsnpC> sumstatsnp;
  /* Déclaration des méthodes */
//  LocusGroupC(LocusGroupC const & source);
  LocusGroupC & operator= (LocusGroupC  const & source);

  ~LocusGroupC(){
    if (not loc.empty()) loc.clear();
    if (not sumstat.empty()) sumstat.clear();
    if (not sumstatsnp.empty()) sumstatsnp.clear();
  };


};

/**
 * Classe ScenarioC :éléments de définition d'un scénario
 */
class ScenarioC
{
public:
  std::vector<double> paramvar;
  double prior_proba;
  int number,popmax,npop,nsamp,nparam,nevent,nn0,nparamvar,nconditions,ipv;
  std::vector<int> time_sample;
  std::vector<std::string> stime_sample;
  std::vector<EventC> event;
  std::vector<Ne0C> ne0;
  std::vector<HistParameterC> histparam;
  std::vector<ConditionC> condition;
  /* Déclaration des méthodes */

 // ScenarioC(ScenarioC const & source);
  ScenarioC & operator= (ScenarioC  const & source);
  ~ScenarioC(){
  	if (not paramvar.empty()) paramvar.clear();
  	if (not time_sample.empty()) time_sample.clear();
  	if (not stime_sample.empty()) stime_sample.clear();
  	if (not event.empty()) event.clear();
  	if (not ne0.empty()) ne0.clear();
  	if (not histparam.empty()) histparam.clear();
  	if (not condition.empty()) condition.clear();
  };
  
  /* détermination du ou des paramètres contenus dans la std::string s */
  void detparam(std::string s,int cat);
  /* lecture/interprétation des lignes d'un scénario */
  std::string read_events(int nl,std::string *ls);
  /* verification d'un scénario*/
  std::string checklogic();
  void ecris(bool simulfile);
  void libere();
}; /* fin classe ScenarioC */



/**
 * Struct SequenceBitC : éléments de définition d'un segment (vertical) de l'arbre de coalescence'
 */
struct SequenceBitC
{
  /* action = "C" (coal), "M" (merge), "S" (split), "A" (adsamp)
   */
  char action;
  int pop,pop1,pop2,sample;
  int N,t0,t1;
  double admixrate;
  std::vector<bool> popfull;

  ~SequenceBitC(){
	  //std::cout<<"passage dans le destructeur de SequenceBitC\n";
	  if (not popfull.empty()) popfull.clear();
	  //std::cout<<"   sortie du destructeur de SequenceBitC\n";
}
  void ecris();
};



/**
 * Struct NodeC : éléments de définition d'un noeud de l'arbre de coalescence
 */
struct NodeC
{
  int pop,sample,state,brhaut,ndat,nref;
  double height;
  std::string dna;
  bool OK;
  ~NodeC(){}
};

/**
 * Struct BranchC : éléments de définition d'une branche de l'arbre de coalescence
 */
struct BranchC
{
  int bottom,top,nmut;
  double length;
  bool OK,OKOK;
  ~BranchC(){}
};

/**
 * Class GeneTreeC : éléments de définition d'un arbre de coalescence
 */
class GeneTreeC
{
public:
  std::vector <NodeC> nodes;
  std::vector <BranchC> branches;
  int nmutot,nnodes,nbranches,ngenes,nbOK,nbOKOK;

  /* Déclaration des méthodes */
  ~GeneTreeC(){
	//std::cout<<"passage dans le destructeur de GeneTreeC\n";
    if(not nodes.empty()) nodes.clear();
    if(not branches.empty()) branches.clear();
  };

  GeneTreeC & operator=(GeneTreeC const & source);
  void ecris();

};



#endif /* HISTORY_H_ */
