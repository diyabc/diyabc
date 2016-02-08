/*
 * particleset.h
 *
 *  Created on: 9 déc. 2011
 *      Author: ppudlo
 */

#ifndef PARTICLESET_H_
#define PARTICLESET_H_

#include "particuleC\.hpp"


/**
* Structure ParticleSet : définition d'un ensemble de particules
*/
struct ParticleSetC
{
	std::vector<ParticleC> particule;
	std::vector<LocusC> locuslist;
	//std::vector<LocusGroupC> grouplist;
	int npart,nloc,ngr;
	int nsample,nscenarios;
	std::vector<int> nind;
	std::vector<std::vector<int> > indivsexe;
	//std::vector<ScenarioC> scenario;
	double sexratio;
    MwcGen *mw;

    // remplissage de la partie data de la particule p
    //void setdata(int p);
	void init_particule(int npart, bool dnatrue, bool simulfile, int seed);
	// remplissage de la partie groupe de la particule p
    void setgroup(int p);
    // remplissage de la partie locus de la particule p
    void setloci(int p);
    // remplissage de la partie scenarios de la particule p
    void setscenarios (int p, bool simulfile);
    // réinitialisation de la particule p
    void resetparticle (int p);
	// liberation des particules
	//void libere(int npart);
	// simulation de particules de paramètres fixées
	void dosimulstat(int debut, int npart, bool dnatrue, bool multithread, int numscen, int seed, double **stat);    
	// simulation de particules de paramètres fixées
	void dosimulstat2(int debut, int npartic, vector < vector <long double> > &phistarOKK, bool dnatrue, bool multithread, int numscen, int seed, double** stat);    
	// simulation des particules utilisées pour le model checking
    void dosimulphistar(int npart, bool dnatrue,
    		            bool multithread,bool firsttime, int numscen,int seed,int nsel);
    // simulation des particules utilisées pour la table de référence, le biais et la confiance
    void dosimultabref(int npart, bool dnatrue,bool multithread,
    		           bool firsttime, int numscen,int seed, int depuis);
    // simulation des particules utilisées pour la création de fichiers genepop
    string* simulgenepop(int npart, bool multithread, int seed);
    // simulation des particules utilisées pour la création de fichiers SNP
    string* simuldataSNP(int npart, bool multithread, int seed);
	void free_particule(int npart);
};


#endif /* PARTICLESET_H_ */
