/*
 * bias.h
 *
 *  Created on: 13 déc. 2011
 *      Author: ppudlo
 */

#ifndef BIAS_H_
#define BIAS_H_
#include <vector>
#include <string>

#include "header.hpp"
#include "matrices.hpp"
#include "mesutils.hpp"
#include "particleset.hpp"


struct enreC {
    int numscen;
    double *paramvv,*paramvvC,*paramvvS,*stat;
	string *name;
};

// modifie les paramètres historiques en accord avec les priors/valeurs des pseudo-observés
void resethistparam(std::string s);
// modifie les conditions sur les paramètres historiques en accord avec les priors des pseudo-observés
void resetcondition(int j,std::string s);
// modifie les paramètres mutationnels en accord avec les priors des pseudo-observés
void resetmutparam(std::string s);
// calcule les différentes statistiques de biais, rmse...
void biaisrel(int ntest,int nsel,int npv);
// Mise en forme des résultats
void ecrires(int ntest,int npv,int nsel);
//Traitement des paramètres composites (microsat et séquences uniquement)
void setcompo(int p);
//Interprête la ligne de paramètres de l'option biais et lance les calculs correspondants
void dobias(std::string opt,  int seed);


#endif /* BIAS_H_ */
