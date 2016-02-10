/*
 * conf.h
 *
 *  Created on: 13 déc. 2011
 *      Author: ppudlo
 */

#pragma once

#include <string>

// Ecriture de l'entete du fichier confidence.txt contenant les résultats
void ecrientete(int nrec, int ntest, int nseld, int nselr, int nlogreg, std::string shist, std::string smut, bool AFD, int nsel0);

// Appelle l'analyse factorielle discriminante et remplace les summary stats par les composantes des jeux simulés
// sur les axes discriminants
float* transfAFD(int nrec, int nsel, int p);

// Interprête la ligne de paramètres de l'option "confiance dans le choix du scénario"
// et lance les calculs correspondants
void doconf(std::string opt, int seed);

