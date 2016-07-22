/*
 * acploc.h
 *
 *  Created on: 13 déc. 2011
 *      Author: ppudlo
 */

#pragma once

#include <string>
#include <vector>

using namespace std;

struct resACPC {
	int index, nlambda; //index=0 pour matrice de corrélation et !=0 pour matrice de covariance
	long double proportion, *lambda, slambda, **vectprop, **princomp, *moy, *sd;


};

resACPC ACP(int nli, int nco, vector<vector<long double>>& X, long double prop, int index);
void cal_acp();
void cal_loc();
void doacpl(string opt);
