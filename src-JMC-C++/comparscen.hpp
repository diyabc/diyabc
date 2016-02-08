/*
 * comparscen.h
 *
 *  Created on: 13 déc. 2011
 *      Author: ppudlo
 */

#ifndef COMPARSCEN_H_
#define COMPARSCEN_H_
#include <string>
#include "header.hpp"
#include "matrices.hpp"
#include "mesutils.hpp"


struct matligneC
{
    long double* x;
};

/**
 *  Structure contenant les résultats de l'estimation de la probabilité a posteriori d'un scénario
 */
struct posteriorscenC
{
    long double x,inf,sup;
	int err;
};
/**
* définit l'opérateur de comparaison de deux lignes de la matrice cmatA
* pour l'utilisation de la fonction sort du module algorithm
*/
struct complignes
{
   bool operator() (const matligneC & lhs, const matligneC & rhs) const
   {
      return lhs.x[0] < rhs.x[0];
   }
};

void allouecmat(int nmodel, int nli, int nco);
void liberecmat(int nmodel, int nli, int nco);
posteriorscenC** comp_direct(int n);
void  save_comp_direct(int n, posteriorscenC** posts, std::string path, std::string ident);
void save_comp_logistic(int nlogreg,int nselr,posteriorscenC** postscenlog,std::string path, std::string ident);
void rempli_mat0(int n, float* stat_obs);
void expbeta(int bsize, long double *b, long double *eb);
void remplimatriceYP(int nli, int nco, int nmodel, long double **cmatP,
		long double *cmatYP, long double *cbeta, long double **cmatX,
		long double *cvecW, long double **cmatY, long double *csmatP);
bool cal_loglik(int nli, int nmodel, int rep, long double *cloglik,
		long double **cmatY, long double **cmatP, long double *vecW,
		long double *csmatY, long double *csmatP);
void calcul_psd(int nmodel, long double *b0, long double **matV, long double *sd, long double *px);
void ordonne(int nmodel,int nli,int nco,long double *vecY,int *numod);
void reordonne(int nmodel,int *numod,long double *px, long double *pxi,long double *pxs);
int polytom_logistic_regression(int nli, int nco, long double **cmatX0,
		long double *vecY, long double *cvecW, long double *px, long double *pxi, long double *pxs);
posteriorscenC* call_polytom_logistic_regression(int nts, float *stat_obs, int nscenutil,int *scenchoisiutil);
posteriorscenC* comp_logistic(int nts,float *stat_obs);
void docompscen(std::string opt);

#endif /* COMPARSCEN_H_ */
