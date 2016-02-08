/*
 * estimparam.h
 *
 *  Created on: 13 déc. 2011
 *      Author: ppudlo
 */

#ifndef ESTIMPARAM_H_
#define ESTIMPARAM_H_
#include <string>
#include "history.hpp"

/**
* définition de la densité d'un paramètre
*/
struct pardensC {
    long double *x,*priord,*postd;
    long double xmin,xmax,xdelta;
    int ncl;
};

struct parstatC {
    double moy,med,mod,q025,q050,q250,q750,q950,q975;
};


/**
* compte les paramètres communs aux scénarios choisis (nparamcom), le nombre de paramètres composites
* remplit le tableau numpar des numéros de paramètres communs par scénario
*/
void det_numpar();
/**
* recalcule les valeurs de paramètres en fonction de la transformation choisie
*/
void recalparamO(int n);
void recalparamC(int n);
void recalparamS(int n);
/**
* effectue le remplissage de la matrice matX0, du vecteur des poids vecW et
* de la matrice des paramètres parsim (éventuellement transformés)
*/
void rempli_mat(int n, float* stat_obs);
void delete_mat(int n);
void rempli_parsim(int n, int npa);
int ecrimatL(std::string nomat, int n, int m, long double **A);
int ecrimatD(std::string nomat, int n, int m, double **A);
/**
* effectue la régression locale à partir de la matrice matX0 et le vecteur des poids vecW
*/
long double** cal_matC(int n);
void local_regression(int n, int npa, long double **matC);
/**D
* calcule les phistars pour les paramètres originaux et composites
*/
void calphistarO(int n, long double **phistar);
void calphistarC(int n, long double **phistarcompo);
void calphistarS(int n, long double **phistarscaled);
void copyphistar(int n, int nparam, long double **phistar);
void det_nomparam();
/**
* effectue la sauvegarde des phistars dans le fichier path/ident/phistar.txt
*/
void savephistar(int n, std::string path, std::string ident, long double **phistar, long double **phistarcompo, long double **phistarscaled) ;
/**
* lit les paramètres des enregistrements simulés pour l'établissement des distributions a priori'
*/
void lisimparO();
void lisimparC();
void lisimparS();
/**
* calcule la densité à partir de la loi fournie dans le prior
*/
long double* caldensexact(int ncl,long double *x,PriorC pr);
/**
* calcule un histogramme à partir de la loi fournie dans le prior
*/
long double* calhistexact(int ncl,long double *x,PriorC pr);
/**
* calcule la densité à partir d'un échantillon de valeurs'
*/
long double* calculdensite(int ncl,int n, long double *x, long double **y,int j,
		bool multithread);
/**
* calcule un histogramme à partir d'un échantillon de valeurs'
*/
long double* calculhisto(int ncl, int n, long double *x, long double **y,int j,
		bool multithread);
/**
* traitement global du calcul de la densité des paramètres variables
* si le parametre est à valeurs entières avec moins de 30 classes, la densité est remplacée par un histogramme
* sinon la densité est évaluée pour 1000 points du min au max
*/
void histodens(int n, bool multithread, std::string progressfilename,int* iprog,int* nprog);
/**
*calcule les statistiques des paramètres
*/
parstatC* calparstat(int n, int nparam, long double **par);
//parstatC* calparstatO(int n, long double **phistar);
//parstatC* calparstatC(int n, long double **phistarcompo);
//parstatC* calparstatS(int n, long double **phistarscaled);
/**
*sauvegarde les statistiques et les densités des paramètres
*/
void saveparstat(int nsel, std::string path, std::string ident);
/**
* Interprète les paramètres de la ligne de commande et
* effectue l'estimation ABC des paramètres (directe + régression locale)
*/
void doestim(std::string opt);

#endif /* ESTIMPARAM_H_ */
