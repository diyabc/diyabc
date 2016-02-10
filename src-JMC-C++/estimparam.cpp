/*
 * estimparam.cpp
 *
 *  Created on: 9 march 2011
 *      Author: cornuet
 */

#include <string>
#include <vector>
#include <algorithm>
#include <cmath>
#include <iomanip>
#include <time.h>

#include "header.hpp"
#include "reftable.hpp"
#include "matrices.hpp"
#include "mesutils.hpp"
#include "estimparam.hpp"

using namespace std;

extern string progressfilename;
extern ReftableC rt;
extern int debuglevel;
extern vector<LocusGroupC> groupe;
extern vector<ScenarioC> scenario;

extern ofstream fprog;

#define c 1.5707963267948966192313216916398
#define co 2.506628274631000502415765284811   //sqrt(pi)
#define coefbw 1.8                            //coefficient multiplicateur de la bandwidth pour le calcul de la densité


extern ReftableC rt;
extern HeaderC header;
long double **matX0, *vecW, **alpsimrat;
long double** parsim;
int nstatOKsel;
//extern char *headerfilename, *reftablefilename,*datafilename,*statobsfilename, *reftablelogfilename,*path,*ident;
extern string headerfilename, reftablefilename, datafilename, statobsfilename, reftablelogfilename, path, ident;
extern bool multithread;
pardensC *pardens, *pardenscompo, *pardensscaled;
long double *var_stat, *parmin, *parmax, *parmincompo, *parmaxcompo, *parminscaled, *parmaxscaled, *diff, *diffcompo, *diffscaled;
long double *parmin0, *parmax0, *parmincompo0, *parmaxcompo0, *parminscaled0, *parmaxscaled0;
long double **beta, **simpar, **simparcompo, **simparscaled;
int nparamcom, nparcompo, nparscaled, **numpar, numtransf = 3, npar, npar0;
long double borne = 10.0, xborne, xbornecompo, xbornescaled;
bool composite = false, scaled = false, original = true;
string enteteO, enteteC, enteteS;
int nsimpar = 100000;
string *nomparamO, *nomparamC, *nomparamS;
parstatC *parstat, *parstatcompo, *parstatscaled;
bool deltanul;


/**
 * compte les paramètres communs aux scénarios choisis (nparamcom), le nombre de paramètres composites
 * remplit le tableau numpar des numéros de paramètres communs par scénario
 */
void det_numpar()
{
    cout << "debut de det_numpar\n";
    vector<string> parname;
    int ii, iscen;
    bool commun, trouve;
    npar = 0;
    npar0 = 0;
    numpar = new int*[rt.nscenchoisi];
    parname.resize(0);
    if (rt.nscenchoisi == 1)
    {
        iscen = rt.scenchoisi[0] - 1;
        //cout<<"rt.nparam[rt.scenchoisi[0]-1] = "<<rt.nparam[iscen]<<"\n";
        numpar[0] = new int[rt.nparam[iscen]];
        //for (int i=0;i<rt.nparam[iscen];i++) {numpar[0][i]=0;cout<<"numpar[0]["<<i<<"]="<<numpar[0][i]<<"\n";}
        for (int i = 0; i < rt.nhistparam[iscen]; i++)
        {
            numpar[0][i] = i;
            npar++;
            parname.push_back(rt.histparam[iscen][i].name);
            if (rt.histparam[iscen][i].category < 2) npar0++;
        }
        //cout<<"avant nparamcom=\n";
        cout << "npar=" << npar << "\n";
        cout << "rt.nparamut=" << rt.nparamut << "\n";
        nparamcom = npar + rt.nparamut;//rt.nparam[rt.scenchoisi[0]-1]-scenario[rt.scenchoisi[0]-1].nparam;
        cout << "nparamcom=" << nparamcom << "\n";
        cout << "rt.nparam[iscen] = " << rt.nparam[iscen] << "\n";
        //for (int i=0;i<rt.nparam[iscen];i++) {cout<<"numpar[0]["<<i<<"]="<<numpar[0][i]<<"\n";}
        for (int i = rt.nhistparam[iscen]; i < nparamcom; i++)
        {
            numpar[0][i] = i;/*cout<<"numpar[0]["<<i<<"]="<<numpar[0][i]<<"\n";*/
        }
        cout << "noms des parametres communs : ";
        for (int i = 0; i < rt.nhistparam[iscen]; i++) cout << rt.histparam[iscen][i].name << "   ";
        for (int i = 0; i < rt.nparamut; i++) cout << rt.mutparam[i].name << "   ";
        cout << "\n";
    }
    else
    {
        iscen = rt.scenchoisi[0] - 1; //on prend les paramètres du premier scénario choisi et on les compare à ceux des autres scénarios
        for (int i = 0; i < rt.nhistparam[iscen]; i++)
        {
            commun = true;
            cout << "parametre teste : " << rt.histparam[rt.scenchoisi[0] - 1][i].name << "\n";
            for (int j = 1; j < rt.nscenchoisi; j++)
            {
                trouve = false;
                cout << "rt.nhistparam[rt.scenchoisi[" << j << "]-1] = " << rt.nhistparam[rt.scenchoisi[j] - 1] << "\n";
                for (int k = 0; k < rt.nhistparam[rt.scenchoisi[j] - 1]; k++)
                {
                    cout << "   " << rt.histparam[rt.scenchoisi[j] - 1][k].name;
                    trouve = (rt.histparam[rt.scenchoisi[0] - 1][i].name.compare(rt.histparam[rt.scenchoisi[j] - 1][k].name) == 0);
                    if (trouve)
                    {/*cout<< "\nparam commun : "<<rt.histparam[rt.scenchoisi[0]-1][i].name<<"\n";*/
                        break;
                    }
                }
                commun = trouve;
                if (not commun) break;
            }
            if (commun)
            {
                if (not rt.histparam[rt.scenchoisi[0] - 1][i].prior.constant)
                {
                    npar++;//cout<<"npar = "<<npar<<"\n";
                    parname.push_back(rt.histparam[rt.scenchoisi[0] - 1][i].name);
                    //cout<<parname[npar-1]<<"\n\n";
                    if (rt.histparam[rt.scenchoisi[0] - 1][i].category < 2) npar0++;
                }
            }
            else cout << "\nnon commun\n\n";
        }
        cout << "avant affichage des parname\n";
        for (int i = 0; i < npar; i++) cout << parname[i] << "   ";
        cout << "\n";
        nparamcom = npar + rt.nparamut;//rt.nparam[rt.scenchoisi[0]-1]-scenario[rt.scenchoisi[0]-1].nparam;
        for (int j = 0; j < rt.nscenchoisi; j++)
        {
            numpar[j] = new int[nparamcom];
            ii = 0;
            for (int i = 0; i < npar; i++)
            {
                for (int k = 0; k < rt.nhistparam[rt.scenchoisi[j] - 1]; k++)
                {
                    if (rt.histparam[rt.scenchoisi[j] - 1][k].name.compare(parname[i]) == 0)
                    {
                        numpar[j][ii] = k;
                        ii++;
                    }
                }
            }
            for (int k = rt.nhistparam[rt.scenchoisi[j] - 1]; k < rt.nparam[rt.scenchoisi[j] - 1]; k++)
            {
                numpar[j][ii] = k;
                ii++;
            }
            /*for (int kk=0;kk<npar+header.nparamut;kk++) {
              cout <<numpar[j][kk]<<"  ";
              if (kk<npar) cout<<"("<<rt.histparam[rt.scenchoisi[j]-1][numpar[j][kk]].name<<")  ";
              else cout<<"("<<rt.mutparam[kk-npar].name<<")  ";
              }
              cout<<"\n";*/
        }
        /*cout << "noms des parametres communs : ";
          for (int i=0;i<npar;i++) cout<<parname[i]<<"   ";
          cout <<"\n";*/
    }
    /*cout << "noms des parametres communs : ";
      for (int j=0;j<rt.nscenchoisi;j++) {
      cout<<"scenario "<<rt.scenchoisi[j]<<"\n";
      for (int i=0;i<nparamcom;i++) {cout<<numpar[j][i]<<" ";
      if (i<npar) cout<<"("<<rt.histparam[rt.scenchoisi[j]-1][numpar[j][i]].name<<")  ";
      else cout<<"("<<rt.mutparam[i-npar].name<<")  ";
      }
      cout <<"\n";
      }*/
    if (composite) nparcompo = npar0 * header.ngroupes; else nparcompo = 0;
    if (scaled) nparscaled = npar0; else nparscaled = 0;
}

/**
 * recalcule les valeurs de paramètres en fonction de la transformation choisie
 * pour les paramètres originaux
 */
void recalparamO(int n)
{
    //long double coefmarge=0.001,marge;
    int jj, k;
    alpsimrat = new long double*[n];
    for (int i = 0; i < n; i++) alpsimrat[i] = new long double[nparamcom];
    if (deltanul)
    {
        for (int i = 0; i < n; i++)
        {
            for (int j = 0; j < nparamcom; j++)
            {
                k = 0;
                while (rt.enrsel[i].numscen != rt.scenchoisi[k])k++;
                alpsimrat[i][j] = rt.enrsel[i].param[numpar[k][j]];
            }
        }
        return;
    }
    switch (numtransf)
    {
    case 1: //no transform
        for (int i = 0; i < n; i++)
        {
            for (int j = 0; j < nparamcom; j++)
            {
                k = 0;
                while (rt.enrsel[i].numscen != rt.scenchoisi[k])k++;
                alpsimrat[i][j] = rt.enrsel[i].param[numpar[k][j]];
            }
        }
        break;
    case 2: // log transform
        for (int i = 0; i < n; i++)
        {
            for (int j = 0; j < nparamcom; j++)
            {
                k = 0;
                while (rt.enrsel[i].numscen != rt.scenchoisi[k])k++;
                if (rt.enrsel[i].param[numpar[k][j]] <= 0.00000000001) rt.enrsel[i].param[numpar[k][j]] = 1E-11;
                alpsimrat[i][j] = log(rt.enrsel[i].param[numpar[k][j]]);
            }
        }
        break;
    case 3: //logit transform
        //cout<<"borne="<<borne<<"\n";
        parmin = new long double[nparamcom];
        parmax = new long double[nparamcom];
        diff = new long double[nparamcom];
        parmin0 = new long double[nparamcom];
        parmax0 = new long double[nparamcom];
        xborne = borne;
        for (int j = 0; j < nparamcom - rt.nparamut; j++)
        {
            parmin0[j] = rt.histparam[rt.scenchoisi[0] - 1][numpar[0][j]].prior.mini;
            parmax0[j] = rt.histparam[rt.scenchoisi[0] - 1][numpar[0][j]].prior.maxi;
        }
        for (int j = nparamcom - rt.nparamut; j < nparamcom; j++)
        {
            jj = j + rt.nparamut - nparamcom;
            parmin0[j] = rt.mutparam[jj].prior.mini;
            parmax0[j] = rt.mutparam[jj].prior.maxi;
        }
        for (int j = 0; j < nparamcom - rt.nparamut; j++)
        {
            //cout<<"j="<<j<<"   nparamcom="<<nparamcom<<"\n";
            if (rt.histparam[rt.scenchoisi[0] - 1][numpar[0][j]].category < 2)
            {
                parmin[j] = parmin0[j] - 0.5;
                parmax[j] = parmax0[j] + 0.5;
            }
            else
            {
                parmin[j] = parmin0[j] - 0.0005;
                parmax[j] = parmax0[j] + 0.0005;
            }
            diff[j] = parmax[j] - parmin[j];
        }
        for (int j = nparamcom - rt.nparamut; j < nparamcom; j++)
        {
            jj = j + rt.nparamut - nparamcom;
            parmin[j] = 0.95 * parmin0[j];
            parmax[j] = 1.05 * parmax0[j];
            diff[j] = parmax[j] - parmin[j];
        }
        for (int i = 0; i < n; i++)
        {
            for (int j = 0; j < nparamcom; j++)
            {
                k = 0;
                while (rt.enrsel[i].numscen != rt.scenchoisi[k])k++;
                if (rt.enrsel[i].param[numpar[k][j]] <= parmin0[j]) alpsimrat[i][j] = -xborne;
                else if (rt.enrsel[i].param[numpar[k][j]] >= parmax0[j]) alpsimrat[i][j] = xborne;
                else alpsimrat[i][j] = log((rt.enrsel[i].param[numpar[k][j]] - parmin[j]) / (parmax[j] - rt.enrsel[i].param[numpar[k][j]]));
            }
        }
        break;
    case 4: //log(tg) transform
        parmin = new long double[nparamcom];
        parmax = new long double[nparamcom];
        diff = new long double[nparamcom];
        parmin0 = new long double[nparamcom];
        parmax0 = new long double[nparamcom];
        xborne = borne;
        for (int j = 0; j < nparamcom - rt.nparamut; j++)
        {
            parmin0[j] = rt.histparam[rt.scenchoisi[0] - 1][numpar[0][j]].prior.mini;
            parmax0[j] = rt.histparam[rt.scenchoisi[0] - 1][numpar[0][j]].prior.maxi;
        }
        for (int j = nparamcom - rt.nparamut; j < nparamcom; j++)
        {
            jj = j + rt.nparamut - nparamcom;
            parmin0[j] = rt.mutparam[jj].prior.mini;
            parmax0[j] = rt.mutparam[jj].prior.maxi;
        }
        for (int j = 0; j < nparamcom - rt.nparamut; j++)
        {
            //cout<<"j="<<j<<"   nparamcom="<<nparamcom<<"\n";
            if (rt.histparam[rt.scenchoisi[0] - 1][numpar[0][j]].category < 2)
            {
                parmin[j] = parmin0[j] - 0.5;
                parmax[j] = parmax0[j] + 0.5;
            }
            else
            {
                parmin[j] = parmin0[j] - 0.0005;
                parmax[j] = parmax0[j] + 0.0005;
            }
            diff[j] = parmax[j] - parmin[j];
        }
        for (int j = nparamcom - rt.nparamut; j < nparamcom; j++)
        {
            jj = j + rt.nparamut - nparamcom;
            parmin[j] = 0.95 * parmin0[j];
            parmax[j] = 1.05 * parmax0[j];
            diff[j] = parmax[j] - parmin[j];
        }
        for (int i = 0; i < n; i++)
        {
            for (int j = 0; j < nparamcom; j++)
            {
                k = 0;
                while (rt.enrsel[i].numscen != rt.scenchoisi[k])k++;
                if (rt.enrsel[i].param[numpar[k][j]] <= parmin0[j]) alpsimrat[i][j] = -xborne;
                else if (rt.enrsel[i].param[numpar[k][j]] >= parmax0[j]) alpsimrat[i][j] = xborne;
                else alpsimrat[i][j] = log(tan((rt.enrsel[i].param[numpar[k][j]] - parmin[j]) * c / diff[j]));
            }
        }
        break;
    }
}

/**
 * recalcule les valeurs de paramètres en fonction de la transformation choisie
 * pour les paramètres composites
 */
void recalparamC(int n)
{
    //long double coefmarge=0.000000001,marge;
    long double **xx, pmut;
    int kk, qq, kscen;
    alpsimrat = new long double*[n];
    for (int i = 0; i < n; i++) alpsimrat[i] = new long double[nparcompo];
    xx = new long double*[n];
    for (int i = 0; i < n; i++) xx[i] = new long double[nparcompo];
    int kp = 0, kp0;
    if (numtransf > 2)
    {
        parmincompo = new long double[nparcompo];
        parmaxcompo = new long double[nparcompo];
        diffcompo = new long double[nparcompo];
        parmincompo0 = new long double[nparcompo];
        parmaxcompo0 = new long double[nparcompo];
        for (int j = 0; j < nparcompo; j++)
        {
            parmincompo0[j] = 1E100;
            parmaxcompo0[j] = -1E100;
        }
        xbornecompo = 1000.0;
    }
    for (int gr = 1; gr < header.ngroupes + 1; gr++)
    {
        kp0 = kp;
        if (groupe[gr].type == 0)
        {
            if (groupe[gr].priormutmoy.constant)
            {
                if (groupe[gr].priorsnimoy.constant)
                {
                    pmut = (long double)(groupe[gr].mutmoy + groupe[gr].snimoy);
                    for (int j = 0; j < npar; j++)
                    {
                        if (scenario[rt.scenchoisi[0] - 1].histparam[numpar[0][j]].category < 2)
                        {
                            for (int i = 0; i < n; i++)
                            {
                                kscen = 0;
                                while (rt.enrsel[i].numscen != rt.scenchoisi[kscen])kscen++;
                                xx[i][kp] = pmut * (long double)rt.enrsel[i].param[numpar[kscen][j]];
                            }
                            kp++;
                        }
                    }
                }
                else
                {
                    kk = 0;
                    while (not ((header.mutparam[kk].groupe == gr)and (header.mutparam[kk].category == 2))) kk++;
                    for (int j = 0; j < npar; j++)
                    {
                        if (scenario[rt.scenchoisi[0] - 1].histparam[numpar[0][j]].category < 2)
                        {
                            for (int i = 0; i < n; i++)
                            {
                                kscen = 0;
                                while (rt.enrsel[i].numscen != rt.scenchoisi[kscen])kscen++;
                                pmut = (long double)groupe[gr].mutmoy + (long double)rt.enrsel[i].param[numpar[kscen][npar + kk]];
                                xx[i][kp] = pmut * (long double)rt.enrsel[i].param[numpar[kscen][j]];
                            }
                            kp++;
                        }
                    }
                }
            }
            else
            {
                if (groupe[gr].priorsnimoy.constant)
                {
                    kk = 0;
                    while (not ((header.mutparam[kk].groupe == gr)and (header.mutparam[kk].category == 0))) kk++;
                    for (int j = 0; j < npar; j++)
                    {
                        if (scenario[rt.scenchoisi[0] - 1].histparam[numpar[0][j]].category < 2)
                        {
                            for (int i = 0; i < n; i++)
                            {
                                kscen = 0;
                                while (rt.enrsel[i].numscen != rt.scenchoisi[kscen])kscen++;
                                pmut = (long double)rt.enrsel[i].param[numpar[kscen][npar + kk]] + (long double) groupe[gr].snimoy;
                                xx[i][kp] = pmut * (long double)rt.enrsel[i].param[numpar[kscen][j]];
                            }
                            kp++;
                        }
                    }
                }
                else
                {
                    kk = 0;
                    while (not ((header.mutparam[kk].groupe == gr)and (header.mutparam[kk].category == 0))) kk++;
                    qq = 0;
                    while (not ((header.mutparam[qq].groupe == gr)and (header.mutparam[qq].category == 2))) qq++;
                    cout << "dans recalparamC  kk=" << kk << "  qq=" << qq << "\n";
                    for (int j = 0; j < npar; j++)
                    {
                        if (scenario[rt.scenchoisi[0] - 1].histparam[numpar[0][j]].category < 2)
                        {
                            for (int i = 0; i < n; i++)
                            {
                                kscen = 0;
                                while (rt.enrsel[i].numscen != rt.scenchoisi[kscen])kscen++;
                                pmut = (long double)(rt.enrsel[i].param[numpar[kscen][npar + kk]] + rt.enrsel[i].param[numpar[kscen][npar + qq]]);
                                xx[i][kp] = pmut * (long double)rt.enrsel[i].param[numpar[kscen][j]];
                            }
                            kp++;
                        }
                    }
                }
            }
        }
        if (groupe[gr].type == 1)
        {
            if (groupe[gr].priormusmoy.constant)
            {
                pmut = (long double)groupe[gr].musmoy;
                for (int j = 0; j < npar; j++)
                {
                    if (scenario[rt.scenchoisi[0] - 1].histparam[numpar[0][j]].category < 2)
                    {
                        for (int i = 0; i < n; i++)
                        {
                            kscen = 0;
                            while (rt.enrsel[i].numscen != rt.scenchoisi[kscen])kscen++;
                            xx[i][kp] = pmut * (long double)rt.enrsel[i].param[numpar[kscen][j]];
                        }
                        kp++;
                    }
                }
            }
            else
            {
                kk = 0;
                while (not ((header.mutparam[kk].groupe == gr)and (header.mutparam[kk].category == 3))) kk++;
                for (int j = 0; j < npar; j++)
                {
                    if (scenario[rt.scenchoisi[0] - 1].histparam[numpar[0][j]].category < 2)
                    {
                        for (int i = 0; i < n; i++)
                        {
                            kscen = 0;
                            while (rt.enrsel[i].numscen != rt.scenchoisi[kscen])kscen++;
                            pmut = (long double)rt.enrsel[i].param[numpar[kscen][npar + kk]];
                            xx[i][kp] = pmut * (long double)rt.enrsel[i].param[numpar[kscen][j]];
                        }
                        kp++;
                    }
                }
            }
        }
        if (deltanul)
        {
            for (int i = 0; i < n; i++)
            {
                for (int j = kp0; j < kp; j++)
                {
                    alpsimrat[i][j] = xx[i][j];
                }
            }
            for (int i = 0; i < n; i++) delete [] xx[i];
            delete [] xx;
            return;
        }
        switch (numtransf)
        {
        case 1: //no transform
            for (int i = 0; i < n; i++)
            {
                for (int j = kp0; j < kp; j++)
                {
                    alpsimrat[i][j] = xx[i][j];
                }
            }
            break;
        case 2: // log transform
            for (int i = 0; i < n; i++)
            {
                for (int j = kp0; j < kp; j++)
                {
                    if (xx[i][j] <= 0.00000000001) xx[i][j] = 1E-11;
                    alpsimrat[i][j] = log(xx[i][j]);
                }
            }
            break;
        case 3: //logit transform
            for (int i = 0; i < n; i++)
            {
                for (int j = kp0; j < kp; j++)
                {
                    if (xx[i][j] < parmincompo0[j]) parmincompo0[j] = xx[i][j];
                    if (xx[i][j] > parmaxcompo0[j]) parmaxcompo0[j] = xx[i][j];
                }
            }
            for (int j = kp0; j < kp; j++)
            {
                parmincompo[j] = 0.95 * parmincompo0[j];
                parmaxcompo[j] = 1.05 * parmaxcompo0[j];
                diffcompo[j] = parmaxcompo[j] - parmincompo[j];
            }
            for (int i = 0; i < n; i++)
            {
                for (int j = kp0; j < kp; j++)
                {
                    if (xx[i][j] <= parmincompo0[j]) alpsimrat[i][j] = -xbornecompo;
                    else if (xx[i][j] >= parmaxcompo0[j]) alpsimrat[i][j] = xbornecompo;
                    else if (diffcompo[j] > 0.0) alpsimrat[i][j] = log((xx[i][j] - parmincompo[j]) / (parmaxcompo[j] - xx[i][j]));
                    else alpsimrat[i][j] = 0.0;
                }
            }
            break;
        case 4: //log(tg) transform
            for (int i = 0; i < n; i++)
            {
                for (int j = kp0; j < kp; j++)
                {
                    if (xx[i][j] < parmincompo0[j]) parmincompo0[j] = xx[i][j];
                    if (xx[i][j] > parmaxcompo0[j]) parmaxcompo0[j] = xx[i][j];
                }
            }
            for (int j = kp0; j < kp; j++)
            {
                parmincompo[j] = 0.95 * parmincompo0[j];
                parmaxcompo[j] = 1.05 * parmaxcompo0[j];
                diffcompo[j] = parmaxcompo[j] - parmincompo[j];
            }
            for (int i = 0; i < n; i++)
            {
                for (int j = kp0; j < kp; j++)
                {
                    if (xx[i][j] <= parmincompo0[j]) alpsimrat[i][j] = -xbornecompo;
                    else if (xx[i][j] >= parmaxcompo0[j]) alpsimrat[i][j] = xbornecompo;
                    else if (diffcompo[j] > 0.0) alpsimrat[i][j] = log(tan((xx[i][j] - parmincompo[j]) * c / diffcompo[j]));
                    else alpsimrat[i][j] = 0.0;
                }
            }
            break;
        }
    }
    for (int i = 0; i < n; i++) delete [] xx[i];
    delete [] xx;
}

/**
 * recalcule les valeurs de paramètres en fonction de la transformation choisie
 * pour les paramètres scaled
 */
void recalparamS(int n)
{
    //long double coefmarge=0.000000001,marge,pmut;
    long double **xx, Ne;
    int kk, kp, kscen, nNe = 0;
    alpsimrat = new long double*[n];
    for (int i = 0; i < n; i++) alpsimrat[i] = new long double[nparscaled];
    xx = new long double*[n];
    for (int i = 0; i < n; i++) xx[i] = new long double[nparscaled];
    for (int i = 0; i < n; i++)
    {
        kscen = 0;
        while (rt.enrsel[i].numscen != rt.scenchoisi[kscen])kscen++;
        ScenarioC& monscen = scenario[rt.scenchoisi[kscen] - 1];
        Ne = 0.0;
        nNe = 0;
        for (int j = 0; j < monscen.npop; j++)
        {
            for (int ievent = 0; ievent < monscen.nevent; ievent++)
            {
                if ((monscen.event[ievent].action == 'E')and (monscen.event[ievent].pop == j + 1))
                {
                    kk = 0;
                    while ((kk < monscen.nn0)and (monscen.histparam[kk].name != monscen.ne0[j].name)) kk++;
                    if (monscen.histparam[kk].prior.constant)
                        Ne += (long double)monscen.histparam[kk].prior.mini;
                    else Ne += (long double)rt.enrsel[i].param[numpar[kscen][kk]];
                    nNe++;
                    //if (i<10) cout<<"Ne = "<<Ne<<"\n";
                }
            }
        }
        Ne = Ne / (long double)nNe;
        kp = 0;
        for (int j = 0; j < npar; j++)
        {
            if (monscen.histparam[numpar[kscen][j]].category < 2)
            {
                xx[i][kp] = (long double)rt.enrsel[i].param[numpar[kscen][j]] / Ne;
                kp++;
            }
        }
        /*if (i<10) {
          cout<<"--->Ne="<<Ne;
          for (int j=0;j<npar;j++) cout<<"  "<<rt.enrsel[i].param[numpar[kscen][j]];
          cout<<"     ";
          for (int j=0;j<nparscaled;j++) cout<<"  "<<xx[i][j];
          cout<<"\n";
          }*/
    }
    if (numtransf > 2)
    {
        parminscaled = new long double[nparamcom];
        parmaxscaled = new long double[nparamcom];
        diffscaled = new long double[nparamcom];
        parminscaled0 = new long double[nparamcom];
        parmaxscaled0 = new long double[nparamcom];
        xbornescaled = 1000.0;
    }
    if (deltanul)
    {
        for (int i = 0; i < n; i++)
        {
            for (int j = 0; j < nparscaled; j++)
            {
                alpsimrat[i][j] = xx[i][j];
            }
        }
        for (int i = 0; i < n; i++) delete [] xx[i];
        delete [] xx;
        return;
    }
    switch (numtransf)
    {
    case 1: //no transform
        for (int i = 0; i < n; i++)
        {
            for (int j = 0; j < nparscaled; j++)
            {
                alpsimrat[i][j] = xx[i][j];
            }
        }
        break;
    case 2: // log transform
        for (int i = 0; i < n; i++)
        {
            for (int j = 0; j < nparscaled; j++)
            {
                if (xx[i][j] <= 0.00000000001) xx[i][j] = 1E-11;
                alpsimrat[i][j] = log(xx[i][j]);
            }
        }
        break;
    case 3: //logit transform
        for (int j = 0; j < nparscaled; j++)
        {
            parminscaled0[j] = 1E100;
            parmaxscaled0[j] = -1E100;
        }
        for (int i = 0; i < n; i++)
        {
            for (int j = 0; j < nparscaled; j++)
            {
                if (xx[i][j] < parminscaled0[j]) parminscaled0[j] = xx[i][j];
                if (xx[i][j] > parmaxscaled0[j]) parmaxscaled0[j] = xx[i][j];
            }
        }
        for (int j = 0; j < nparscaled; j++)
        {
            parminscaled[j] = 0.95 * parminscaled0[j];
            parmaxscaled[j] = 1.05 * parmaxscaled0[j];
            diffscaled[j] = parmaxscaled[j] - parminscaled[j];
        }
        for (int i = 0; i < n; i++)
        {
            for (int j = 0; j < nparscaled; j++)
            {
                if (xx[i][j] <= parminscaled0[j]) alpsimrat[i][j] = -xbornescaled;
                else if (xx[i][j] >= parmaxscaled0[j]) alpsimrat[i][j] = xbornescaled;
                else if (diffscaled[j] > 0.0) alpsimrat[i][j] = log((xx[i][j] - parminscaled[j]) / (parmaxscaled[j] - xx[i][j]));
                else alpsimrat[i][j] = 0.0;
            }
        }
        break;
    case 4: //log(tg) transform
        for (int j = 0; j < nparscaled; j++)
        {
            parminscaled0[j] = 1E100;
            parmaxscaled0[j] = -1E100;
        }
        for (int i = 0; i < n; i++)
        {
            for (int j = 0; j < nparscaled; j++)
            {
                if (xx[i][j] < parminscaled0[j]) parminscaled0[j] = xx[i][j];
                if (xx[i][j] > parmaxscaled0[j]) parmaxscaled0[j] = xx[i][j];
            }
        }
        for (int j = 0; j < nparscaled; j++)
        {
            parminscaled[j] = 0.95 * parminscaled0[j];
            parmaxscaled[j] = 1.05 * parmaxscaled0[j];
            diffscaled[j] = parmaxscaled[j] - parminscaled[j];
        }
        for (int i = 0; i < n; i++)
        {
            for (int j = 0; j < nparscaled; j++)
            {
                if (xx[i][j] <= parminscaled0[j]) alpsimrat[i][j] = -xbornescaled;
                else if (xx[i][j] >= parmaxscaled0[j]) alpsimrat[i][j] = xbornescaled;
                else if (diffscaled[j] > 0.0) alpsimrat[i][j] = log(tan((xx[i][j] - parminscaled[j]) * c / diffscaled[j]));
                else alpsimrat[i][j] = 0.0;
            }
        }
        break;
    }
    for (int i = 0; i < n; i++) delete [] xx[i];
    delete [] xx;
}

/**
 * effectue le remplissage de la matrice matX0, du vecteur des poids vecW et
 * de la matrice des paramètres parsim (éventuellement transformés)
 */
void rempli_mat(int n, float* stat_obs)
{
    int icc;
    long double delta, som, x, *var_statsel, nn;
    long double *sx, *sx2, *mo;
    nn = (long double)n;
    delta = rt.enrsel[n - 1].dist;
    deltanul = (delta == 0.0);
    if (deltanul)
    {
        cout << "delta = " << delta << "    Local regression is impossible and useless anyway...\n";
        return;
    }
    //printf("delta = %12.8Lf    (%12.8Lf,%12.8Lf)\n",delta,rt.enrsel[n-2].dist,rt.enrsel[n].dist);
    sx = new long double[rt.nstat];
    sx2 = new long double[rt.nstat];
    mo = new long double[rt.nstat];
    var_statsel = new long double[rt.nstat];
    for (int i = 0; i < rt.nstat; i++)
    {
        sx[i] = 0.0;
        sx2[i] = 0.0;
        mo[i] = 0.0;
    }
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < rt.nstat; j++)
        {
            x = (long double)rt.enrsel[i].stat[j];
            sx[j] += x;
            sx2[j] += x * x;
        }
    }
    nstatOKsel = 0;
    for (int j = 0; j < rt.nstat; j++)
    {
        var_statsel[j] = (sx2[j] - sx[j] / nn * sx[j]) / (nn - 1.0);
        if (var_statsel[j] > 1E-20) nstatOKsel++;
        mo[j] = sx[j] / nn;
    }
    //for (int j=0;j<rt.nstat;j++) printf("stat[%3d]  moy=%12.8Lf   var=%16Le\n",j,mo[j],var_statsel[j]);
    //cout<<"nstat="<<rt.nstat<<"     nstatOK="<<nstatOKsel<<"\n";
    matX0 = new long double*[n];
    for (int i = 0; i < n; i++)matX0[i] = new long double[nstatOKsel];
    vecW = new long double[n];
    //cout <<"hello\n";
    som = 0.0;
    for (int i = 0; i < n; i++)
    {
        icc = -1;
        for (int j = 0; j < rt.nstat; j++)
        {
            if (var_statsel[j] > 1E-20)
            {
                icc++;
                matX0[i][icc] = (long double)(rt.enrsel[i].stat[j] - stat_obs[j]) / sqrt(var_statsel[j]);
                //if ((i<10)) printf("enrsel[%d] stat[%d] %16.10f  %16.10f  %16.10Lf    %16.10Lf\n",i,j,rt.enrsel[i].stat[j],stat_obs[j],var_statsel[j],matX0[i][icc]);
            }
        }
        //cout<<"\n";
        x = rt.enrsel[i].dist / delta;
        vecW[i] = (1.5 / delta) * (1.0 - x * x);
        if (vecW[i] != vecW[i])
        {
            cout << "vecW[" << i << "]=" << vecW[i] << "\n";
            exit(1);
        }
        som = som + vecW[i];
    }
    for (int i = 0; i < n; i++) vecW[i] /= som;
    delete [] sx;
    delete [] sx2;
    delete [] mo;
    delete [] var_statsel;

    //for (int i=0;i<10;i++) cout<<vecW[i]<<"  ";
    //cout<<"\n";
}

void delete_mat(int n)
{
    for (int i = 0; i < n; i++) delete [] matX0[i];
    delete [] matX0;
    delete [] vecW;
}

void rempli_parsim(int n, int npa)
{
    parsim = new long double*[n];
    for (int i = 0; i < n; i++) parsim[i] = new long double[npa];
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < npa; j++)parsim[i][j] = (long double)alpsimrat[i][j];
    }
}

int ecrimatL(string nomat, int n, int m, long double** A)
{
    cout << "\n" << nomat << "\n";
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < m; j++) cout << setiosflags(ios::fixed) << setw(13) << setprecision(9) << A[i][j] << "  ";
        cout << "\n";
    }
    cout << "\n";
    return 0;
}

int ecrimatD(string nomat, int n, int m, double** A)
{
    cout << "\n" << nomat << "\n";
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < m; j++) cout << setiosflags(ios::fixed) << setw(10) << setprecision(6) << A[i][j] << "  ";
        cout << "\n";
    }
    cout << "\n";
    return 0;
}

long double** cal_matC(int n)
{
    long double **matX, **matXT, **matA, **matB, **matAA, **matC, **matBB, kap, coeff;
    int err;
    matA = new long double*[nstatOKsel + 1];
    for (int j = 0; j < nstatOKsel + 1; j++) matA[j] = new long double[n];
    matX = new long double*[n];
    for (int i = 0; i < n; i++)
    {
        matX[i] = new long double[nstatOKsel + 1];
        matX[i][0] = 1.0;
        for (int j = 1; j < nstatOKsel + 1; j++) matX[i][j] = (long double)matX0[i][j - 1];
    }
    matB = new long double*[nstatOKsel + 1];
    for (int j = 0; j < nstatOKsel + 1; j++) matB[j] = new long double[nstatOKsel + 1];
    matBB = new long double*[nstatOKsel + 1];
    for (int j = 0; j < nstatOKsel + 1; j++) matBB[j] = new long double[nstatOKsel + 1];
    matXT = transposeL(n, nstatOKsel + 1, matX);
    for (int i = 0; i < n; i++)
    {
        if (vecW[i] != vecW[i])
        {
            cout << "vecW[" << i << "] = " << vecW[i] << "\n";
            exit(1);
        }
        for (int j = 0; j < nstatOKsel + 1; j++)
        {
            if (matXT[j][i] != matXT[j][i])
            {
                cout << "matXT[" << j << "][" << i << "] = " << matXT[j][i] << "\n";
                exit(1);
            }
            matA[j][i] = matXT[j][i] * vecW[i];
        }
    }
    for (int i = 0; i < nstatOKsel + 1; i++)
    {
        for (int j = 0; j < nstatOKsel + 1; j++)
            if (matA[i][j] != matA[i][j])
            {
                cout << "matA[" << i << "][" << j << "] = " << matA[i][j] << "\n";
                exit(1);
            }
    }
    //cout<<"matA OK\n";
    matAA = prodML(nstatOKsel + 1, n, nstatOKsel + 1, matA, matX);
    kap = kappa(nstatOKsel + 1, matAA);
    printf("kappa (AA) = %16Le", kap);
    if (kap > 1.0E99) cout << "   MATRICE SINGULIERE\n";
    cout << "\n";
    if (kap > 1.0E99) coeff = 1.0E-15; else coeff = 1.0E-21;
    err = inverse_Tik2(nstatOKsel + 1, matAA, matB, coeff);
    do
    {
        err = inverse_Tik2(nstatOKsel + 1, matAA, matB, coeff);
        if (err != 0) coeff *= sqrt(10.0);
    }
    while ((err != 0)and (coeff < 0.01));
    if (debuglevel == 11) cout << "cal_matC   err=" << err << "    coeff=" << coeff << "\n";
    for (int i = 0; i < nstatOKsel + 1; i++)
    {
        for (int j = 0; j < nstatOKsel + 1; j++)
            if (matB[i][j] != matB[i][j])
            {
                cout << "matB[" << i << "][" << j << "] = " << matB[i][j] << "\n";
                exit(1);
            }
    }
    //cout<<"matB OK\n";
    matC = prodML(nstatOKsel + 1, nstatOKsel + 1, n, matB, matA);
    libereL(nstatOKsel + 1, matA);
    libereL(n, matX);
    libereL(nstatOKsel + 1, matXT);
    libereL(nstatOKsel + 1, matB);
    libereL(nstatOKsel + 1, matBB);
    libereL(nstatOKsel + 1, matAA);
    return matC;
}

void local_regression(int n, int npa, long double** matC)
{
    beta = new long double*[nstatOKsel + 1];
    for (int j = 0; j < nstatOKsel + 1; j++) beta[j] = new long double[npa];
    for (int i = 0; i < nstatOKsel + 1; i++)
    {
        for (int j = 0; j < nstatOKsel + 1; j++)
            if (matC[i][j] != matC[i][j])
            {
                cout << "matC[" << i << "][" << j << "] = " << matC[i][j] << "\n";
                exit(1);
            }
    }
    beta = prodML(nstatOKsel + 1, n, npa, matC, parsim);
    /*cout<<"regression locale\n";
      for (int j=0;j<nstatOKsel+1;j++) {
      for (int k=0;k<npa;k++) cout<<"beta["<<j<<"]["<<k<<"] = "<<beta[j][k]<<"  ";
      cout<<"\n";
      }
      cout<<"\n";*/
    libereL(n, parsim);
}

/**
 * effectue la régression locale à partir de la matrice matX0 et le vecteur des poids vecW
 */
/*    void local_regression(int n, int npa) {
      long double **matX,**matXT,**matA,**matB,**matAA,**matC,**bet = NULL,**matBB,kap,mdiff,coeff;
      int err;
      matA = new long double*[nstatOKsel+1];
      for (int j=0;j<nstatOKsel+1;j++) matA[j] = new long double[n];
      matX = new long double*[n];
      for (int i=0;i<n;i++) {
      matX[i] = new long double[nstatOKsel+1];
      matX[i][0] = 1.0;
      for (int j=1;j<nstatOKsel+1;j++) matX[i][j] = (long double)matX0[i][j-1];
      }
      matB = new long double*[nstatOKsel+1];
      for (int j=0;j<nstatOKsel+1;j++) matB[j] = new long double[nstatOKsel+1];
      matBB = new long double*[nstatOKsel+1];
      for (int j=0;j<nstatOKsel+1;j++) matBB[j] = new long double[nstatOKsel+1];
      beta = new long double*[nstatOKsel+1];
      for (int j=0;j<nstatOKsel+1;j++) beta[j] = new long double[npa];
      //ecrimatL("matX0",10,nstatOKsel,matX0);
      //ecrimat("matX",10,10,matX);

      matXT = transposeL(n,nstatOKsel+1,matX);
      for (int i=0;i<n;i++) {
      for (int j=0;j<nstatOKsel+1;j++) matA[j][i] = matXT[j][i]*vecW[i];
      }
      //ecrimat("matA",10,10,matA);
      matAA = prodML(nstatOKsel+1,n,nstatOKsel+1,matA,matX);
      //ecrimatL("matAA",nstatOKsel+1,nstatOKsel+1,matAA);
      kap = kappa(nstatOKsel+1,matAA);
      //printf("kappa (AA) = %16Le",kap);
      if (kap>1.0E99) cout <<"   MATRICE SINGULIERE\n";
      //cout<<"\n";
      if (kap>1.0E99) coeff=1.0E-15; else coeff=1.0E-21;
      err = inverse_Tik2(nstatOKsel+1,matAA,matB,coeff);
      if (err==0) {
      matC = prodML(nstatOKsel+1,nstatOKsel+1,n,matB,matA);
      bet = prodML(nstatOKsel+1,n,npa,matC,parsim);
      }
                
                
      do {
      for (int i=0;i<nstatOKsel+1;i++) {for (int j=0;j<npa;j++) beta[i][j] = bet[i][j];}
      coeff *=sqrt(10.0);
      err = inverse_Tik2(nstatOKsel+1,matAA,matB,coeff);
//              cout<<"\n\n err="<<err<<"\n";
//              ecrimatL("matB = inv(matAA)",nstatOKsel+1,nstatOKsel+1,matB);
matC = prodML(nstatOKsel+1,nstatOKsel+1,n,matB,matA);
//              cout<<"apres matC\n";
bet = prodML(nstatOKsel+1,n,npa,matC,parsim);
mdiff = 0.0;
for (int i=0;i<nstatOKsel+1;i++) {for (int j=0;j<npa;j++) mdiff += fabs(beta[i][j] - bet[i][j])/fabs(beta[i][j] + bet[i][j]);}
mdiff /=(nstatOKsel+1)*(nstatOKsel+1);
//cout<<"coeff = 1E"<<log10(coeff)<<"    mdiff = "<<mdiff<<"\n";
//printf("coeff = %8Le   mdiff = %8.5Lf\n",coeff,mdiff);
} while ((mdiff>0.0002)and(coeff<0.01));
                
//for (int i=0;i<nstatOKsel+1;i++) {for (int j=0;j<npa;j++) beta[i][j] = bet[i][j];}
//ecrimatL("beta",nstatOKsel+1,npa,beta);
libereL(nstatOKsel+1,matA);
libereL(n,matX);
libereL(nstatOKsel+1,matXT);
libereL(nstatOKsel+1,matB);
libereL(nstatOKsel+1,matBB);
libereL(nstatOKsel+1,matAA);
libereL(nstatOKsel+1,matC);
libereL(n,parsim);
libereL(nstatOKsel+1,bet);
}*/

/**
 * calcule les phistars pour les paramètres originaux
 */
void calphistarO(int n, long double** phistar)
{
    //cout<<"debut de calphistarO\n";
    //int k,kk,qq;
    //long double pmut;
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < nparamcom; j++)
        {
            phistar[i][j] = alpsimrat[i][j];
            //if ((i<100)and(j==0)) cout<< phistar[i][j]<<"   ";
            for (int k = 0; k < nstatOKsel; k++) phistar[i][j] -= matX0[i][k] * beta[k + 1][j];
            //if((i<100)and(j==0)) cout<< phistar[i][j]<<"\n";
            switch (numtransf)
            {
            case 1: break;
            case 2: if (phistar[i][j] < 100.0) phistar[i][j] = exp(phistar[i][j]); else phistar[i][j] = exp(100.0);
                break;
            case 3: if (phistar[i][j] <= -xborne) phistar[i][j] = parmin0[j];
                else if (phistar[i][j] >= xborne) phistar[i][j] = parmax0[j];
                else phistar[i][j] = (exp(phistar[i][j]) * parmax[j] + parmin[j]) / (1.0 + exp(phistar[i][j]));
                if (phistar[i][j] < parmin0[j]) phistar[i][j] = parmin0[j];
                if (phistar[i][j] > parmax0[j]) phistar[i][j] = parmax0[j];
                break;
            case 4: if (phistar[i][j] <= -xborne) phistar[i][j] = parmin0[j];
                else if (phistar[i][j] >= xborne) phistar[i][j] = parmax0[j];
                else phistar[i][j] = parmin[j] + (diff[j] / c * atan(exp(phistar[i][j])));
                if (phistar[i][j] < parmin0[j]) phistar[i][j] = parmin0[j];
                if (phistar[i][j] > parmax0[j]) phistar[i][j] = parmax0[j];
                break;
            }
        }
    }
    if (numtransf > 2)
    {
        delete []parmin;
        delete []parmax;
        delete []diff;
        delete []parmin0;
        delete []parmax0;
    }
    //cout<<"nparcompo = "<<nparcompo<<"   k="<<k<<"\n";
    for (int i = 0; i < n; i++) delete []alpsimrat[i];
    delete []alpsimrat;
    //for (int i=0;i<n;i++) delete []matX0[i]; delete []matX0;
    for (int i = 0; i < nstatOKsel + 1; i++) delete []beta[i];
    delete []beta;
}


void copyphistar(int n, int nparam, long double** phistar)
{
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < nparam; j++)
        {
            phistar[i][j] = alpsimrat[i][j];
        }
    }
}

/**
 * calcule les phistars pour les paramètres composites
 */
void calphistarC(int n, long double** phistarcompo)
{
    //cout<<"debut de calphistarC\n";
    //int k,kk,qq;
    //long double ppp;
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < nparcompo; j++)
        {
            phistarcompo[i][j] = alpsimrat[i][j];
            //if (i<100) cout<< phistarcompo[i][j]<<"   ";
            for (int k = 0; k < nstatOKsel; k++) phistarcompo[i][j] -= matX0[i][k] * beta[k + 1][j];
            //ppp=phistarcompo[i][j];
            //for (int k=0;k<nstatOKsel;k++) if (beta[k+1][j]!=beta[k+1][j])cout<<"beta=nan  k="<<k<<"   j="<<j<<"    "<<beta[k+1][j]<<"\n";
            //if(i<100) cout<< phistarcompo[i][j]<<"   ";
            switch (numtransf)
            {
            case 1: break;
            case 2: if (phistarcompo[i][j] < 100.0) phistarcompo[i][j] = exp(phistarcompo[i][j]); else phistarcompo[i][j] = exp(100.0);
                break;
            case 3: if (phistarcompo[i][j] <= -xbornecompo) phistarcompo[i][j] = parmincompo0[j];
                else if (phistarcompo[i][j] >= xbornecompo) phistarcompo[i][j] = parmaxcompo0[j];
                else phistarcompo[i][j] = (exp(phistarcompo[i][j]) * parmaxcompo[j] + parmincompo[j]) / (1.0 + exp(phistarcompo[i][j]));
                //if(i<100) cout<< phistarcompo[i][j]<<"   ("<<parmincompo[j]<<","<<parmaxcompo[j]<<")\n";
                if (phistarcompo[i][j] < parmincompo0[j]) phistarcompo[i][j] = parmincompo0[j];
                if (phistarcompo[i][j] > parmaxcompo0[j]) phistarcompo[i][j] = parmaxcompo0[j];
                break;
            case 4: if (phistarcompo[i][j] <= -xbornecompo) phistarcompo[i][j] = parmincompo0[j];
                else if (phistarcompo[i][j] >= xbornecompo) phistarcompo[i][j] = parmaxcompo0[j];
                else phistarcompo[i][j] = parmincompo[j] + (diffcompo[j] / c * atan(exp(phistarcompo[i][j])));
                if (phistarcompo[i][j] < parmincompo0[j]) phistarcompo[i][j] = parmincompo0[j];
                if (phistarcompo[i][j] > parmaxcompo0[j]) phistarcompo[i][j] = parmaxcompo0[j];
                break;
            }
        }
    }
    if (numtransf > 2)
    {
        delete []parmincompo;
        delete []parmaxcompo;
        delete []diffcompo;
        delete []parmincompo0;
        delete []parmaxcompo0;
    }
    //cout<<"nparcompo = "<<nparcompo<<"   k="<<k<<"\n";
    for (int i = 0; i < n; i++) delete []alpsimrat[i];
    delete []alpsimrat;
    //for (int i=0;i<n;i++) delete []matX0[i]; delete []matX0;
    for (int i = 0; i < nstatOKsel + 1; i++) delete []beta[i];
    delete []beta;
}

/**
 * calcule les phistars pour les paramètres scaled
 */
void calphistarS(int n, long double** phistarscaled)
{
    //cout<<"debut de calphistarS\n";
    //int k,kk,qq;
    //long double pmut;
    //long double **phista;
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < nparscaled; j++)
        {
            phistarscaled[i][j] = alpsimrat[i][j];
            //if ((i<10)and(j==0)) cout<< phistarscaled[i][j]<<"   ";
            for (int k = 0; k < nstatOKsel; k++) phistarscaled[i][j] -= matX0[i][k] * beta[k + 1][j];
            //if ((i<10)and(j==0)) cout<< phistarscaled[i][j]<<"\n";
            switch (numtransf)
            {
            case 1: break;
            case 2: if (phistarscaled[i][j] < 100.0) phistarscaled[i][j] = exp(phistarscaled[i][j]); else phistarscaled[i][j] = exp(100.0);
                break;
            case 3: if (phistarscaled[i][j] < -xbornescaled) phistarscaled[i][j] = parminscaled0[j];
                else if (phistarscaled[i][j] > xbornescaled) phistarscaled[i][j] = parmaxscaled0[j];
                else phistarscaled[i][j] = (exp(phistarscaled[i][j]) * parmaxscaled[j] + parminscaled[j]) / (1.0 + exp(phistarscaled[i][j]));
                //if(i<10) cout<<"      "<< phistarscaled[i][j]<<"   ("<<parminscaled[j]<<","<<parmaxscaled[j]<<")\n";
                if (phistarscaled[i][j] < parminscaled0[j]) phistarscaled[i][j] = parminscaled0[j];
                if (phistarscaled[i][j] > parmaxscaled0[j]) phistarscaled[i][j] = parmaxscaled0[j];
                break;
            case 4: if (phistarscaled[i][j] < -xbornescaled) phistarscaled[i][j] = parminscaled0[j];
                else if (phistarscaled[i][j] > xbornescaled) phistarscaled[i][j] = parmaxscaled0[j];
                else phistarscaled[i][j] = parminscaled[j] + (diffscaled[j] / c * atan(exp(phistarscaled[i][j])));
                if (phistarscaled[i][j] < parminscaled0[j]) phistarscaled[i][j] = parminscaled0[j];
                if (phistarscaled[i][j] > parmaxscaled0[j]) phistarscaled[i][j] = parmaxscaled0[j];
                break;
            }
        }
        //if(i<100) cout<<"\n";
    }
    if (numtransf > 2)
    {
        delete []parminscaled;
        delete []parmaxscaled;
        delete []diffscaled;
        delete []parminscaled0;
        delete []parmaxscaled0;
    }
    //cout<<"nparcompo = "<<nparcompo<<"   k="<<k<<"\n";
    for (int i = 0; i < n; i++) delete []alpsimrat[i];
    delete []alpsimrat;
    //for (int i=0;i<n;i++) delete []matX0[i]; delete []matX0;
    for (int i = 0; i < nstatOKsel + 1; i++) delete []beta[i];
    delete []beta;
}

void det_nomparam()
{
    int k;
    string pp;
    if (original)
    {
        nomparamO = new string[nparamcom];
        enteteO = "scenario";
        for (int j = 0; j < npar; j++)
        {
            nomparamO[j] = string(rt.histparam[rt.scenchoisi[0] - 1][numpar[0][j]].name);
            enteteO += centre(nomparamO[j], 15);
        }
        for (int j = npar; j < nparamcom; j++)
        {
            nomparamO[j] = string(rt.mutparam[j - npar].name);
            enteteO += centre(nomparamO[j], 15);
        }
    }
    if (composite)
    {
        nomparamC = new string[nparcompo];
        enteteC = "scenario";
        k = 0;
        for (int gr = 1; gr < header.ngroupes + 1; gr++)
        {
            for (int j = 0; j < npar; j++)
            {
                if (rt.histparam[rt.scenchoisi[0] - 1][numpar[0][j]].category < 2)
                {
                    pp = rt.histparam[rt.scenchoisi[0] - 1][numpar[0][j]].name;
                    if (groupe[gr].type == 0) pp = pp + "(u+sni)_" + IntToString(gr) + " ";
                    if (groupe[gr].type == 1) pp = pp + "useq_" + IntToString(gr) + " ";
                    nomparamC[k] = pp;
                    k++;
                    enteteC += centre(pp, 16);
                }
            }
        }
    }
    if (scaled)
    {
        nomparamS = new string[nparscaled];
        enteteS = "scenario";
        k = 0;
        for (int j = 0; j < npar; j++)
        {
            if (rt.histparam[rt.scenchoisi[0] - 1][numpar[0][j]].category < 2)
            {
                pp = rt.histparam[rt.scenchoisi[0] - 1][numpar[0][j]].name;
                pp = pp + "/Mean(N)";
                nomparamS[k] = pp;
                k++;
                enteteS += centre(pp, 16);
            }
        }
    }
}

string fmLB(long double x, int largeur, int precision)
{
    string s, s1;
    size_t pos;
    int lon, ifacteur = 0, ix;
    long double x1, facteur = 1.0, y = 0.5;
    for (int i = 0; i < precision; i++) y /= 10.0;
    //bool retenue;
    //cout<<"fmLB        x="<<x<<"`\n";
    if (x > 100.0)
    {
        ix = (int)x;
        s = LongDoubleToString(ix);
    }
    else if (x >= 1.0)
    {
        s = LongDoubleToString(x + y);
        if (s.find(".") == string::npos) s = s + ".0";
        s = s + "000000000000000";
        pos = s.find(".");
        //s1=s.substr(pos+precision+1,1);
        //retenue = (s1=="5")or(s1=="5")
        s = s.substr(0, pos + precision + 1);
    }
    if (x < 1.0)
    {
        if (x == 0.0)
        {
            s = "0.";
            for (int i = 0; i < precision; i++) s = +"0";
        }
        else
        {
            do
            {
                facteur *= 10.0;
                ifacteur++;
                x1 = facteur * x;
            }
            while (x1 < 1.0);
            s = LongDoubleToString(x1 + y);
            s = s + "000000000000000";
            pos = s.find(".");
            s = s.substr(0, pos + precision + 1);
            s = s + "e-" + IntToString(ifacteur);
        }
    }
    lon = s.length();
    while (lon < largeur)
    {
        s = " " + s;
        lon++;
    }
    //cout<<"                             s final="<<s<<"\n";
    return s;
}

/**
 * effectue la sauvegarde des phistars dans le fichier path/ident/phistar.txt
 */
void savephistar(int n, string path, string ident, long double** phistar, long double** phistarcompo, long double** phistarscaled)
{
    string nomphistar;
    if (original)
    {
        nomphistar = path + ident + "_phistar_original.txt";
        ofstream f1;
        f1.open(nomphistar.c_str());
        f1 << enteteO << "\n";
        for (int i = 0; i < n; i++)
        {
            f1 << setiosflags(ios::fixed) << setw(4) << rt.enrsel[i].numscen;
            for (int j = 0; j < nparamcom; j++) f1 << fmLB(phistar[i][j], 15, 3);
            f1 << "\n";
        }
        f1.close();
    }
    if (composite)
    {
        nomphistar = path + ident + "_phistar_composite.txt";
        ofstream f1;
        f1.open(nomphistar.c_str());
        f1 << enteteC << "\n";
        for (int i = 0; i < n; i++)
        {
            f1 << setiosflags(ios::fixed) << setw(4) << rt.enrsel[i].numscen;
            for (int j = 0; j < nparcompo; j++)
                if (phistarcompo[i][j] < 0.0)
                {
                    cout << "i=" << i << "   j=" << j << "   " << phistarcompo[i][j] << "\n";
                    exit(1);
                }
            for (int j = 0; j < nparcompo; j++) f1 << fmLB(phistarcompo[i][j], 16, 3);
            f1 << "\n";
        }
        f1.close();
    }
    if (scaled)
    {
        nomphistar = path + ident + "_phistar_scaled.txt";
        ofstream f1;
        f1.open(nomphistar.c_str());
        f1 << enteteS << "\n";
        for (int i = 0; i < n; i++)
        {
            cout << "i=" << i << "\n";
            f1 << setiosflags(ios::fixed) << setw(4) << rt.enrsel[i].numscen;
            cout << "1\n";
            cout << "nparscaled=" << nparscaled << "\n";
            for (int j = 0; j < nparscaled; j++)
            {
                cout << "j=" << j << "\n";
                if (phistarscaled[i][j] < 0.0)
                {
                    cout << "i=" << i << "   j=" << j << "   " << phistarscaled[i][j] << "\n";
                    exit(1);
                }
            }
            cout << "2\n";
            for (int j = 0; j < nparscaled; j++) f1 << fmLB(phistarscaled[i][j], 16, 3);
            cout << "3\n";
            f1 << "\n";
        }
        f1.close();
        cout << "fin scaled\n";
    }
}

/**
 * calcul de la densité par noyau gaussien
 * x : vecteur des abcisses
 * y vecteur des
 */


/**
 * lit les paramètres originaux des enregistrements simulés pour l'établissement des distributions a priori'
 */
void lisimparO()
{
    int bidon, m, nr;
    bool scenOK;
    //long double pmut;
    if (nsimpar > rt.nrec) nsimpar = rt.nrec;
    nr = 0;
    for (int i = 0; i < rt.nscenchoisi; i++) nr += rt.nrecscen[rt.scenchoisi[i] - 1];
    if (nsimpar > nr) nsimpar = nr;
    simpar = new long double*[nsimpar];
    enregC enr;
    enr.stat = vector<float>(rt.nstat);
    enr.param = vector<float>(rt.nparamax);
    int i = 0;
    rt.openfile2();
    while (i < nsimpar)
    {
        bidon = rt.readrecord(&enr);
        if (bidon != 0) cout << "problème dans la lecture du reftable\n";
        m = 0;
        scenOK = false;
        while ((not scenOK)and (m < rt.nscenchoisi))
        {
            scenOK = (enr.numscen == rt.scenchoisi[m]);
            if (not scenOK) m++;
        }
        if (scenOK)
        {
            simpar[i] = new long double[nparamcom];
            for (int j = 0; j < nparamcom; j++) simpar[i][j] = enr.param[numpar[m][j]];
            i++;
        }
    }
    rt.closefile();
}

/**
 * calcule les paramètres composites des enregistrements simulés pour l'établissement des distributions a priori'
 */
void lisimparC()
{
    int bidon, m, k, kk, qq, nr;
    bool scenOK;
    long double pmut;
    if (nsimpar > rt.nrec) nsimpar = rt.nrec;
    nr = 0;
    for (int i = 0; i < rt.nscenchoisi; i++) nr += rt.nrecscen[rt.scenchoisi[i] - 1];
    if (nsimpar > nr) nsimpar = nr;
    simparcompo = new long double*[nsimpar];
    enregC enr;
    enr.stat = vector<float>(rt.nstat);
    enr.param = vector<float>(rt.nparamax);
    int i = 0;
    rt.openfile2();
    while (i < nsimpar)
    {
        bidon = rt.readrecord(&enr);
        if (bidon != 0) cout << "problème dans la lecture du reftable\n";
        m = 0;
        scenOK = false;
        while ((not scenOK)and (m < rt.nscenchoisi))
        {
            scenOK = (enr.numscen == rt.scenchoisi[m]);
            if (not scenOK) m++;
        }
        if (scenOK)
        {
            simparcompo[i] = new long double[nparcompo];
            k = 0;
            for (int gr = 1; gr < header.ngroupes + 1; gr++)
            {
                if (groupe[gr].type == 0)
                {
                    if (groupe[gr].priormutmoy.constant)
                    {
                        if (groupe[gr].priorsnimoy.constant)
                        {
                            pmut = (long double)(groupe[gr].mutmoy + groupe[gr].snimoy);
                            for (int j = 0; j < npar; j++)
                            {
                                if (scenario[rt.scenchoisi[m] - 1].histparam[numpar[m][j]].category < 2)
                                {
                                    simparcompo[i][k] = pmut * (long double)enr.param[numpar[m][j]];
                                    k++;
                                }
                            }
                        }
                        else
                        {
                            kk = 0;
                            while (not ((header.mutparam[kk].groupe == gr)and (header.mutparam[kk].category == 2))) kk++;
                            for (int j = 0; j < npar; j++)
                            {
                                if (scenario[rt.scenchoisi[m] - 1].histparam[numpar[m][j]].category < 2)
                                {
                                    pmut = (long double)groupe[gr].mutmoy + (long double)enr.param[numpar[m][npar + kk]];
                                    simparcompo[i][k] = pmut * (long double)enr.param[numpar[m][j]];
                                    k++;
                                }
                            }
                        }
                    }
                    else
                    {
                        if (groupe[gr].priorsnimoy.constant)
                        {
                            kk = 0;
                            while (not ((header.mutparam[kk].groupe == gr)and (header.mutparam[kk].category == 0))) kk++;
                            for (int j = 0; j < npar; j++)
                            {
                                if (scenario[rt.scenchoisi[m] - 1].histparam[numpar[m][j]].category < 2)
                                {
                                    pmut = (long double)(enr.param[numpar[m][npar + kk]] + groupe[gr].snimoy);
                                    simparcompo[i][k] = pmut * (long double)enr.param[numpar[m][j]];
                                    k++;
                                }
                            }
                        }
                        else
                        {
                            kk = 0;
                            while (not ((header.mutparam[kk].groupe == gr)and (header.mutparam[kk].category == 0))) kk++;
                            qq = 0;
                            while (not ((header.mutparam[qq].groupe == gr)and (header.mutparam[qq].category == 2))) qq++;
                            for (int j = 0; j < npar; j++)
                            {
                                if (scenario[rt.scenchoisi[m] - 1].histparam[numpar[m][j]].category < 2)
                                {
                                    pmut = (long double)enr.param[numpar[m][npar + kk]] + (long double)enr.param[numpar[m][npar + qq]];
                                    simparcompo[i][k] = pmut * (long double)enr.param[numpar[m][j]];
                                    k++;
                                }
                            }
                        }
                    }
                }
                if (groupe[gr].type == 1)
                {
                    if (groupe[gr].priormusmoy.constant)
                    {
                        pmut = (long double)groupe[gr].musmoy;
                        for (int j = 0; j < npar; j++)
                        {
                            if (scenario[rt.scenchoisi[m] - 1].histparam[numpar[m][j]].category < 2)
                            {
                                simparcompo[i][k] = pmut * (long double)enr.param[numpar[m][j]];
                                k++;
                            }
                        }
                    }
                    else
                    {
                        kk = 0;
                        while (not ((header.mutparam[kk].groupe == gr)and (header.mutparam[kk].category == 3))) kk++;
                        for (int j = 0; j < npar; j++)
                        {
                            if (scenario[rt.scenchoisi[m] - 1].histparam[numpar[m][j]].category < 2)
                            {
                                pmut = (long double)enr.param[numpar[m][npar + kk]];
                                simparcompo[i][k] = pmut * (long double)enr.param[numpar[m][j]];
                                k++;
                            }
                        }
                    }
                }
            }
            i++;
        }
    }
    rt.closefile();
}

/**
 * calcule les paramètres scaled des enregistrements simulés pour l'établissement des distributions a priori'
 */
void lisimparS()
{
    int bidon, m, k, nr, kk, nNe = 0;
    bool scenOK;
    long double Ne;
    if (nsimpar > rt.nrec) nsimpar = rt.nrec;
    nr = 0;
    for (int i = 0; i < rt.nscenchoisi; i++) nr += rt.nrecscen[rt.scenchoisi[i] - 1];
    if (nsimpar > nr) nsimpar = nr;
    simparscaled = new long double*[nsimpar];
    enregC enr;
    enr.stat = vector<float>(rt.nstat);
    enr.param = vector<float>(rt.nparamax);
    int i = 0;
    rt.openfile2();
    while (i < nsimpar)
    {
        bidon = rt.readrecord(&enr);
        if (bidon != 0) cout << "problème à la lecture de la reftable\n";
        m = 0;
        scenOK = false;
        while ((not scenOK)and (m < rt.nscenchoisi))
        {
            scenOK = (enr.numscen == rt.scenchoisi[m]);
            if (not scenOK) m++;
        }
        if (scenOK)
        {
            //cout<<"scenOK   m="<<m<<"\n";
            simparscaled[i] = new long double[nparscaled];
            Ne = 0.0;
            nNe = 0;
            for (int j = 0; j < scenario[rt.scenchoisi[m] - 1].npop; j++)
            {
                for (int ievent = 0; ievent < scenario[rt.scenchoisi[m] - 1].nevent; ievent++)
                {
                    if ((scenario[rt.scenchoisi[m] - 1].event[ievent].action == 'E')and (scenario[rt.scenchoisi[m] - 1].event[ievent].pop == j + 1))
                    {
                        //cout<<"nn0="<<scenario[rt.scenchoisi[m]-1].nn0<<"\n";
                        kk = 0;
                        while ((kk < scenario[rt.scenchoisi[m] - 1].nn0)and (scenario[rt.scenchoisi[m] - 1].histparam[kk].name != scenario[rt.scenchoisi[m] - 1].ne0[j].name))
                        {
                            //cout<<"scenario[rt.scenchoisi[m]-1].ne0[j].name : ---"<<scenario[rt.scenchoisi[m]-1].ne0[j].name<<"---\n";
                            //cout<<"scenario[rt.scenchoisi[m]-1].histparam["<<kk<<"].name : ---"<<scenario[rt.scenchoisi[m]-1].histparam[kk].name<<"---\n";
                            kk++;
                        }
                        //cout<<"kk="<<kk<<"\n";
                        if (scenario[rt.scenchoisi[m] - 1].histparam[kk].prior.constant)
                            Ne += (long double)scenario[rt.scenchoisi[m] - 1].histparam[kk].prior.mini;
                        else Ne += (long double)enr.param[numpar[m][kk]];
                        nNe++;
                        //cout<<"Ne = "<<Ne<<"  "<<((int)Ne % 3)<<"\n";
                    }
                }
            }
            Ne = Ne / (long double)nNe;
            //cout<<"Ne = "<<Ne<<"  nNe="<<nNe<<"\n";
            k = 0; //cout<<"Ne = "<<Ne<<"  "<<((int)Ne % 3)<<"\n";
            for (int j = 0; j < npar; j++)
            {
                if (scenario[rt.scenchoisi[m] - 1].histparam[numpar[m][j]].category < 2)
                {
                    //cout<<"i="<<i<<"   nsimpar="<<nsimpar<<"   k="<<k<<"\n";
                    simparscaled[i][k] = (long double)enr.param[numpar[m][j]] / Ne;
                    k++;
                }
            }
            i++;
        }
    }
    rt.closefile();
}

/**
 * calcule la densité à partir de la loi fournie dans le prior
 */
long double* caldensexact(int ncl, long double* x, PriorC pr)
{
    long double *dens, xb, som;
    dens = new long double[ncl];
    //cout<<"caldensexact\n";
    if (pr.loi == "UN") for (int i = 0; i < ncl; i++) dens[i] = 1.0 / (pr.maxi - pr.mini);
    if (pr.loi == "LU") for (int i = 0; i < ncl; i++) dens[i] = 1.0 / (pr.maxi - pr.mini) / x[i];
    if (pr.loi == "NO") for (int i = 0; i < ncl; i++) dens[i] = exp(-(x[i] - pr.mean) * (x[i] - pr.mean) / 2.0 / (pr.sdshape * pr.sdshape)) / (pr.sdshape * co);
    if (pr.loi == "LN") for (int i = 0; i < ncl; i++) dens[i] = exp(-(log(x[i]) / pr.mean) * (log(x[i]) / pr.mean) / 2.0 / (pr.sdshape * pr.sdshape)) / (pr.sdshape * co * x[i]);
    if (pr.loi == "GA")
        for (int i = 0; i < ncl; i++)
        {
            xb = x[i] * pr.sdshape / pr.mean;
            dens[i] = exp((pr.sdshape - 1.0) * log(xb) - xb);
        }
    som = 0.0;
    for (int i = 0; i < ncl; i++) som += dens[i];
    for (int i = 0; i < ncl; i++) dens[i] /= som;
    return dens;
}

/**
 * calcule un histogramme à partir de la loi fournie dans le prior
 */
long double* calhistexact(int ncl, long double* x, PriorC pr)
{
    long double *dens, xb, som, z, pas = x[1] - x[0];
    int ncl2 = 11;
    dens = new long double[ncl];
    if (pr.loi == "UN") for (int i = 0; i < ncl; i++) dens[i] = 1.0 / (pr.maxi - pr.mini);
    if (pr.loi == "LU")
        for (int i = 0; i < ncl; i++)
        {
            dens[i] = 0.0;
            for (int j = 0; j < ncl2 - 1; j++)
            {
                z = x[i] + (long double)j * pas / (long double)(ncl2 - 1);
                dens[i] += 1.0 / (pr.maxi - pr.mini) / z / (ncl2 - 1);
            }
        }
    if (pr.loi == "NO")
        for (int i = 0; i < ncl; i++)
        {
            dens[i] = 0.0;
            for (int j = 0; j < ncl2 - 1; j++)
            {
                z = x[i] + (long double)j * pas / (long double)(ncl2 - 1);
                dens[i] += exp(-(z - pr.mean) * (z - pr.mean) / 2.0 / (pr.sdshape * pr.sdshape)) / (pr.sdshape * co);
            }
        }
    if (pr.loi == "LN")
        for (int i = 0; i < ncl; i++)
        {
            dens[i] = 0.0;
            for (int j = 0; j < ncl2 - 1; j++)
            {
                z = x[i] + (long double)j * pas / (long double)(ncl2 - 1);
                dens[i] += exp(-(log(z) / pr.mean) * (log(z) / pr.mean) / 2.0 / (pr.sdshape * pr.sdshape)) / (pr.sdshape * co * x[i]);
            }
        }
    if (pr.loi == "GA")
        for (int i = 0; i < ncl; i++)
        {
            dens[i] = 0.0;
            for (int j = 0; j < ncl2 - 1; j++)
            {
                z = x[i] + (long double)j * pas / (long double)(ncl2 - 1);
                xb = z * pr.sdshape / pr.mean;
                dens[i] += exp((pr.sdshape - 1.0) * log(xb) - xb);
            }
        }
    som = 0.0;
    for (int i = 0; i < ncl; i++) som += dens[i];
    for (int i = 0; i < ncl; i++) dens[i] /= som;
    return dens;
}


/**
 * calcule la densité à partir d'un échantillon de valeurs'
 */
long double* calculdensite(int ncl, int n, long double* x, long double** y, int j, bool multithread)
{
    //cout<<"calculdensite\n";
    long double bw, *dens, sd, *z, denom, som;
    dens = new long double[ncl];
    z = new long double[n];
    for (int i = 0; i < n; i++) z[i] = y[i][j];
    sd = cal_sdL(n, z);
    if (sd > 0.0) bw = coefbw * exp(-0.2 * log((double)n)) * sd;
    else bw = 1.0;
    //cout<<"calculdensite  sd="<<sd<<"    bw="<<bw<<"\n";
#pragma omp parallel for shared(dens,z,x,bw) private(denom) if(multithread)
    for (int j = 0; j < ncl; j++)
    {
        dens[j] = 0.0;
        for (int i = 0; i < n; i++) dens[j] += exp(lnormal_dens(z[i], x[j], bw));
        //if(2*j==ncl-1) cout<<" avant denom dens["<<j<<"] = "<< dens[j]<<"\n";
        denom = pnorm5((x[ncl - 1] - x[j]) / bw, 0.0, 1.0) - pnorm5((x[0] - x[j]) / bw, 0.0, 1.0);
        if (denom > 0.0) dens[j] /= denom;
        //if(2*j==ncl-1) cout<<" apres denom dens["<<j<<"] = "<< dens[j]<<"   (denom="<<denom<<")\n";
    }
    som = 0.0;
    for (int j = 0; j < ncl; j++)
    {
        som += dens[j];
    }
    if (som > 0.0) for (int j = 0; j < ncl; j++) dens[j] /= som;
    delete []z;
    return dens;
}

/**
 * calcule un histogramme à partir d'un échantillon de valeurs'
 */
long double* calculhisto(int ncl, int n, long double* x, long double** y, int j, bool multithread)
{
    int k;
    long double som, *dens, *z, d, demipas = 0.5 * (x[1] - x[0]);
    dens = new long double[ncl];
    z = new long double[n];
    for (int i = 0; i < n; i++) z[i] = y[i][j];
    //for (int i=0;i<n;i++) cout <<"  "<<z[i];cout<<"\n";
    for (int i = 0; i < ncl; i++) dens[i] = 0.0;
    d = 1.0 / (long double)n;
#pragma omp parallel for shared(dens,z,x) private(k) if(multithread)
    for (int i = 0; i < n; i++)
    {
        for (k = 0; k < ncl; k++) if (fabs(z[i] - x[k]) <= demipas) break;
        dens[k] += d;
        if ((k == 0)or (k == ncl - 1)) dens[k] += d;
    }
    som = 0.0;
    for (int k = 0; k < ncl; k++)
    {
        som += dens[k];
    }
    if (som > 0.0) for (int k = 0; k < ncl; k++) dens[k] /= som;
    delete []z;
    return dens;
}

/**
 * traitement global du calcul de la densité des paramètres variables originaux
 * si le parametre est à valeurs entières avec moins de 30 classes, la densité est remplacée par un histogramme
 * sinon la densité est évaluée pour 1000 points du min au max
 */
void histodensO(int n, bool multithread, string progressfilename, int* iprog, int* nprog, long double** phistar)
{
    bool condition;
    long double *densprior, *denspost, *x;
    int ncl, ii;
    pardens = new pardensC[nparamcom];
    //cout<<"dans histodens npar="<<npar<<"\n";
    for (int j = 0; j < nparamcom; j++)
    {
        pardens[j].ncl = 1001;
        condition = false;
        if (j < npar)
        { //histparam
            if (rt.histparam[rt.scenchoisi[0] - 1][numpar[0][j]].category < 2)
            {
                if (rt.histparam[rt.scenchoisi[0] - 1][numpar[0][j]].prior.maxi - rt.histparam[rt.scenchoisi[0] - 1][numpar[0][j]].prior.mini <= 30)
                    pardens[j].ncl = rt.histparam[rt.scenchoisi[0] - 1][numpar[0][j]].prior.maxi - rt.histparam[rt.scenchoisi[0] - 1][numpar[0][j]].prior.mini + 1;
            }
            pardens[j].xmin = rt.histparam[rt.scenchoisi[0] - 1][numpar[0][j]].prior.mini;
            pardens[j].xmax = rt.histparam[rt.scenchoisi[0] - 1][numpar[0][j]].prior.maxi;
            if (header.nconditions > 0)
            {
                for (int i = 0; i < header.nconditions; i++)
                {
                    condition = ((rt.histparam[rt.scenchoisi[0] - 1][numpar[0][j]].name == header.condition[i].param1)or (rt.histparam[rt.scenchoisi[0] - 1][numpar[0][j]].name == header.condition[i].param2));
                    if (condition) break;
                }
            }
        }
        else if (j < nparamcom)
        { //mutparam
            pardens[j].xmin = rt.mutparam[j - npar].prior.mini;
            pardens[j].xmax = rt.mutparam[j - npar].prior.maxi;
        }
        pardens[j].xdelta = (pardens[j].xmax - pardens[j].xmin) / (long double)(pardens[j].ncl - 1);
        //cout<<nomparam[j]<<"   xmin="<<pardens[j].xmin<<"   xmax="<<pardens[j].xmax<<"   xdelta="<<pardens[j].xdelta<<"   ncl="<<pardens[j].ncl<<"\n";
        pardens[j].x = new long double[pardens[j].ncl];
        for (int k = 0; k < pardens[j].ncl; k++) pardens[j].x[k] = pardens[j].xmin + k * pardens[j].xdelta;
        //cout<<"1\n";
        if (pardens[j].ncl > 31)
        {
            if ((condition)or (j >= nparamcom))
            {
                pardens[j].priord = calculdensite(pardens[j].ncl, nsimpar, pardens[j].x, simpar, j, multithread);
            }
            else
            {
                if (j < npar) pardens[j].priord = caldensexact(pardens[j].ncl, pardens[j].x, rt.histparam[rt.scenchoisi[0] - 1][numpar[0][j]].prior);
                else pardens[j].priord = caldensexact(pardens[j].ncl, pardens[j].x, rt.mutparam[j - npar].prior);
            }
            pardens[j].postd = calculdensite(pardens[j].ncl, n, pardens[j].x, phistar, j, multithread);
        }
        else
        {
            /*if ((condition)or(j>=npar)) pardens[j].priord = calculhisto(pardens[j].ncl,nsimpar,pardens[j].x,simpar,j,multithread);
              else pardens[j].priord =calhistexact(pardens[j].ncl,pardens[j].x,rt.histparam[rt.scenchoisi[0]-1][numpar[0][j]].prior);
              pardens[j].postd = calculhisto(pardens[j].ncl,n,pardens[j].x,phistar,j,multithread);*/
            if ((condition)or (j >= nparamcom)) densprior = calculhisto(pardens[j].ncl, nsimpar, pardens[j].x, simpar, j, multithread);
            else
            {
                if (j < npar) densprior = calhistexact(pardens[j].ncl, pardens[j].x, rt.histparam[rt.scenchoisi[0] - 1][numpar[0][j]].prior);
                else densprior = calhistexact(pardens[j].ncl, pardens[j].x, rt.mutparam[j - npar].prior);
            }
            denspost = calculhisto(pardens[j].ncl, n, pardens[j].x, phistar, j, multithread);
            ncl = pardens[j].ncl;
            //for (int k=0;k<ncl;k++) cout<<"   "<<pardens[j].x[k];cout<<"\n";
            //for (int k=0;k<ncl;k++) cout<<"   "<<densprior[k];cout<<"\n";
            //for (int k=0;k<ncl;k++) cout<<"   "<<denspost[k];cout<<"\n";
            x = new long double[ncl];
            for (int k = 0; k < ncl; k++) x[k] = pardens[j].x[k];
            delete []pardens[j].x;
            pardens[j].ncl = 1001;
            pardens[j].x = new long double[pardens[j].ncl];
            pardens[j].priord = new long double[pardens[j].ncl];
            pardens[j].postd = new long double[pardens[j].ncl];
            pardens[j].xdelta = (x[ncl - 1] - x[0]) / (long double)(pardens[j].ncl - 1);
            //delta = x[1]-x[0];
            for (int k = 0; k < pardens[j].ncl; k++)
            {
                pardens[j].x[k] = x[0] + k * pardens[j].xdelta;
                for (ii = 0; ii < ncl - 1; ii++)
                {
                    if ((pardens[j].x[k] >= x[ii])and (pardens[j].x[k] <= x[ii + 1])) break;
                }
                pardens[j].priord[k] = 0.5 * (densprior[ii] + densprior[ii + 1]);
                pardens[j].postd[k] = 0.5 * (denspost[ii] + denspost[ii + 1]);
                //cout<<"x="<<pardens[j].x[k]<<"   priord="<<pardens[j].priord[k]<<"   postd="<<pardens[j].postd[k]<<"   ii="<<ii<<"\n";
            }
            delete [] densprior;
            delete [] denspost;
        }
        cout << "fin du calcul du parametre " << j + 1 << "  sur " << nparamcom << "\n";
        *iprog += 10;
        fprog.open(progressfilename.c_str());
        fprog << *iprog << "   " << *nprog << "\n";
        fprog.close();
    }
}

/**
 * traitement global du calcul de la densité des paramètres composites
 * si le parametre est à valeurs entières avec moins de 30 classes, la densité est remplacée par un histogramme
 * sinon la densité est évaluée pour 1000 points du min au max
 */
void histodensC(int n, bool multithread, string progressfilename, int* iprog, int* nprog, long double** phistarcompo)
{
    //bool condition;
    long double *densprior, *denspost, *x;
    int ncl, ii;
    pardenscompo = new pardensC[nparcompo];
    //cout<<"dans histodens npar="<<npar<<"\n";
    for (int j = 0; j < nparcompo; j++)
    {
        pardenscompo[j].ncl = 1001;
        //condition=false;
        pardenscompo[j].xmin = 1E100;
        pardenscompo[j].xmax = 0;
        for (int i = 0; i < nsimpar; i++)
        {
            if (pardenscompo[j].xmin > simparcompo[i][j]) pardenscompo[j].xmin = simparcompo[i][j];
            if (pardenscompo[j].xmax < simparcompo[i][j]) pardenscompo[j].xmax = simparcompo[i][j];
        }
        for (int i = 0; i < n; i++)
        {
            if (pardenscompo[j].xmin > phistarcompo[i][j]) pardenscompo[j].xmin = phistarcompo[i][j];
            if (pardenscompo[j].xmax < phistarcompo[i][j]) pardenscompo[j].xmax = phistarcompo[i][j];
        }
        pardenscompo[j].xdelta = (pardenscompo[j].xmax - pardenscompo[j].xmin) / (long double)(pardenscompo[j].ncl - 1);
        //cout<<nomparam[j]<<"   xmin="<<pardenscompo[j].xmin<<"   xmax="<<pardenscompo[j].xmax<<"   xdelta="<<pardenscompo[j].xdelta<<"   ncl="<<pardenscompo[j].ncl<<"\n";
        pardenscompo[j].x = new long double[pardenscompo[j].ncl];
        for (int k = 0; k < pardenscompo[j].ncl; k++) pardenscompo[j].x[k] = pardenscompo[j].xmin + k * pardenscompo[j].xdelta;
        //cout<<"1\n";
        if (pardenscompo[j].ncl > 31)
        {
            pardenscompo[j].priord = calculdensite(pardenscompo[j].ncl, nsimpar, pardenscompo[j].x, simparcompo, j, multithread);
            pardenscompo[j].postd = calculdensite(pardenscompo[j].ncl, n, pardenscompo[j].x, phistarcompo, j, multithread);
        }
        else
        {
            densprior = calculhisto(pardenscompo[j].ncl, nsimpar, pardenscompo[j].x, simparcompo, j, multithread);
            denspost = calculhisto(pardenscompo[j].ncl, n, pardenscompo[j].x, phistarcompo, j, multithread);
            ncl = pardenscompo[j].ncl;
            x = new long double[ncl];
            for (int k = 0; k < ncl; k++) x[k] = pardenscompo[j].x[k];
            delete []pardenscompo[j].x;
            pardenscompo[j].ncl = 1001;
            pardenscompo[j].x = new long double[pardenscompo[j].ncl];
            pardenscompo[j].priord = new long double[pardenscompo[j].ncl];
            pardenscompo[j].postd = new long double[pardenscompo[j].ncl];
            pardenscompo[j].xdelta = (x[ncl - 1] - x[0]) / (long double)(pardenscompo[j].ncl - 1);
            //delta = x[1]-x[0];
            for (int k = 0; k < pardenscompo[j].ncl; k++)
            {
                pardenscompo[j].x[k] = x[0] + k * pardenscompo[j].xdelta;
                for (ii = 0; ii < ncl - 1; ii++)
                {
                    if ((pardenscompo[j].x[k] >= x[ii])and (pardenscompo[j].x[k] <= x[ii + 1])) break;
                }
                pardenscompo[j].priord[k] = 0.5 * (densprior[ii] + densprior[ii + 1]);
                pardenscompo[j].postd[k] = 0.5 * (denspost[ii] + denspost[ii + 1]);
                //cout<<"x="<<pardenscompo[j].x[k]<<"   priord="<<pardenscompo[j].priord[k]<<"   postd="<<pardenscompo[j].postd[k]<<"   ii="<<ii<<"\n";
            }
            delete [] densprior;
            delete [] denspost;
        }
        cout << "fin du calcul du parametre " << j + 1 << "  sur " << nparcompo << "\n";
        *iprog += 10;
        ofstream fprog(progressfilename.c_str(), ios::out);
        fprog << *iprog << "   " << *nprog << "\n";
        fprog.close();
    }
}

/**
 * traitement global du calcul de la densité des paramètres scaled
 * si le parametre est à valeurs entières avec moins de 30 classes, la densité est remplacée par un histogramme
 * sinon la densité est évaluée pour 1000 points du min au max
 */
void histodensS(int n, bool multithread, string progressfilename, int* iprog, int* nprog, long double** phistarscaled)
{
    //bool condition;
    long double *densprior, *denspost, *x;
    int ncl, ii;
    pardensscaled = new pardensC[nparscaled];
    //cout<<"dans histodens npar="<<npar<<"\n";
    for (int j = 0; j < nparscaled; j++)
    {
        pardensscaled[j].ncl = 1001;
        //condition=false;
        pardensscaled[j].xmin = 1E100;
        pardensscaled[j].xmax = 0;
        for (int i = 0; i < nsimpar; i++)
        {
            if (pardensscaled[j].xmin > simparscaled[i][j]) pardensscaled[j].xmin = simparscaled[i][j];
            if (pardensscaled[j].xmax < simparscaled[i][j]) pardensscaled[j].xmax = simparscaled[i][j];
        }
        for (int i = 0; i < n; i++)
        {
            if (pardensscaled[j].xmin > phistarscaled[i][j]) pardensscaled[j].xmin = phistarscaled[i][j];
            if (pardensscaled[j].xmax < phistarscaled[i][j]) pardensscaled[j].xmax = phistarscaled[i][j];
        }
        pardensscaled[j].xdelta = (pardensscaled[j].xmax - pardensscaled[j].xmin) / (long double)(pardensscaled[j].ncl - 1);
        cout << nomparamS[j] << "   xmin=" << pardensscaled[j].xmin << "   xmax=" << pardensscaled[j].xmax << "   xdelta=" << pardensscaled[j].xdelta << "   ncl=" << pardensscaled[j].ncl << "\n";
        pardensscaled[j].x = new long double[pardensscaled[j].ncl];
        for (int k = 0; k < pardensscaled[j].ncl; k++) pardensscaled[j].x[k] = pardensscaled[j].xmin + k * pardensscaled[j].xdelta;
        //cout<<"1\n";
        if (pardensscaled[j].ncl > 31)
        {
            pardensscaled[j].priord = calculdensite(pardensscaled[j].ncl, nsimpar, pardensscaled[j].x, simparscaled, j, multithread);
            pardensscaled[j].postd = calculdensite(pardensscaled[j].ncl, n, pardensscaled[j].x, phistarscaled, j, multithread);
        }
        else
        {
            densprior = calculhisto(pardensscaled[j].ncl, nsimpar, pardensscaled[j].x, simparscaled, j, multithread);
            denspost = calculhisto(pardensscaled[j].ncl, n, pardensscaled[j].x, phistarscaled, j, multithread);
            ncl = pardensscaled[j].ncl;
            x = new long double[ncl];
            for (int k = 0; k < ncl; k++) x[k] = pardensscaled[j].x[k];
            delete []pardensscaled[j].x;
            pardensscaled[j].ncl = 1001;
            pardensscaled[j].x = new long double[pardensscaled[j].ncl];
            pardensscaled[j].priord = new long double[pardensscaled[j].ncl];
            pardensscaled[j].postd = new long double[pardensscaled[j].ncl];
            pardensscaled[j].xdelta = (x[ncl - 1] - x[0]) / (long double)(pardensscaled[j].ncl - 1);
            //delta = x[1]-x[0];
            for (int k = 0; k < pardensscaled[j].ncl; k++)
            {
                pardensscaled[j].x[k] = x[0] + k * pardensscaled[j].xdelta;
                for (ii = 0; ii < ncl - 1; ii++)
                {
                    if ((pardensscaled[j].x[k] >= x[ii])and (pardensscaled[j].x[k] <= x[ii + 1])) break;
                }
                pardensscaled[j].priord[k] = 0.5 * (densprior[ii] + densprior[ii + 1]);
                pardensscaled[j].postd[k] = 0.5 * (denspost[ii] + denspost[ii + 1]);
                //cout<<"x="<<pardensscaled[j].x[k]<<"   priord="<<pardensscaled[j].priord[k]<<"   postd="<<pardensscaled[j].postd[k]<<"   ii="<<ii<<"\n";
            }
            delete [] densprior;
            delete [] denspost;
        }
        cout << "fin du calcul du parametre " << j + 1 << "  sur " << nparscaled << "\n";
        *iprog += 10;
        fprog.open(progressfilename.c_str());
        fprog << *iprog << "   " << *nprog << "\n";
        fprog.close();
    }
}

/**
 *calcule les statistiques de paramètres quelconques
 */
parstatC* calparstat(int n, int nparam, long double** par)
{
    long double* x;
    parstatC* parst;
    parst = new parstatC[nparam];
    x = new long double[n];
    for (int j = 0; j < nparam; j++)
    {
        //cout<<"calcul des stat du parametre "<<j<<" sur "<<nparam<<"\n";
        for (int i = 0; i < n; i++) x[i] = par[i][j];
        sort(&x[0], &x[n]);
        //cout<<"apres le sort\n";
        parst[j].q025 = x[(int)floor(0.025 * n + 0.5) - 1];
        parst[j].q050 = x[(int)floor(0.050 * n + 0.5) - 1];
        parst[j].q250 = x[(int)floor(0.250 * n + 0.5) - 1];
        parst[j].med = x[(int)floor(0.500 * n + 0.5) - 1];
        parst[j].q750 = x[(int)floor(0.750 * n + 0.5) - 1];
        parst[j].q950 = x[(int)floor(0.950 * n + 0.5) - 1];
        parst[j].q975 = x[(int)floor(0.975 * n + 0.5) - 1];
        //cout<<"apres les quantiles\n";
        parst[j].moy = cal_moyL(n, x);
        //cout<<"apres la moyenne\n";
        parst[j].mod = cal_modeL(n, x);
        //cout<<"apres le mode\n";
    }
    delete []x;
    return parst;
}


/**
 *sauvegarde les statistiques et les densités des paramètres
 */
void saveparstat(int nsel, string path, string ident)
{
    time_t rawtime;
    struct tm* timeinfo;
    time(&rawtime);
    timeinfo = localtime(&rawtime);
    string nomfiparstat;
    ofstream f1;
    if (original)
    {
        nomfiparstat = path + ident + "_paramstatdens_original.txt";
        cout << nomfiparstat << "\n";
        f1.open(nomfiparstat.c_str());
        for (int j = 0; j < nparamcom; j++)
        {
            f1 << nomparamO[j] << "\n";//cout<<"1\n";
            f1 << setiosflags(ios::scientific) << setw(11) << setprecision(5) << parstat[j].moy << "  " << parstat[j].med << "  ";
            f1 << parstat[j].mod << "  " << parstat[j].q025 << "  " << parstat[j].q050 << "  " << parstat[j].q250 << "  ";
            f1 << parstat[j].q750 << "  " << parstat[j].q950 << "  " << parstat[j].q975 << "\n";
            //f1<<setiosflags(ios::fixed)<<pardens[j].ncl<<"\n  ";
            //f1<<setiosflags(ios::scientific)<<setw(11)<<setprecision(5);
            f1 << pardens[j].ncl << "\n";
            for (int i = 0; i < pardens[j].ncl; i++) f1 << "  " << pardens[j].x[i];
            f1 << "\n";
            for (int i = 0; i < pardens[j].ncl; i++) f1 << "  " << pardens[j].priord[i];
            f1 << "\n";
            for (int i = 0; i < pardens[j].ncl; i++) f1 << "  " << pardens[j].postd[i];
            f1 << "\n";
        }
        f1.close();
        nomfiparstat = path + ident + "_mmmq_original.txt";
        cout << nomfiparstat << "\n";
        f1.open(nomfiparstat.c_str());
        f1 << "DIYABC :                      ABC parameter estimation                         " << asctime(timeinfo) << "\n";
        f1 << "Data file       : " << header.datafilename << "\n";
        f1 << "Reference table : " << rt.filename << "\n";
        switch (numtransf)
        {
        case 1: f1 << "No transformation of parameters\n";
            break;
        case 2: f1 << "Transformation LOG of parameters\n";
            break;
        case 3: f1 << "Transformation LOGIT of parameters\n";
            break;
        case 4: f1 << "Transformation LOG(TG) of parameters\n";
            break;
        }
        f1 << "Chosen scenario(s) : ";
        for (int i = 0; i < rt.nscenchoisi; i++) f1 << rt.scenchoisi[i];
        f1 << "\n";
        f1 << "Number of simulated data sets : " << rt.nreclus << "\n";
        f1 << "Number of selected data sets  : " << nsel << "\n\n";
        f1 << "Parameter          mean     median    mode      q025      q050      q250      q750      q950      q975\n";
        f1 << "------------------------------------------------------------------------------------------------------\n";
        f1 << setiosflags(ios::scientific) << setprecision(2);
        for (int j = 0; j < nparamcom; j++)
        {
            f1 << nomparamO[j];
            for (int i = 0; i < 17 - (int)(nomparamO[j].length()); i++) f1 << " ";
            f1 << parstat[j].moy << "  " << parstat[j].med << "  " << parstat[j].mod << "  " << parstat[j].q025 << "  " << parstat[j].q050;
            f1 << "  " << parstat[j].q250 << "  " << parstat[j].q750 << "  " << parstat[j].q950 << "  " << parstat[j].q975 << "\n";
        }
        f1.close();
    }
    if (composite)
    {
        nomfiparstat = path + ident + "_paramstatdens_composite.txt";
        cout << nomfiparstat << "\n";
        f1.open(nomfiparstat.c_str());
        for (int j = 0; j < nparcompo; j++)
        {
            f1 << nomparamC[j] << "\n";//cout<<"1\n";
            f1 << setiosflags(ios::scientific) << setw(11) << setprecision(5) << parstatcompo[j].moy << "  " << parstatcompo[j].med << "  ";
            f1 << parstatcompo[j].mod << "  " << parstatcompo[j].q025 << "  " << parstatcompo[j].q050 << "  " << parstatcompo[j].q250 << "  ";
            f1 << parstatcompo[j].q750 << "  " << parstatcompo[j].q950 << "  " << parstatcompo[j].q975 << "\n";
            f1 << pardenscompo[j].ncl << "\n";
            for (int i = 0; i < pardenscompo[j].ncl; i++) f1 << "  " << pardenscompo[j].x[i];
            f1 << "\n";
            for (int i = 0; i < pardenscompo[j].ncl; i++) f1 << "  " << pardenscompo[j].priord[i];
            f1 << "\n";
            for (int i = 0; i < pardenscompo[j].ncl; i++) f1 << "  " << pardenscompo[j].postd[i];
            f1 << "\n";
        }
        f1.close();
        nomfiparstat = path + ident + "_mmmq_composite.txt";
        cout << nomfiparstat << "\n";
        f1.open(nomfiparstat.c_str());
        f1 << "DIYABC :                      ABC parameter estimation                         " << asctime(timeinfo) << "\n";
        f1 << "Data file       : " << header.datafilename << "\n";
        f1 << "Reference table : " << rt.filename << "\n";
        switch (numtransf)
        {
        case 1: f1 << "No transformation of parameters\n";
            break;
        case 2: f1 << "Transformation LOG of parameters\n";
            break;
        case 3: f1 << "Transformation LOGIT of parameters\n";
            break;
        case 4: f1 << "Transformation LOG(TG) of parameters\n";
            break;
        }
        f1 << "Chosen scenario(s) : ";
        for (int i = 0; i < rt.nscenchoisi; i++) f1 << rt.scenchoisi[i];
        f1 << "\n";
        f1 << "Number of simulated data sets : " << rt.nreclus << "\n";
        f1 << "Number of selected data sets  : " << nsel << "\n\n";
        f1 << "Parameter          mean     median    mode      q025      q050      q250      q750      q950      q975\n";
        f1 << "------------------------------------------------------------------------------------------------------\n";
        f1 << setiosflags(ios::scientific) << setprecision(2);
        for (int j = 0; j < nparcompo; j++)
        {
            f1 << nomparamC[j];
            for (int i = 0; i < 17 - (int)(nomparamC[j].length()); i++) f1 << " ";
            f1 << parstatcompo[j].moy << "  " << parstatcompo[j].med << "  " << parstatcompo[j].mod << "  " << parstatcompo[j].q025 << "  " << parstatcompo[j].q050;
            f1 << "  " << parstatcompo[j].q250 << "  " << parstatcompo[j].q750 << "  " << parstatcompo[j].q950 << "  " << parstatcompo[j].q975 << "\n";
        }
        f1.close();
    }
    if (scaled)
    {
        nomfiparstat = path + ident + "_paramstatdens_scaled.txt";
        cout << nomfiparstat << "\n";
        f1.open(nomfiparstat.c_str());
        for (int j = 0; j < nparscaled; j++)
        {
            f1 << nomparamS[j] << "\n";//cout<<"1\n";
            f1 << setiosflags(ios::scientific) << setw(11) << setprecision(5) << parstatscaled[j].moy << "  " << parstatscaled[j].med << "  ";
            f1 << parstatscaled[j].mod << "  " << parstatscaled[j].q025 << "  " << parstatscaled[j].q050 << "  " << parstatscaled[j].q250 << "  ";
            f1 << parstatscaled[j].q750 << "  " << parstatscaled[j].q950 << "  " << parstatscaled[j].q975 << "\n";
            f1 << pardensscaled[j].ncl << "\n";
            for (int i = 0; i < pardensscaled[j].ncl; i++) f1 << "  " << pardensscaled[j].x[i];
            f1 << "\n";
            for (int i = 0; i < pardensscaled[j].ncl; i++) f1 << "  " << pardensscaled[j].priord[i];
            f1 << "\n";
            for (int i = 0; i < pardensscaled[j].ncl; i++) f1 << "  " << pardensscaled[j].postd[i];
            f1 << "\n";
        }
        f1.close();
        nomfiparstat = path + ident + "_mmmq_scaled.txt";
        cout << nomfiparstat << "\n";
        f1.open(nomfiparstat.c_str());
        f1 << "DIYABC :                      ABC parameter estimation                         " << asctime(timeinfo) << "\n";
        f1 << "Data file       : " << header.datafilename << "\n";
        f1 << "Reference table : " << rt.filename << "\n";
        switch (numtransf)
        {
        case 1: f1 << "No transformation of parameters\n";
            break;
        case 2: f1 << "Transformation LOG of parameters\n";
            break;
        case 3: f1 << "Transformation LOGIT of parameters\n";
            break;
        case 4: f1 << "Transformation LOG(TG) of parameters\n";
            break;
        }
        f1 << "Chosen scenario(s) : ";
        for (int i = 0; i < rt.nscenchoisi; i++) f1 << rt.scenchoisi[i];
        f1 << "\n";
        f1 << "Number of simulated data sets : " << rt.nreclus << "\n";
        f1 << "Number of selected data sets  : " << nsel << "\n\n";
        f1 << "Parameter          mean     median    mode      q025      q050      q250      q750      q950      q975\n";
        f1 << "------------------------------------------------------------------------------------------------------\n";
        f1 << setiosflags(ios::scientific) << setprecision(2);
        for (int j = 0; j < nparscaled; j++)
        {
            f1 << nomparamS[j];
            for (int i = 0; i < 17 - (int)(nomparamS[j].length()); i++) f1 << " ";
            f1 << parstatscaled[j].moy << "  " << parstatscaled[j].med << "  " << parstatscaled[j].mod << "  " << parstatscaled[j].q025 << "  " << parstatscaled[j].q050;
            f1 << parstatscaled[j].q250 << "  " << parstatscaled[j].q750 << "  " << parstatscaled[j].q950 << "  " << parstatscaled[j].q975 << "\n";
        }
        f1.close();
    }
}

/**
 * Interprète les paramètres de la ligne de commande et 
 * effectue l'estimation ABC des paramètres (directe + régression locale)
 */
void doestim(string opt)
{
    int nstatOK, iprog, nprog;
    int nrec, nsel = 0, ns, nrecpos = 0;
    string s, s0, s1;
    vector<string> ss, ss1;
    float* stat_obs;
    long double** matC;
    long double **phistar, **phistarcompo, **phistarscaled;
    progressfilename = path + ident + "_progress.txt";
    cout << "debut doestim  options : " << opt << "\n";
    splitwords(opt, ";", ss);
    ns = ss.size();
    for (int i = 0; i < ns; i++)
    { //cout<<ss[i]<<"\n";
        s0 = ss[i].substr(0, 2);
        s1 = ss[i].substr(2);
        if (s0 == "s:")
        {
            splitwords(s1, ",", ss1);
            rt.nscenchoisi = ss1.size();
            rt.scenchoisi = new int[rt.nscenchoisi];
            for (int j = 0; j < rt.nscenchoisi; j++) rt.scenchoisi[j] = atoi(ss1[j].c_str());
            nrecpos = 0;
            for (int j = 0; j < rt.nscenchoisi; j++) nrecpos += rt.nrecscen[rt.scenchoisi[j] - 1];
            cout << "scenario(s) choisi(s) : ";
            for (int j = 0; j < rt.nscenchoisi; j++)
            {
                cout << rt.scenchoisi[j];
                if (j < rt.nscenchoisi - 1) cout << ",";
            }
            cout << "\n";
        }
        else if (s0 == "n:")
        {
            nrec = atoi(s1.c_str());
            if (nrec > nrecpos) nrec = nrecpos;
            cout << "nombre total de jeux de données considérés (pour le(s) scénario(s) choisi(s) )= " << nrec << "\n";
        }
        else if (s0 == "m:")
        {
            nsel = atoi(s1.c_str());
            cout << "nombre de jeux de données considérés pour la régression locale = " << nsel << "\n";
        }
        else if (s0 == "t:")
        {
            numtransf = atoi(s1.c_str());
            switch (numtransf)
            {
            case 1: cout << " pas de transformation des paramètres\n";
                break;
            case 2: cout << " transformation log des paramètres\n";
                break;
            case 3: cout << " transformation logit des paramètres\n";
                break;
            case 4: cout << " transformation log(tg) des paramètres\n";
                break;
            }
        }
        else if (s0 == "p:")
        {
            original = (s1.find("o") != string::npos);
            composite = (s1.find("c") != string::npos);
            scaled = (s1.find("s") != string::npos);
            if ((original)and (not composite)and (not scaled)) cout << "paramètres  originaux  ";
            if ((not original)and (composite)and (not scaled)) cout << "paramètres  composites  ";
            if ((not original)and (not composite)and (scaled)) cout << "paramètres  scaled  ";
            if ((original)and (composite)and (not scaled)) cout << "paramètres  originaux et composites ";
            if ((original)and (not composite)and (scaled)) cout << "paramètres  originaux et scaled ";
            if ((not original)and (composite)and (scaled)) cout << "paramètres  composites et scaled ";
            if ((original)and (composite)and (scaled)) cout << "paramètres  originaux, composites et scaled ";
            cout << "\n";
        }
    }

    nstatOK = rt.cal_varstat();
    cout << "apres cal_varstat   nstatOK = " << nstatOK << "\n";
    //header.calstatobs(statobsfilename);  cout<<"apres read_statobs\n";
    stat_obs = new float[header.nstat];
    for (int i = 0; i < header.nstat; i++) stat_obs[i] = header.stat_obs[i];
    nprog = 200;
    iprog = 2;
    fprog.open(progressfilename.c_str());
    fprog << iprog << "   " << nprog << "\n";
    fprog.close();
    cout << "--->" << iprog << "   sur " << nprog << "\n";
    rt.alloue_enrsel(nsel);
    rt.cal_dist(nrec, nsel, stat_obs, false, false); //cout<<"apres cal_dist\n";
    iprog += 8;
    fprog.open(progressfilename.c_str());
    fprog << iprog << "   " << nprog << "\n";
    fprog.close();
    cout << "--->" << iprog << "   sur " << nprog << "\n";
    det_numpar();
    cout << "apres det_numpar\n";
    nprog = 12 + 10 * (nparamcom + nparcompo + nparscaled);
    if (original) nprog += 6;
    if (composite) nprog += 6;
    if (scaled) nprog += 6;
    rempli_mat(nsel, stat_obs);
    cout << "apres rempli_mat\n";
    if (not deltanul) matC = cal_matC(nsel);
    cout << "cal_matC\n";
    cout << "original\n";
    if (original)
    {
        cout << "\ntraitement des parametres originaux\n";
        recalparamO(nsel);
        cout << "apres recalparamO\n";
        if (not deltanul)
        {
            rempli_parsim(nsel, nparamcom);
            cout << "apres rempli_parsim(O)\n";
            local_regression(nsel, nparamcom, matC);
            cout << "apres local_regression\n";
        }
        iprog += 2;
        fprog.open(progressfilename.c_str());
        fprog << iprog << "   " << nprog << "\n";
        fprog.close();
        cout << "original--->" << iprog << "   sur " << nprog << "\n";
        phistar = new long double*[nsel];
        for (int i = 0; i < nsel; i++) phistar[i] = new long double[nparamcom];
        if (not deltanul) calphistarO(nsel, phistar); else copyphistar(nsel, nparamcom, phistar);
        cout << "apres calphistarO\n";
    }
    if (composite)
    {
        cout << "\ntraitement des parametres composites\n";
        recalparamC(nsel);
        cout << "apres recalparamC\n";
        if (not deltanul)
        {
            rempli_parsim(nsel, nparcompo);
            cout << "apres rempli_parsim(C)\n";
            local_regression(nsel, nparcompo, matC);
            cout << "apres local_regression\n";
        }
        iprog += 2;
        fprog.open(progressfilename.c_str());
        fprog << iprog << "   " << nprog << "\n";
        fprog.close();
        cout << "composite--->" << iprog << "   sur " << nprog << "\n";
        phistarcompo = new long double*[nsel];
        for (int i = 0; i < nsel; i++) phistarcompo[i] = new long double[nparcompo];
        if (not deltanul) calphistarC(nsel, phistarcompo); else copyphistar(nsel, nparcompo, phistarcompo);
        cout << "apres calphistarcompo\n";
    }
    if (scaled)
    {
        cout << "\ntraitement des parametres scaled\n";
        recalparamS(nsel);
        cout << "apres recalparamS\n";
        if (not deltanul)
        {
            rempli_parsim(nsel, nparscaled);
            cout << "apres rempli_parsim(S)\n";
            local_regression(nsel, nparscaled, matC);
            cout << "apres local_regression\n";
        }
        iprog += 2;
        fprog.open(progressfilename.c_str());
        fprog << iprog << "   " << nprog << "\n";
        fprog.close();
        cout << "scaled--->" << iprog << "   sur " << nprog << "\n";
        phistarscaled = new long double*[nsel];
        for (int i = 0; i < nsel; i++) phistarscaled[i] = new long double[nparscaled];
        if (not deltanul) calphistarS(nsel, phistarscaled); else copyphistar(nsel, nparscaled, phistarscaled);
        cout << "apres calphistarscaled\n";
    }
    for (int i = 0; i < nsel; i++) delete []matX0[i];
    delete []matX0;
    cout << "apres delete matX0\n";
    det_nomparam();
    cout << "apres det_nomparam\n";
    savephistar(nsel, path, ident, phistar, phistarcompo, phistarscaled);
    cout << "apres savephistar\n";
    if (original)
    {
        lisimparO();
        cout << "apres lisimparO\n";
        iprog += 2;
        fprog.open(progressfilename.c_str());
        fprog << iprog << "   " << nprog << "\n";
        fprog.close();
        cout << "--->" << iprog << "   sur " << nprog << "\n";
    }
    if (composite)
    {
        lisimparC();
        cout << "apres lisimparC\n";
        iprog += 2;
        fprog.open(progressfilename.c_str());
        fprog << iprog << "   " << nprog << "\n";
        fprog.close();
        cout << "--->" << iprog << "   sur " << nprog << "\n";
    }
    if (scaled)
    {
        lisimparS();
        cout << "apres lisimparS\n";
        iprog += 2;
        fprog.open(progressfilename.c_str());
        fprog << iprog << "   " << nprog << "\n";
        fprog.close();
        cout << "--->" << iprog << "   sur " << nprog << "\n";
    }
    if (original)
    {
        histodensO(nsel, multithread, progressfilename, &iprog, &nprog, phistar);
        cout << "apres histodensO\n";
        cout << "--->" << iprog << "   sur " << nprog << "\n";
    }
    if (composite)
    {
        histodensC(nsel, multithread, progressfilename, &iprog, &nprog, phistarcompo);
        cout << "apres histodensC\n";
        cout << "--->" << iprog << "   sur " << nprog << "\n";
    }
    if (scaled)
    {
        histodensS(nsel, multithread, progressfilename, &iprog, &nprog, phistarscaled);
        cout << "apres histodensS\n";
        cout << "--->" << iprog << "   sur " << nprog << "\n";
    }
    if (original)
    {
        parstat = calparstat(nsel, nparamcom, phistar);
        cout << "apres calparstatO\n";
        iprog += 2;
        fprog.open(progressfilename.c_str());;
        fprog << iprog << "   " << nprog << "\n";
        fprog.close();
        cout << "--->" << iprog << "   sur " << nprog << "\n";
    }
    if (composite)
    {
        parstatcompo = calparstat(nsel, nparcompo, phistarcompo);
        cout << "apres calparstatC\n";
        iprog += 2;
        fprog.open(progressfilename.c_str());
        fprog << iprog << "   " << nprog << "\n";
        fprog.close();
        cout << "--->" << iprog << "   sur " << nprog << "\n";
    }
    if (scaled)
    {
        parstatscaled = calparstat(nsel, nparscaled, phistarscaled);
        cout << "apres calparstatS\n";
        iprog += 2;
        fprog.open(progressfilename.c_str());
        fprog << iprog << "   " << nprog << "\n";
        fprog.close();
        cout << "--->" << iprog << "   sur " << nprog << "\n";
    }
    saveparstat(nsel, path, ident);
    rt.desalloue_enrsel(nsel);
    iprog += 2;
    fprog.open(progressfilename.c_str());
    fprog << iprog << "   " << nprog << "\n";
    fprog.close();
    cout << "--->" << iprog << "   sur " << nprog << "\n";
}

