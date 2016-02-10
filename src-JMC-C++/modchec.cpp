/*
 * modchec.cpp
 *
 *  Created on: 4 april 2011
 *      Author: cornuet
 */

#include <string>
#include <vector>
#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <cmath>
#include <time.h>

#include "header.hpp"
#include "reftable.hpp"
#include "mesutils.hpp"
#include "particleset.hpp"
#include "acploc.hpp"
#include "estimparam.hpp"

using namespace std;

extern int debuglevel;
extern ParticleSetC ps;
extern HeaderC header;
extern ReftableC rt;
extern enregC* enreg;
extern int nenr, numtransf;
extern int nparamcom, nparcompo, nparscaled;
extern string progressfilename, statobsfilename, headerfilename, path, ident;
extern bool original, composite, scaled;
long double** phistarOK;
extern string* stat_type;
extern int* stat_num;
extern bool multithread;
extern string scurfile;
extern bool deltanul;
extern vector<LocusGroupC> groupe;
extern vector<ScenarioC> scenario;

extern ofstream fprog;

long double **ssphistar, **ssref;

bool resetstats(string s) {
    cout << "debut de resetstats\n";
    cout << s << "\n";
    int j, ns, nq, gr, lonentstat, *groupmodif, ngroupmodif, *nustat;
    vector<string> ss, qq, ss1;
    vector<vector<StatsnpC>> statsnp;
    vector<vector<string>> statname;
    StatsnpC stsnp;
    bool trouve;
    int catsnp, gg, gs, gt;
    statsnp.resize(header.ngroupes + 1);
    for (gr = 0; gr <= header.ngroupes; gr++) statsnp[gr].resize(0);
    //recopiage des anciennes stat
    statname.resize(header.ngroupes + 1);
    for (gr = 1; gr <= header.ngroupes; gr++) statname[gr].resize(groupe[gr].nstat);
    cout << "apres statname.resize\n";
    gs = 0;
    for (gr = 1; gr <= header.ngroupes; gr++) {
        for (int i = 0; i < groupe[gr].nstat; i++) {
            //cout<<"header.statname["<<gs<<"] = \n";
            //cout<<"                         "<< header.statname[gs]<<"\n";
            statname[gr][i] = header.statname[gs];
            gs++;
        }
    }
    nustat = new int[header.ngroupes + 1];
    for (gr = 1; gr <= header.ngroupes; gr++) nustat[gr] = 0;
    //repérage des groupes modifiés
    groupmodif = new int[header.ngroupes];
    for (gr = 0; gr < header.ngroupes; gr++) groupmodif[gr] = 0;
    ngroupmodif = 0;
    splitwords(s, " ", ss);
    ns = ss.size();
    for (int i = 0; i < ns; i++) {
        splitwords(ss[i], "_", qq);
        nq = qq.size();
        if (nq != 3) return false;
        gr = atoi(qq[1].c_str());
        if (ngroupmodif == 0) {
            groupmodif[ngroupmodif] = gr;
            ngroupmodif++;
        } else {
            gg = 0;
            while ((gg < ngroupmodif)and (gr != groupmodif[gg])) gg++;
            if (gg == ngroupmodif) {
                groupmodif[gg] = gr;
                ngroupmodif++;
            }
        }
    }
    //comptage des stat
    if ((ngroupmodif > 0)and (debuglevel == 2)) {
        cout << "ngroupmodif = " << ngroupmodif << "\n";
        for (gg = 0; gg < ngroupmodif; gg++) cout << "groupmodif[" << gg << "] = " << groupmodif[gg] << "\n";
    }
    for (gg = 0; gg < ngroupmodif; gg++) groupe[groupmodif[gg]].nstat = 0;
    //for (gr=1;gr<=header.ngroupes;gr++) {groupe[gr].nstat=0;} 
    cout << header.ngroupes << " groupe(s)\n";
    //delete [] ss;
    //ss =splitwords(s," ",&ns);
    if (debuglevel == 2) cout << "ns=" << ns << "\n";
    for (int i = 0; i < ns; i++) {
        splitwords(ss[i], "_", qq);
        nq = qq.size();
        if (nq != 3) return false;
        if (debuglevel == 2) cout << "nq=" << nq << "   qq[1]=" << qq[1] << "\n";
        gr = atoi(qq[1].c_str());
        if (debuglevel == 2) cout << "gr=" << gr << "\n";
        groupe[gr].nstat++;
    }
    header.nstat = 0;
    for (gr = 1; gr <= header.ngroupes; gr++) header.nstat += groupe[gr].nstat;
    if (debuglevel == 2) {
        cout << "dans resetstat nstat = " << header.nstat << "\n";
        for (gr = 1; gr <= header.ngroupes; gr++) cout << "      groupe " << gr << "   nstat=" << groupe[gr].nstat << "\n";
    }

    for (gg = 0; gg < ngroupmodif; gg++) groupe[groupmodif[gg]].sumstat = vector<StatC>(groupe[groupmodif[gg]].nstat);

    if (debuglevel == 2) {
        cout << "dans resetstat nstat = " << header.nstat << "\n";
        for (gr = 1; gr <= header.ngroupes; gr++) cout << "      groupe " << gr << "   nstat=" << groupe[gr].nstat << "\n";
    }

    //delete [] ss;
    //ss =splitwords(s," ",&ns);
    //delete []header.statname;
    lonentstat = header.entetestat.length();
    header.entetestat = "";
    header.statname = vector<string>(header.nstat);
    gs = 0;
    gt = 0;
    for (gr = 1; gr <= header.ngroupes; gr++) {
        gg = 0;
        while ((gg < ngroupmodif)and (gr != groupmodif[gg])) gg++;
        if (gg == ngroupmodif) {
            for (int i = 0; i < groupe[gr].nstat; i++) {
                header.statname[gs] = statname[gr][i];
                header.entetestat += centre(header.statname[gs], 14);
                gs++;
            }
        } else {
            for (int i = 0; i < groupe[gr].nstat; i++) {
                header.statname[gs] = ss[gt];
                gt++;
                header.entetestat += centre(header.statname[gs], 14);
                gs++;
            }
        }
    }
    header.entete.erase(header.entete.length() - lonentstat, lonentstat);
    header.entete += header.entetestat;
    cout << "\n" << header.entete << "\n";
    for (int i = 0; i < ns; i++) {
        if (debuglevel == 2) cout << "\nss[" << i << "] = " << ss[i] << "\n";
        splitwords(ss[i], "_", qq);
        nq = qq.size();
        gr = atoi(qq[1].c_str());
        j = 0;
        while (qq[0] != stat_type[j]) j++;
        if (debuglevel == 2) cout << "ss[" << i << "] = " << ss[i] << "   j=" << j << "\n";
        if (groupe[gr].type == 0) { //MICROSAT
            if (stat_num[j] < 5) {
                groupe[gr].sumstat[nustat[gr]].cat = stat_num[j];
                groupe[gr].sumstat[nustat[gr]].samp = atoi(qq[2].c_str());
                nustat[gr]++;
            } else if (stat_num[j] < 12) {
                groupe[gr].sumstat[nustat[gr]].cat = stat_num[j];
                splitwords(qq[2], "&", ss1);
                groupe[gr].sumstat[nustat[gr]].samp = atoi(ss1[0].c_str());
                groupe[gr].sumstat[nustat[gr]].samp1 = atoi(ss1[1].c_str());
                nustat[gr]++;
            } else if (stat_num[j] == 12) {
                groupe[gr].sumstat[nustat[gr]].cat = stat_num[j];
                splitwords(qq[2], "&", ss1);
                groupe[gr].sumstat[nustat[gr]].samp = atoi(ss1[0].c_str());
                groupe[gr].sumstat[nustat[gr]].samp1 = atoi(ss1[1].c_str());
                groupe[gr].sumstat[nustat[gr]].samp2 = atoi(ss1[2].c_str());
                nustat[gr]++;
            }
        } else if (groupe[gr].type == 1) { //DNA SEQUENCE
            if (debuglevel == 2) cout << "DNA SEQUENCE stat_num[" << j << "] = " << stat_num[j] << "   qq[2]=" << qq[2] << "\n";
            if (debuglevel == 2) cout << "gr=" << gr << "   nustat[gr]=" << nustat[gr] << "\n";
            if (stat_num[j] > -9) {
                groupe[gr].sumstat[nustat[gr]].cat = stat_num[j];
                groupe[gr].sumstat[nustat[gr]].samp = atoi(qq[2].c_str());
                nustat[gr]++;
            } else if (stat_num[j] > -14) {
                groupe[gr].sumstat[nustat[gr]].cat = stat_num[j];
                splitwords(qq[2], "&", ss1);
                groupe[gr].sumstat[nustat[gr]].samp = atoi(ss1[0].c_str());
                groupe[gr].sumstat[nustat[gr]].samp1 = atoi(ss1[1].c_str());
                nustat[gr]++;
            } else if (stat_num[j] == -14) {
                groupe[gr].sumstat[nustat[gr]].cat = stat_num[j];
                splitwords(qq[2], "&", ss1);
                groupe[gr].sumstat[nustat[gr]].samp = atoi(ss1[0].c_str());
                groupe[gr].sumstat[nustat[gr]].samp1 = atoi(ss1[1].c_str());
                groupe[gr].sumstat[nustat[gr]].samp2 = atoi(ss1[2].c_str());
                nustat[gr]++;
            }
            //cout<<"fin\n";
        } else if (groupe[gr].type == 2) { //SNP
            catsnp = (stat_num[j] - 21) / 4;
            if (debuglevel == 2) cout << "snp catsnp=" << catsnp << "   stat_num[" << j << "]=" << stat_num[j] << "\n";
            if (stat_num[j] < 25) {
                groupe[gr].sumstat[nustat[gr]].cat = stat_num[j];
                groupe[gr].sumstat[nustat[gr]].samp = atoi(qq[2].c_str());
                trouve = false;
                if (debuglevel == 2) cout << "statsnp[" << gr << "].size=" << statsnp[gr].size() << "\n";
                if (statsnp[gr].size() > 0) {
                    for (size_t jj = 0; jj < statsnp[gr].size(); jj++) {
                        trouve = ((statsnp[gr][jj].cat == catsnp)and (statsnp[gr][jj].samp == groupe[gr].sumstat[nustat[gr]].samp));
                        if (trouve) {
                            groupe[gr].sumstat[nustat[gr]].numsnp = jj;
                            break;
                        }
                    }
                }
                if (debuglevel == 2) cout << "trouve=" << trouve << "\n";
                if (not trouve) {
                    stsnp.cat = catsnp;
                    stsnp.samp = groupe[gr].sumstat[nustat[gr]].samp;
                    stsnp.defined = false;
                    groupe[gr].sumstat[nustat[gr]].numsnp = statsnp[gr].size();
                    cout << "statsnp[" << statsnp[gr].size() << "]   cat=" << stsnp.cat << "   samp=" << stsnp.samp << "\n";
                    statsnp[gr].push_back(stsnp);
                }
                nustat[gr]++;
            } else if ((stat_num[j] > 24)and (stat_num[j] < 33)) {
                groupe[gr].sumstat[nustat[gr]].cat = stat_num[j];
                splitwords(qq[2], "&", ss1);
                groupe[gr].sumstat[nustat[gr]].samp = atoi(ss1[0].c_str());
                groupe[gr].sumstat[nustat[gr]].samp1 = atoi(ss1[1].c_str());
                trouve = false;
                if (debuglevel == 2) cout << "statsnp[" << gr << "].size=" << statsnp[gr].size() << "\n";
                if (statsnp[gr].size() > 0) {
                    for (size_t jj = 0; jj < statsnp[gr].size(); jj++) {
                        trouve = ((statsnp[gr][jj].cat == catsnp)and (statsnp[gr][jj].samp == groupe[gr].sumstat[nustat[gr]].samp)and (statsnp[gr][jj].samp1 == groupe[gr].sumstat[nustat[gr]].samp1));
                        if (trouve) {
                            groupe[gr].sumstat[nustat[gr]].numsnp = jj;
                            break;
                        }
                    }
                }
                //cout<<"statsnp.size = "<<statsnp.size()<<"   trouve = "<<trouve<<"\n";
                if (debuglevel == 2) cout << "trouve=" << trouve << "\n";
                if (not trouve) {
                    stsnp.cat = catsnp;
                    stsnp.samp = groupe[gr].sumstat[nustat[gr]].samp;
                    stsnp.samp1 = groupe[gr].sumstat[nustat[gr]].samp1;
                    stsnp.defined = false;
                    groupe[gr].sumstat[nustat[gr]].numsnp = statsnp[gr].size();
                    cout << "statsnp[" << statsnp[gr].size() << "]   cat=" << stsnp.cat << "   samp=" << stsnp.samp << "   samp1=" << stsnp.samp1 << "\n";
                    statsnp[gr].push_back(stsnp);
                }
                //cout<<"numsnp = "<<statsnp.size()<<"\n";
                nustat[gr]++;
            } else if (stat_num[j] > 32) {
                groupe[gr].sumstat[nustat[gr]].cat = stat_num[j];
                splitwords(qq[2], "&", ss1);
                groupe[gr].sumstat[nustat[gr]].samp = atoi(ss1[0].c_str());
                groupe[gr].sumstat[nustat[gr]].samp1 = atoi(ss1[1].c_str());
                groupe[gr].sumstat[nustat[gr]].samp2 = atoi(ss1[2].c_str());
                trouve = false;
                if (debuglevel == 2) cout << "statsnp.size=" << statsnp.size() << "\n";
                if (statsnp[gr].size() > 0) {
                    for (size_t jj = 0; jj < statsnp[gr].size(); jj++) {
                        trouve = ((statsnp[gr][jj].cat == catsnp)and (statsnp[gr][jj].samp == groupe[gr].sumstat[nustat[gr]].samp)and (statsnp[gr][jj].samp1 == groupe[gr].sumstat[nustat[gr]].samp1)and (statsnp[gr][jj].samp2 == groupe[gr].sumstat[nustat[gr]].samp2));
                        if (trouve) {
                            groupe[gr].sumstat[nustat[gr]].numsnp = jj;
                            break;
                        }
                    }
                }
                if (debuglevel == 2) cout << "trouve=" << trouve << "\n";
                if (not trouve) {
                    stsnp.cat = catsnp;
                    stsnp.samp = groupe[gr].sumstat[nustat[gr]].samp;
                    stsnp.samp1 = groupe[gr].sumstat[nustat[gr]].samp1;
                    stsnp.samp2 = groupe[gr].sumstat[nustat[gr]].samp2;
                    stsnp.defined = false;
                    groupe[gr].sumstat[nustat[gr]].numsnp = statsnp[gr].size();
                    cout << "statsnp[" << gr << "][" << statsnp[gr].size() << "]   cat=" << stsnp.cat << "   samp=" << stsnp.samp << "   samp1=" << stsnp.samp1 << "   samp2=" << stsnp.samp2 << "\n";
                    statsnp[gr].push_back(stsnp);
                }
                nustat[gr]++;
            }
        }
    }
    for (gr = 1; gr <= header.ngroupes; gr++) {
        if (groupe[gr].type == 2) {
            groupe[gr].nstatsnp = statsnp[gr].size();
            cout << "groupe[gr].nstatsnp=" << groupe[gr].nstatsnp << "\n";
            if (groupe[gr].nstatsnp > 0) {
                groupe[gr].sumstatsnp = vector<StatsnpC>(groupe[gr].nstatsnp);
                for (int i = 0; i < groupe[gr].nstatsnp; i++) {
                    groupe[gr].sumstatsnp[i].cat = statsnp[gr][i].cat;
                    groupe[gr].sumstatsnp[i].samp = statsnp[gr][i].samp;
                    groupe[gr].sumstatsnp[i].samp1 = statsnp[gr][i].samp1;
                    groupe[gr].sumstatsnp[i].samp2 = statsnp[gr][i].samp2;
                    groupe[gr].sumstatsnp[i].defined = statsnp[gr][i].defined;
                }
                statsnp[gr].clear();
            }
            statsnp.clear();
        }
    }
    cout << "fin de resetstats\n";
    return true;
}

string pseudoprior2(long double x) {
    string spr;
    long double mini = 0.99999 * x, maxi = 1.00001 * x;
    spr = "UN[" + LongDoubleToString(mini) + "," + LongDoubleToString(maxi) + ",0.0,0.0";
    return spr;
}

int detphistarOK(int nsel, long double** phistar) {
    bool OK;
    int npv, ip1, ip2, nphistarOK = 0, scen = rt.scenteste - 1, k1;
    npv = rt.nparam[scen];
    for (int i = 0; i < nsel; i++) {
        OK = true;
        if (numtransf < 3) {
            for (int k = 0; k < rt.nhistparam[scen]; k++) {
                OK = ((scenario[scen].histparam[k].prior.mini < phistar[i][k])and (scenario[scen].histparam[k].prior.maxi > phistar[i][k]));
                if (not OK) break;
            }
            if (OK) {
                k1 = rt.nhistparam[scen];
                for (int k = 0; k < rt.nparamut; k++) {
                    OK = ((rt.mutparam[k].prior.mini < phistar[i][k1 + k])and (rt.mutparam[k].prior.maxi > phistar[i][k1 + k]));
                    if (not OK) break;
                }
            }
        }
        if ((scenario[scen].nconditions > 0)and (OK)) {
            for (int j = 0; j < scenario[scen].nconditions; j++) {
                ip1 = 0;
                for (int k = 0; k < rt.nparam[scen]; ++k) {
                    if (scenario[scen].condition[j].param1 == scenario[scen].histparam[k].name) {
                        k1 = k;
                        break;
                    } else {
                        if (not scenario[scen].histparam[k].prior.constant) ip1++;
                    }
                }

                ip2 = 0;
                for (int k = 0; k < rt.nparam[scen]; ++k) {
                    if (scenario[scen].condition[j].param2 == scenario[scen].histparam[k].name) {
                        break;
                    } else {
                        if (not scenario[scen].histparam[k].prior.constant) ip2++;
                    }
                }
                if (scenario[scen].histparam[k1].category < 2) {
                    if (scenario[scen].condition[j].operateur == ">") OK = (floor(0.5 + phistar[i][ip1]) > floor(0.5 + phistar[i][ip2]));
                    else if (scenario[scen].condition[j].operateur == "<") OK = (floor(0.5 + phistar[i][ip1]) < floor(0.5 + phistar[i][ip2]));
                    else if (scenario[scen].condition[j].operateur == ">=") OK = (floor(0.5 + phistar[i][ip1]) >= floor(0.5 + phistar[i][ip2]));
                    else if (scenario[scen].condition[j].operateur == "<=") OK = (floor(0.5 + phistar[i][ip1]) <= floor(0.5 + phistar[i][ip2]));
                } else {
                    if (scenario[scen].condition[j].operateur == ">") OK = (phistar[i][ip1] > phistar[i][ip2]);
                    else if (scenario[scen].condition[j].operateur == "<") OK = (phistar[i][ip1] < phistar[i][ip2]);
                    else if (scenario[scen].condition[j].operateur == ">=") OK = (phistar[i][ip1] >= phistar[i][ip2]);
                    else if (scenario[scen].condition[j].operateur == "<=") OK = (phistar[i][ip1] <= phistar[i][ip2]);
                }
                if (not OK) break;
            }
        }
        if (OK) { //cout<<nphistarOK<<"    ";
            for (int j = 0; j < npv; j++) {
                phistarOK[nphistarOK][j] = phistar[i][j];
                //cout <<phistarOK[nphistarOK][j]<<"  ";
            }
            //cout<<"\n";
            nphistarOK++;
            //if ((nphistarOK % 10)==0) cin >>c;
        }
    }
    return nphistarOK;
}

void call_loc(int npart, int nrec, int nsel, long double** ss, float* stat_obs) {
    int *avant, *apres, *egal;
    long double* qobs;
    float diff;
    string* star;
    if (debuglevel == 2) cout << "début de call_loc  header.nstat = " << header.nstat << "\n";
    qobs = new long double[header.nstat];
    star = new string[header.nstat];
    avant = new int[header.nstat];
    for (int i = 0; i < header.nstat; i++) avant[i] = 0;
    apres = new int[header.nstat];
    for (int i = 0; i < header.nstat; i++) apres[i] = 0;
    egal = new int[header.nstat];
    for (int i = 0; i < header.nstat; i++) egal[i] = 0;
    for (int p = 0; p < npart; p++) {
        for (int j = 0; j < header.nstat; j++) {
            diff = stat_obs[j] - (float)ss[p][j];
            if (diff > 0.001) avant[j]++;
            else if (diff < -0.001) apres[j]++; else egal[j]++;
        }
    }
    for (int j = 0; j < header.nstat; j++) {
        qobs[j] = (long double)(avant[j] + apres[j] + egal[j]);
        if (qobs[j] > 0.0) qobs[j] = (0.5 * (long double)egal[j] + (long double)avant[j]) / qobs[j]; else qobs[j] = -1;
        star[j] = "      ";
        if ((qobs[j] > 0.95)or (qobs[j] < 0.05)) star[j] = " (*)  ";
        if ((qobs[j] > 0.99)or (qobs[j] < 0.01)) star[j] = " (**) ";
        if ((qobs[j] > 0.999)or (qobs[j] < 0.001)) star[j] = " (***)";
        cout << setiosflags(ios::left) << setw(15) << header.statname[j] << "    (" << setiosflags(ios::fixed) << setw(8) << setprecision(4) << stat_obs[j] << ")   ";
        cout << setiosflags(ios::fixed) << setw(8) << setprecision(4) << qobs[j] << star[j] << "  ";
        cout << "\n";
    }
    string nomfiloc;
    nomfiloc = path + ident + "_mclocate.txt";
    //strcpy(nomfiloc,path);
    //strcat(nomfiloc,ident);
    //strcat(nomfiloc,"_mclocate.txt");
    cout << nomfiloc << "\n";
    time_t rawtime;
    struct tm* timeinfo;
    time(&rawtime);
    timeinfo = localtime(&rawtime);
    ofstream f12(nomfiloc.c_str(), ios::out);
    f12 << "DIYABC :                   POSTERIOR CHECKING                         " << asctime(timeinfo) << "\n";
    f12 << "Data file       : " << header.datafilename << "\n";
    f12 << "Reference table : " << rt.filename << "\n";
    f12 << "Chosen scenario : " << rt.scenteste << "\n";
    f12 << "Number of simulated data sets used to compute posterior    : " << nrec << "\n";
    f12 << "Number of simulated data sets used in the local regression : " << nsel << "\n";
    f12 << "Number of data sets simulated from the posterior           : " << npart << "\n";
    switch (numtransf) {
    case 1: f12 << "No transformation of parameters\n\n";
        break;
    case 2: f12 << "Transformation of parameters : Log\n\n";
        break;
    case 3: f12 << "Transformation of parameters : Logit\n\n";
        break;
    case 4: f12 << "Transformation of parameters : Log(tg)\n\n";
        break;
    }

    f12 << " summary           observed       proportion   \n";
    f12 << "statistics           value    (simulated<observed)\n";
    for (int j = 0; j < header.nstat; j++) {
        f12 << setiosflags(ios::left) << setw(15) << header.statname[j] << "     " << setiosflags(ios::fixed) << setw(8) << setprecision(4) << stat_obs[j] << "       ";
        f12 << setiosflags(ios::fixed) << setw(6) << setprecision(4) << qobs[j] << star[j] << "   ";
        f12 << "\n";
    }
    f12.close();
}

void call_acp(int nr, int ns, int nstat, int* numscen, long double** ssref, long double** ssphistar, float* stat_obs) {
    resACPC rACP;
    long double *pca_statobs, **pca_ss;
    pca_statobs = new long double[nstat];
    pca_ss = new long double*[ns];
    for (int i = 0; i < ns; i++) pca_ss[i] = new long double[nstat];
    rACP = ACP(nr, nstat, ssref, 1.0, 0);
    //for (int k=0;k<nstat;k++) cout<<stat_obs[k]<<"   "<<rACP.moy[k]<<"   "<<rACP.sd[k]<<"\n";
    for (int j = 0; j < rACP.nlambda; j++) {
        pca_statobs[j] = 0.0;
        for (int k = 0; k < nstat; k++) if (rACP.sd[k] > 0) pca_statobs[j] += (stat_obs[k] - rACP.moy[k]) / rACP.sd[k] * rACP.vectprop[k][j];
        cout << "  " << pca_statobs[j];
    }
    cout << "\n";
    for (int i = 0; i < ns; i++) {
        for (int j = 0; j < rACP.nlambda; j++) {
            pca_ss[i][j] = 0.0;
            for (int k = 0; k < nstat; k++) if (rACP.sd[k] > 0) pca_ss[i][j] += (ssphistar[i][k] - rACP.moy[k]) / rACP.sd[k] * rACP.vectprop[k][j];
        }
    }
    string nomfiACP;
    nomfiACP = path + ident + "_mcACP.txt";
    //strcpy(nomfiACP,path);
    //strcat(nomfiACP,ident);
    //strcat(nomfiACP,"_mcACP.txt");
    cout << nomfiACP << "\n";
    ofstream f1;
    f1.open(nomfiACP.c_str());
    f1 << setiosflags(ios::fixed) << nr << " " << rACP.nlambda;
    f1 << setprecision(3);
    for (int i = 0; i < rACP.nlambda; i++) f1 << " " << (rACP.lambda[i] / rACP.slambda);
    f1 << "\n";
    f1 << "0";
    for (int i = 0; i < rACP.nlambda; i++) f1 << " " << pca_statobs[i];
    f1 << "\n";
    for (int i = 0; i < nr; i++) {
        f1 << setiosflags(ios::fixed) << setprecision(0) << numscen[i];
        f1 << setprecision(3);
        for (int j = 0; j < rACP.nlambda; j++) f1 << " " << rACP.princomp[i][j];
        f1 << "\n";
    }
    for (int i = 0; i < ns; i++) {
        f1 << setiosflags(ios::fixed) << setprecision(0) << rt.scenteste;
        f1 << setprecision(3);
        for (int j = 0; j < rACP.nlambda; j++) f1 << " " << pca_ss[i][j];
        f1 << "\n";
    }
    f1.close();
}

void domodchec(string opt, int seed) {
    int nstatOK, nphistarOK, iprog, nprog;
    int nrec = 0, nsel = 0, ns, nrecpos = 0, newsspart = 0, npv, npvmax, nss, nsr, newrefpart, *numscen, nparamax, bidon;
    string s, s0, s1, snewstat;
    vector<string> ss, ss1;
    //float  *stat_obs;
    bool usestats, firsttime, dopca = false, doloc = false, newstat = false;
    long double **matC, **phistar, **phistarcompo, **phistarscaled;

    progressfilename = path + ident + "_progress.txt";
    scurfile = path + "first_records_of_the_reference_table_" + ident + ".txt";
    cout << "debut domodchec  options : " << opt << "\n";
    original = true;
    composite = false;
    scaled = false;
    splitwords(opt, ";", ss);
    ns = ss.size();
    numtransf = 3;
    for (int i = 0; i < ns; i++) {
        cout << ss[i] << "\n";
        s0 = ss[i].substr(0, 2);
        s1 = ss[i].substr(2);
        if (s0 == "s:") {
            splitwords(s1, ",", ss1);
            rt.nscenchoisi = ss1.size();
            if (rt.nscenchoisi > 1) rt.nscenchoisi = 1;
            rt.scenchoisi = new int[rt.nscenchoisi];
            for (int j = 0; j < rt.nscenchoisi; j++) rt.scenchoisi[j] = atoi(ss1[j].c_str());
            nrecpos = rt.nrecscen[rt.scenchoisi[0] - 1];
            cout << "scenario choisi : " << rt.scenchoisi[0] << "\n";
            rt.scenteste = rt.scenchoisi[0];
        } else if (s0 == "n:") {
            nrec = atoi(s1.c_str());
            if (nrec > nrecpos) nrec = nrecpos;
            cout << "nombre total de jeux de données considérés (pour le scénario choisi )= " << nrec << "\n";
        } else if (s0 == "m:") {
            nsel = atoi(s1.c_str());
            cout << "nombre de jeux de données considérés pour la régression locale = " << nsel << "\n";
        } else if (s0 == "t:") {
            numtransf = atoi(s1.c_str());
            switch (numtransf) {
            case 1: cout << " pas de transformation des paramètres\n";
                break;
            case 2: cout << " transformation log des paramètres\n";
                break;
            case 3: cout << " transformation logit des paramètres\n";
                break;
            case 4: cout << " transformation log(tg) des paramètres\n";
                break;
            }
        } else if (s0 == "v:") {
            cout << "" << "\n";
            snewstat = s1;
            newstat = (s1.length() > 0);
        } else if (s0 == "q:") {
            newsspart = atoi(s1.c_str());
            cout << "nombre de particules à simuler à partir du posterior = " << newsspart << "\n";
        } else if (s0 == "a:") {
            dopca = (s1.find("p") != string::npos);
            doloc = (s1.find("l") != string::npos);
            if (dopca) cout << "Perform ACP  ";
            if ((s1 == "pl")or (s1 == "lp")) cout << "et ";
            if (doloc) cout << "locate  ";
            cout << "\n";
        }
    }
    nprog = newsspart + 100;
    if ((newstat)and (dopca)) nprog += scenario.size() * 10000;
    iprog = 10;
    fprog.open(progressfilename.c_str());
    fprog << iprog << "   " << nprog << "\n";
    fprog.close();
    nstatOK = rt.cal_varstat();
    cout << "apres cal_varstat   nstatOK=" << nstatOK << "\n";
    //header.calstatobs(statobsfilename);
    //stat_obs = header.stat_obs;  cout<<"apres read_statobs\n";
    cout << "nrec=" << nrec << "     nsel=" << nsel << "\n";
    rt.alloue_enrsel(nsel);
    rt.cal_dist(nrec, nsel, &header.stat_obs[0], false, false);
    cout << "apres cal_dist\n";
    iprog += 40;
    fprog.open(progressfilename.c_str());
    fprog << iprog << "   " << nprog << "\n";
    fprog.close();
    det_numpar();
    cout << "apres det_numpar\n";
    rempli_mat(nsel, &header.stat_obs[0]);
    cout << "apres rempli_mat\n";
    if (not deltanul) matC = cal_matC(nsel);
    recalparamO(nsel);
    cout << "apres recalparam\n";
    if (not deltanul) {
        rempli_parsim(nsel, nparamcom);
        cout << "apres rempli_parsim(O)\n";
        local_regression(nsel, nparamcom, matC);
        cout << "apres local_regression\n";
    }
    iprog += 20;
    fprog.open(progressfilename.c_str());
    fprog << iprog << "   " << nprog << "\n";
    fprog.close();
    phistar = new long double* [nsel];
    for (int i = 0; i < nsel; i++) phistar[i] = new long double[nparamcom];
    if (not deltanul) calphistarO(nsel, phistar);else copyphistar(nsel, nparamcom, phistar);
    cout << "apres calphistar\n";
    det_nomparam();
    phistarcompo = new long double*[nsel];
    phistarscaled = new long double*[nsel];
    for (int i = 0; i < nsel; i++) {
        phistarcompo[i] = new long double[0];
        phistarscaled[i] = new long double[0];
    }
    savephistar(nsel, path, ident, phistar, phistarcompo, phistarscaled);
    cout << "apres savephistar\n";
    iprog += 20;
    fprog.open(progressfilename.c_str());
    fprog << iprog << "   " << nprog << "\n";
    fprog.close();
    phistarOK = new long double*[nsel];
    for (int i = 0; i < nsel; i++) phistarOK[i] = new long double[rt.nparam[rt.scenteste - 1]];
    cout << "scenario[rt.scenteste-1].nparam = " << scenario[rt.scenteste - 1].nparam << "\n";
    nphistarOK = detphistarOK(nsel, phistar);
    cout << "apres detphistarOK  nphistarOK=" << nphistarOK << "\n";
    //for (int i=0;i<nsel;i++) delete [] phistar[i];delete phistar;

    cout//<<"naparamcom="<<nparamcom<<"   nparcompo="<<nparcompo<<"   nenr="<<nenr
        << "   nphistarOK=" << nphistarOK << "   nstat=" << header.nstat << "\n";
    //cout <<"DEBUG: j'arrête là." << endl; exit(1);
    if (nphistarOK < 100) {
        cout << "Not enough suitable particles (" << nphistarOK << ")to perform model checking. Stopping computations." << endl;
        exit(1);
    }
    npv = rt.nparam[rt.scenteste - 1];
    npvmax = npv;
    for (int i = 0; i < (int)scenario.size(); i++) if (npvmax < rt.nparam[i]) npvmax = rt.nparam[i];
    cout << "npv=" << npv << "    nenr=" << nenr << "   nstat=" << header.nstat << "\n";
    cout << "npvmax = " << npvmax << "\n";
    //delete []enreg;
    nss = 0;
    firsttime = true;
    cout << "ns=" << ns << "\n";
    //cout<<phistarOK[0][0]<<"\n";
    if (newstat) usestats = resetstats(snewstat);
    if (not usestats) cout << "problème dans la construction des nouvelles stats\n";
    cout << "header.nstat = " << header.nstat << "\n";
    enreg = new enregC[nenr];
    cout << "apres new enregC[nenr]\n";
    for (int p = 0; p < nenr; p++) {
        enreg[p].stat = vector<float>(header.nstat);
        enreg[p].param = vector<float>(npvmax);
        enreg[p].numscen = rt.scenteste;
    }
    cout << "apres dimensionnement des enreg nstat=" << header.nstat << "   nparam=" << npv << "   nenr=" << nenr << "\n";
    ssphistar = new long double*[newsspart];
    for (int i = 0; i < newsspart; i++) ssphistar[i] = new long double[header.nstat];
    cout << "newstat =" << newstat << "   newsspart=" << newsspart << "     nenr=" << nenr << "\n";
    if (nenr > newsspart) nenr = newsspart;
    while (nss < newsspart) {
        if (debuglevel == 2) cout << "avant dosimulphistar\n";
        ps.dosimulphistar(nenr, false, multithread, firsttime, rt.scenteste, seed, nphistarOK);
        if (debuglevel == 2) cout << "apres dosimulphistar\n";
        for (int i = 0; i < nenr; i++) {
            for (int j = 0; j < header.nstat; j++) ssphistar[i + nss][j] = enreg[i].stat[j];
            for (int j = 0; j < header.nstat; j++) cout << ssphistar[i + nss][j] << "   ";
            cout << "\n";
        }
        firsttime = false;
        nss += nenr;
        iprog += nenr;
        fprog.open(progressfilename.c_str());
        fprog << iprog << "   " << nprog << "\n";
        fprog.close();
        //cout<<nss<<"\n";
    }
    if (newstat) {
        //cout<<"\n\n\nAVANT CALSTATOBS\n";
        header.calstatobs(statobsfilename);/*stat_obs = header.stat_obs;*/
        //cout<<"\nAPRES CALSTATOBS\n\n\n";
    }
    if (doloc) call_loc(newsspart, nrec, nsel, ssphistar, &header.stat_obs[0]);
    if (dopca) {
        if (newstat) {
            header.readHeader(headerfilename);
            cout << "apres readHeader scenario.size()= " << scenario.size() << "\n";
            usestats = resetstats(snewstat);
            newrefpart = 10000;
            cout << "newrefparts=" << newrefpart << "\n";
            numscen = new int[newrefpart];
            ssref = new long double*[newrefpart];
            for (int i = 0; i < newrefpart; i++) ssref[i] = new long double[header.nstat];
            nsr = 0;
            firsttime = true;
            cout << "avant le while (nsr<newrefpart)\n";
            while (nsr < newrefpart) {
                ps.dosimultabref(nenr, false, multithread, firsttime, 0, seed, 3);
                for (int i = 0; i < nenr; i++) {
                    numscen[i + nsr] = enreg[i].numscen;
                    for (int j = 0; j < header.nstat; j++) ssref[i + nsr][j] = enreg[i].stat[j];
                    if (debuglevel == 2) {
                        for (int j = 0; j < header.nstat; j++) cout << ssref[i + nsr][j] << "   ";
                        cout << "\n";
                    }
                }
                firsttime = false;
                nsr += nenr;
                iprog += nenr;
                fprog.open(progressfilename.c_str());
                fprog << iprog << "   " << nprog << "\n";
                fprog.close();
                cout << nsr << "\n";
            }
        } else {
            rt.openfile2();
            enregC enr;
            nparamax = 0;
            for (int i = 0; i < rt.nscen; i++) if (rt.nparam[i] > nparamax) nparamax = rt.nparam[i];
            enr.stat = vector<float>(rt.nstat);
            enr.param = vector<float>(nparamax);
            newrefpart = scenario.size() * 10000;
            cout << "newrefparts=" << newrefpart << "\n";
            numscen = new int[newrefpart];
            ssref = new long double*[newrefpart];
            for (int i = 0; i < newrefpart; i++) ssref[i] = new long double[header.nstat];
            nsr = 0;
            while (nsr < newrefpart) {
                bidon = rt.readrecord(&enr);
                if (bidon != 0) cout << "probleme à la lecture du reftable\n";
                numscen[nsr] = enr.numscen;
                for (int j = 0; j < rt.nstat; j++) ssref[nsr][j] = enr.stat[j];
                nsr++;
            }
        }
        call_acp(newrefpart, newsspart, header.nstat, numscen, ssref, ssphistar, &header.stat_obs[0]);
    }
    iprog += 10;
    fprog.open(progressfilename.c_str());
    fprog << iprog << "   " << nprog << "\n";
    fprog.close();
}

