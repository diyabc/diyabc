/*
 * simfile.cpp
 *
 *  Created on: 20 september 2011
 *      Author: cornuet
 */
#include <vector>
#include <string>

#include "header.hpp"
#include "reftable.hpp"
#include "matrices.hpp"
#include "mesutils.hpp"
#include "particleset.hpp"
#include "modchec.hpp"
/*
#ifndef HEADER
#include "header.cpp"
#define HEADER
#endif

#ifndef VECTOR
#include <vector>
#define VECTOR
#endif

#ifndef MATRICES
#include "matrices.cpp"
#define MATRICES
#endif

#ifndef MESUTILS
#include "mesutils.cpp"
#define MESUTILS
#endif


#ifndef PARTICLESET
#include "particleset.cpp"
#define PARTICLESET
#endif
*/
#include "simfile.hpp"
#include <stdlib.h>
#include <math.h>
#include <iomanip>

extern ParticleSetC ps;
extern enregC* enreg;
extern bool multithread;
extern string progressfilename, path, paramfilename,statfilename;
extern HeaderC header;
extern ReftableC rt;
extern int nenr;
long double **phistar;
extern long double **phistarOK;
extern DataC dataobs,datasim;
extern vector <ScenarioC> scenario;

void dosimfile(int seed){
	string nomfigp;
	int  ntest;
	string opt,s,s0,s1,sg;
	string *sgp;
	string nomfisim;
	FILE *flog, *fgp;
	cout<<"path="<<path<<"\n";
	progressfilename =path + "progress.txt";
	nomfisim=header.datafilename;
	ntest = header.nsimfile;cout<<"SIMFILE filetype="<<datasim.filetype<<"\n";
	if (datasim.filetype==0) cout<<"nombre de fichiers genepop à simuler = "<<ntest<<"\n";
	if (datasim.filetype>0) cout<<"nombre de fichiers snp à simuler = "<<ntest<<"\n";
	if (datasim.filetype==0) sgp = ps.simulgenepop(ntest, multithread, seed);
	if (datasim.filetype==1) sgp = ps.simuldataSNP(ntest, multithread, seed);
	if (datasim.filetype==2) {sgp = ps.simuldataSNPool(ntest, multithread, seed);cout<<"simfile POOLSEQ\n";}
	//cout<<"apres les simulations\n";
	for (int i=0;i<ntest;i++) {
		string sn=IntToString(i+1);
		if (i<9)  sn ="0"+sn;
		if (i<99) sn="0"+sn;
		sn="_"+sn;
		if (sgp[i]=="") cout<<"une erreur s'est produite lors de la simulation du fichier numero "<<i+1<<"\n";
		else {
			nomfigp = path + nomfisim + sn;
			if (dataobs.filetype==0) nomfigp +=".mss";
			else nomfigp +=".snp";
			cout<<"écriture du fichier "<<nomfigp<<"\n";
			sgp[i] +="\n";
			fgp = fopen(nomfigp.c_str(),"w");fprintf(fgp,"%s", sgp[i].c_str());fclose(fgp);
		}
	}
	flog=fopen(progressfilename.c_str(),"w");fprintf(flog,"OK");fclose(flog);
}

    int detpstarOK(int nsel, int scen, int npv, long double **phistar) {
		//char c;
        bool OK;
        int ip1,ip2,nphistarOK=0,k1,k0;
		//cout<<"depstarOK 0  npv="<<npv<<"   nsel="<<nsel<<"   scen="<<scen<<"  nparamhist="<<scenario[scen].nparam<<"   nparamut="<<header.nparamut<<"\n";
		//cout<<"nparamvar="<<scenario[scen].nparamvar<<"\n";
		for (int i=0;i<nsel;i++) {
			//for(int j=0;j<scenario[scen].nparam+header.nparamut;j++) cout<<phistar[i][j]<<"  ";cout<<"\n";
            OK=true;k0=0;
			for (int k=0;k<scenario[scen].nparam;k++) {
				while (scenario[scen].histparam[k].prior.constant) k++;
				//cout<<"k="<<k<<"  k0="<<k0<<"  "<<scenario[scen].histparam[k].name<<"\n";
				OK=((scenario[scen].histparam[k].prior.mini<=phistar[i][k0])and(scenario[scen].histparam[k].prior.maxi>=phistar[i][k0]));
				if (not OK){
					cout<<scenario[scen].histparam[k].name<<"  ["<<scenario[scen].histparam[k].prior.mini<<","<<scenario[scen].histparam[k].prior.maxi<<"]  "<<phistar[i][k0]<<"\n";
					break;
				} else k0++;
				//cout<<"fin\n";
			}
			//cout<<"depstarOK 1\n";
			//if (OK) cout<<"  histparam OK\n";else {cout<<" histparam["<<k0-1<<"] not OK\n";}
			if (OK){
				k1=scenario[scen].nparamvar;//cout<<"k1="<<k1<<"\n";
				k1=k0;
				for (int k=0;k<header.nparamut;k++) {
					OK=((header.mutparam[k].prior.mini<=phistar[i][k1+k])and(header.mutparam[k].prior.maxi>=phistar[i][k1+k]));
					if (not OK ){
						cout<<header.mutparam[k].name<<"   phistar["<<i<<"]["<<k1+k<<"]="<<phistar[i][k1+k]<<"   mini="<<header.mutparam[k].prior.mini<<"   maxi="<<header.mutparam[k].prior.maxi;
						//for(int j=0;j<scenario[scen].nparam+header.nparamut;j++) cout<<phistar[i][j]<<"  ";cout<<"\n";
						break;
					}
				}
			}
			//cout<<"depstarOK 1";if (OK) cout<<"  mutparam OK\n";else cout<<"mutparam not OK\n";
			//cout<<"nconditions="<<scenario[scen].nconditions<<"\n";
            if ((scenario[scen].nconditions>0)and(OK)) {
                for (int j=0;j<scenario[scen].nconditions;j++) {
                    ip1=0;
                    for(int k =0; k < scenario[scen].nparam; ++k){
                    	if(scenario[scen].condition[j].param1 == scenario[scen].histparam[k].name){
							k1=k;
                    		break;
                    	} else {
                    		if(not scenario[scen].histparam[k].prior.constant) ip1++;
                    	}
                    }

                    ip2=0;
                    for(int k =0; k < scenario[scen].nparam; ++k){
                    	if(scenario[scen].condition[j].param2 == scenario[scen].histparam[k].name){
                    		break;
                    	} else {
                    		if(not scenario[scen].histparam[k].prior.constant) ip2++;
                    	}
                    }
                    if (scenario[scen].histparam[k1].category<2){
						if (scenario[scen].condition[j].operateur==">")       OK=(floor(0.5+phistar[i][ip1]) >  floor(0.5+phistar[i][ip2]));
						else if (scenario[scen].condition[j].operateur=="<")  OK=(floor(0.5+phistar[i][ip1]) <  floor(0.5+phistar[i][ip2]));
						else if (scenario[scen].condition[j].operateur==">=") OK=(floor(0.5+phistar[i][ip1]) >= floor(0.5+phistar[i][ip2]));
						else if (scenario[scen].condition[j].operateur=="<=") OK=(floor(0.5+phistar[i][ip1]) <= floor(0.5+phistar[i][ip2]));
					} else {
						if (scenario[scen].condition[j].operateur==">")       OK=(phistar[i][ip1] >  phistar[i][ip2]);
						else if (scenario[scen].condition[j].operateur=="<")  OK=(phistar[i][ip1] <  phistar[i][ip2]);
						else if (scenario[scen].condition[j].operateur==">=") OK=(phistar[i][ip1] >= phistar[i][ip2]);
						else if (scenario[scen].condition[j].operateur=="<=") OK=(phistar[i][ip1] <= phistar[i][ip2]);
					}
					if (not OK) cout<<"condition "<<j<<"non remplie\n";
                    if (not OK) break;
                }
            }
			//cout<<"depstarOK 2";if (OK) cout<<"  OK\n";else cout<<"  not OK\n";
            if (OK) { //cout<<nphistarOK<<"    ";
                for (int j=0;j<npv;j++) {
                    phistarOK[nphistarOK][j] = phistar[i][j];
                    //cout <<phistarOK[nphistarOK][j]<<"  ";
                }
                //cout<<"\n";
                nphistarOK++;
				//if ((nphistarOK % 10)==0) cin >>c;
            } else cout <<" -->rejecting particle "<<i<<"\n";
        }
        return nphistarOK;
    }
    
void dosimstat(int seed) {
	int nsel,nsel0=0,nphistarOK,*numscen,nss,*nrecscen,*nparam,npart,nn,j0;
	double **stat;
	long double **phistar0;
	string ligne,*ss;
	cout<<"path="<<path<<"\n";
	cout<<"paramfilename ="<<paramfilename<<"\n";
	cout<<header.nscenarios<<" scenarios  "<<header.nparamtot<<" historical and "<<header.nparamut<<" mutation parameters    and "<<header.nstat<<" sumstat\n";
	nrecscen = new int[header.nscenarios];
	nparam = new int[header.nscenarios];
	for (int i=0;i<header.nscenarios;i++) {nrecscen[i]=0;nparam[i]=0;}
	ifstream fp0,fp;
	fp0.open(paramfilename.c_str());
	while (not fp0.eof()){
		getline(fp0,ligne);if (ligne.length()>5) nsel0++;
	}
	fp0.close();
	cout<<"File "<<paramfilename<<" contains "<<nsel0<<" lignes\n";
	numscen = new int[nsel0];
	phistar0 = new long double*[nsel0];
	for (int i=0;i<nsel0;i++) phistar0[i] = new long double[header.nparamtot+header.nparamut];
	phistar = new long double*[nsel0];
	for (int i=0;i<nsel0;i++) phistar[i] = new long double[header.nparamtot+header.nparamut];
	phistarOK = new long double*[nsel0];
	for (int i=0;i<nsel0;i++) phistarOK[i] = new long double[header.nparamtot+header.nparamut];
	fp.open(paramfilename.c_str());
	nsel0=0;
	while (not fp.eof()){
		getline(fp,ligne);if (ligne.length()>5) { 
			ss = splitwords(ligne," ",&nss);//cout<<"nss="<<nss<<"\n";
			numscen[nsel0] = atoi(ss[0].c_str());//cout<<"numscen["<<nsel0<<"]="<<numscen[nsel0]<<"\n";
			nrecscen[numscen[nsel0]-1]++;
			if (nparam[numscen[nsel0]-1]==0){for (int j=1;j<nss;j++){if (ss[j] !="NA")nparam[numscen[nsel0]-1]++;}}
			//cout<<"nparam["<<numscen[nsel0]-1<<"]="<<nparam[numscen[nsel0]-1]<<"\n";
			/*if (nsel0==0){
				for (int j=1;j<nss;j++) cout<<ss[j]<<"  "<<(long double)atof(ss[j].c_str())<<"\n";
				cout<<"nparam="<<nparam[numscen[nsel0]-1]<<"\n\n";}*/
			j0=0;
			for (int j=1;j<nss;j++) {
				if (ss[j] !="NA") {phistar0[nsel0][j0]=(long double)atof(ss[j].c_str());j0++;}
			}
			/*if (nsel0==0){
				for (int j=0;j<nparam[numscen[nsel0]-1];j++) cout<< phistar0[nsel0][j]<<"\n";
				cout<<"nparam="<<nparam[numscen[nsel0]-1]<<"\n\n";}*/
			nsel0++;
		}
	}
	fp.close();
	//cout<<"fin de la lecture du fichier "<<paramfilename<<"\n";
	ofstream fs;
	fs.open(statfilename.c_str());
	for (int iscen=0;iscen<header.nscenarios;iscen++) {
		nsel = nrecscen[iscen];cout<<"\nscenario "<<iscen+1<<"    nsel="<<nsel<<"   nparam="<<nparam[iscen]<<"\n";
		nsel=0;
		for (int i=0;i<nsel0;i++) {
			if (numscen[i]-1==iscen) {
				for (int j=0;j<nparam[iscen];j++) phistar[nsel][j]=phistar0[i][j];
				//if (nsel==0) for (int j=0;j<nparam[iscen];j++) cout<<phistar0[i][j]<<"  "<<phistar[nsel][j]<<"\n";
				nsel++;
			}
		}
		rt.scenteste=iscen+1;
		//cout<<"avant detphistar\n";
		nphistarOK=detpstarOK(nsel,iscen,nparam[iscen],phistar);
		cout<<"nphistarOK="<<nphistarOK<<"\n";
		
		stat=new double*[nsel];
		for (int i=0;i<nsel;i++) {
			stat[i]=new double[header.nstat];
			for (int j=0;j<header.nstat;j++) stat[i][j]=-99.0;
		}
		npart=nenr;nn=0;
		while (nphistarOK>npart+nn) {
			ps.dosimulstat(nn,npart,false,multithread,iscen+1,seed,stat);
			nn +=npart;
		}
		//cout<<"nsel="<<nsel<<"   nphistarOK="<<nphistarOK<<"   nn="<<nn<<"\n";
		if (nn<nphistarOK) ps.dosimulstat(nn,nphistarOK-nn,false,multithread,iscen+1,seed,stat);
		for (int i=0;i<nsel;i++) {
			if (stat[i][0] !=-99.0) {
				fs<<iscen+1<<"  ";//cout<<iscen+1<<"  ";
				for (int j=0;j<header.nstat;j++) fs<<fixed<<setw(10)<<setprecision(4)<<stat[i][j];
				fs<<"\n";
				//for (int j=0;j<header.nstat;j++) cout<<fixed<<setw(14)<<setprecision(6)<<"  "<<stat[i][j];
				//cout<<"\n";
			}
		}
	}
	fs.close();
	cout<<"\nSimulated summary statistics are in file "<<statfilename<<"\n";
}

