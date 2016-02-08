//// ***** Programme de test ****** //////

#ifndef IOSTREAM
#include <iostream>
#define IOSTREAM
#endif
#ifndef STRINGH
#include <string.h>
#define STRINGH
#endif
#ifndef UNISTDH
#include <unistd.h>
#define UNISTDH
#endif
#ifndef PARTICLESET
#include "particleset.cpp"
#define PARTICLESET
#endif
#ifndef REFTABLE
#include "reftable.cpp"
#define REFTABLE
#endif
#ifndef ESTIMPARAM
#include "estimparam.cpp"
#define ESTIMPARAM
#endif
#ifndef COMPARSCEN
#include "comparscen.cpp"
#define COMPARSCEN
#endif
#ifndef BIAS
#include "bias.cpp"
#define BIAS
#endif
#ifndef CONF
#include "conf.cpp"
#define CONF
#endif
#ifndef SYS_TIMEH
#include <sys/time.h>
#define SYS_TIMEH
#endif
#ifndef ACPLOC
#include "acploc.cpp"
#define ACPLOC
#endif
#ifndef MODCHEC
#include "modchec.cpp"
#define MODCHEC
#endif
#ifndef SYSSTATH
#include <sys/stat.h>
#define SYSSTATH
#endif

ReftableC rt;
HeaderC header;    
ParticleSetC ps;
struct stat stFileInfo; 

char *headerfilename, *reftablefilename,*datafilename,*statobsfilename, *reftablelogfilename,*path,*ident,*stopfilename;
bool multithread=false;
int nrecneeded,nenr=1000,niter=1000;
int debuglevel=0;
string sremtime;
double clock_zero=0.0,debut,duree,debutf,dureef,time_file=0.0,time_reftable=0.0,debutr,dureer,remtime;

/**
* lecture du fichier header.txt, calcul des stat_obs et lecture de l'entête de reftable.bin' 
*/

int readheaders() {
    int k;
    header.readHeader(headerfilename);                                   if (debuglevel==1) cout<<"apres header.readHeader\n";
    header.calstatobs(statobsfilename);                                  if (debuglevel==1) cout <<"apres header.calstatobs\n";
    datafilename=strdup(header.datafilename.c_str());                    if (debuglevel==1) cout<<"datafile name : "<<header.datafilename<<"\n";
    k=rt.readheader(reftablefilename,reftablelogfilename,datafilename);  if (debuglevel==1) cout<<"apres rt.readheader k="<<k<<"\n";
    if (k==0) {rt.sethistparamname(header);cout<<"sethistparamname"<<"\n";}
    return k;
}

/**
* programme principal : lecture et interprétation des options et lancement des calculs choisis
*/

int main(int argc, char *argv[]){
    char *dataname,*headername,*reftablename,aaa,*winref,*linref;
    bool firsttime;
	int k,seed,numtest,***n,scen,j0,bidon;
	double **paramstat,hetmoy,hetheo,pidmoy,pidtheo,***sx,***sx2,x,x0,x1,x2,x3,dx,dx2;
	int optchar;
    char action='a';
    bool flagp=false,flagi=false,flags=false,simOK,stoprun=false;
    string message,***paramname;
    FILE *flog;
    enregC enr;
    string star;
        
    debut=walltime(&clock_zero);
	while((optchar = getopt(argc,argv,"p:n:s:m")) !=-1) {
         
	  switch (optchar) {

        case 'n' : numtest = atoi(optarg);
                   break;
        
        case 'p' :
                  headerfilename      = new char[strlen(optarg)+13];
                  reftablefilename    = new char[strlen(optarg)+15];
                  reftablelogfilename = new char[strlen(optarg)+15];
                  statobsfilename     = new char[strlen(optarg)+14];
                  stopfilename        = new char[strlen(optarg)+13];
                  winref              = new char[strlen(optarg)+20];
                  linref              = new char[strlen(optarg)+20];
                  path                = new char[strlen(optarg)+1];
                  strcpy(path,optarg);
                  strcpy(headerfilename,optarg);
                  strcpy(reftablefilename,optarg);
                  strcpy(reftablelogfilename,optarg);
                  strcpy(statobsfilename,optarg);
                  strcpy(stopfilename,optarg);
                  strcpy(winref,optarg);
                  strcpy(linref,optarg);
                  strcat(headerfilename,"header.txt");
                  strcat(reftablefilename,"reftable.bin");
                  strcat(reftablelogfilename,"reftable.log");
                  strcat(statobsfilename,"statobs.txt");
                  strcat(stopfilename,".stop");
                  flagp=true;
                  if (stat(stopfilename,&stFileInfo)==0) remove(stopfilename);
                  break;
		   
        case 's' :seed=atoi(optarg);
                  flags=true;
                  break;
            
        case 'm' :multithread=true;
                  break;

	    }
	}
	 if (not flagp) {cout << "option -p is compulsory\n";exit(1);}
     if (not flags) seed=time(NULL);	
	switch (numtest) {
    
      case 1 :   k=readheaders();
                 if (k!=0)  {cout<<"cannot find reftable file\n"; exit(1);}
                 enreg = new enregC[nenr];
                 for (int p=0;p<nenr;p++) {
                     enreg[p].stat = new float[header.nstat];
                     enreg[p].param = new float[header.nparamtot+3*header.ngroupes];
                     enreg[p].numscen = 1;
                                  }
                 cout<<"TEST DE l'HETEROZYGOTIE MOYENNE ET DE LA PROBABILITE MOYENNE D'IDENTITE\n D'UN LOCUS SOUS SRICT SMM\n";
                 
                 hetheo=1.0-1.0/sqrt(1.04);pidtheo=0.9805595780769629;
                 cout<<"\ntheta=0.02"<<"\n";
                 firsttime=true;stoprun=false;
                 header.scenario[0].histparam[0].prior.mini=9.9;header.scenario[0].histparam[0].prior.maxi=10.1;
                 hetmoy=0.0;pidmoy=0.0;
                 for (k=0;k<niter;k++) { 
                     ps.dosimultabref(header,nenr,false,multithread,firsttime,0,seed,true,true);
                     firsttime=false;
                     seed +=nenr;
                     for (int i=0;i<nenr;i++) {hetmoy += enreg[i].stat[0];pidmoy += enreg[i].stat[1];}
                 }
                 hetmoy /= (double)(nenr*niter);pidmoy /= (double)(nenr*niter);
                 cout<<"hétérozygotie attendue = "<<hetheo<<"    hétérozygotie moyenne = "<<hetmoy;
                 if(fabs(hetmoy-hetheo)<0.001) cout <<"   test OK\n"; else cout<<"   PROBLEME !!!\n";
                 cout<<"prob(identite) attendue = "<<pidtheo<<"    prob(identite) moyenne = "<<pidmoy;
                 if(fabs(pidmoy-pidtheo)<0.001) cout <<"   test OK\n"; else cout<<"   PROBLEME !!!\n";
                 
                 hetheo=1.0-1.0/sqrt(1.4);pidtheo=0.845026245742127;
                 cout<<"\ntheta=0.2"<<"\n";
                 firsttime=true;stoprun=false;
                 header.scenario[0].histparam[0].prior.mini=99.9;header.scenario[0].histparam[0].prior.maxi=100.1;
                 hetmoy=0.0;pidmoy=0.0;
                 for (k=0;k<niter;k++) { 
                     ps.dosimultabref(header,nenr,false,multithread,firsttime,0,seed,true,true);
                     firsttime=false;
                     seed +=nenr;
                     for (int i=0;i<nenr;i++) {hetmoy += enreg[i].stat[0];pidmoy += enreg[i].stat[1];}
                 }
                 hetmoy /= (double)(nenr*niter);pidmoy /= (double)(nenr*niter);
                 cout<<"hétérozygotie attendue = "<<hetheo<<"    hétérozygotie moyenne = "<<hetmoy;
                 if(fabs(hetmoy-hetheo)<0.001) cout <<"   test OK\n"; else cout<<"   PROBLEME !!!\n";
                 cout<<"prob(identite) attendue = "<<pidtheo<<"    prob(identite) moyenne = "<<pidmoy;
                 if(fabs(pidmoy-pidtheo)<0.001) cout <<"   test OK\n"; else cout<<"   PROBLEME !!!\n";
                 
                 hetheo=1.0-1.0/sqrt(5);pidtheo=0.44707336192772973;
                 cout<<"\ntheta=2.0"<<"\n";
                 firsttime=true;stoprun=false;
                 header.scenario[0].histparam[0].prior.mini=999.9;header.scenario[0].histparam[0].prior.maxi=1000.1;
                 hetmoy=0.0;pidmoy=0.0;
                 for (k=0;k<niter;k++) {  
                     ps.dosimultabref(header,nenr,false,multithread,firsttime,0,seed,true,true);
                     firsttime=false;
                     seed +=nenr;
                     for (int i=0;i<nenr;i++) {hetmoy += enreg[i].stat[0];pidmoy += enreg[i].stat[1];}
                 }
                 hetmoy /=(double) (nenr*niter);pidmoy /= (double)(nenr*niter);
                 cout<<"hétérozygotie attendue = "<<hetheo<<"    hétérozygotie moyenne = "<<hetmoy;
                 if(fabs(hetmoy-hetheo)<0.001) cout <<"   test OK\n"; else cout<<"   PROBLEME !!!\n";
                 cout<<"prob(identite) attendue = "<<pidtheo<<"    prob(identite) moyenne = "<<pidmoy;
                 if(fabs(pidmoy-pidtheo)<0.001) cout <<"   test OK\n"; else cout<<"   PROBLEME !!!\n";
                 
                  hetheo=1.0-1.0/sqrt(41);pidtheo=0.1561449771153935;
                 cout<<"\ntheta=20.0"<<"\n";
                 firsttime=true;stoprun=false;
                 header.scenario[0].histparam[0].prior.mini=9999.9;header.scenario[0].histparam[0].prior.maxi=10000.1;
                 hetmoy=0.0;pidmoy=0.0;
                 for (k=0;k<niter;k++) {  
                     ps.dosimultabref(header,nenr,false,multithread,firsttime,0,seed,true,true);
                     firsttime=false;
                     seed +=nenr;
                     for (int i=0;i<nenr;i++) {hetmoy += enreg[i].stat[0];pidmoy += enreg[i].stat[1];}
                 }
                 hetmoy /= (double)(nenr*niter);pidmoy /= (double)(nenr*niter);
                 cout<<"hétérozygotie attendue = "<<hetheo<<"    hétérozygotie moyenne = "<<hetmoy;
                 if(fabs(hetmoy-hetheo)<0.001) cout <<"   test OK\n"; else cout<<"   PROBLEME !!!\n";
                 cout<<"prob(identite) attendue = "<<pidtheo<<"    prob(identite) moyenne = "<<pidmoy;
                 if(fabs(pidmoy-pidtheo)<0.001) cout <<"   test OK\n"; else cout<<"   PROBLEME !!!\n";
              
                 break;
                 
      case 2 :   strcat(winref,"reftable.win");strcat(linref,"reftable.lin");
                 cout<<"winref = "<<winref<<"\n";
                 cout<<"linref = "<<linref<<"\n";
                 paramname = new string**[2];
                 rename(linref,reftablefilename);
                 k=readheaders();
                 if (k!=0)  {cout<<"cannot find reftable file\n"; exit(1);}
                 
                     enr.stat = new float[rt.nstat];
                     enr.param = new float[header.nparamtot+3*header.ngroupes];
                      
                 cout<<"COMPARAISON DE REFTABLE version 1 (Delphi) / version 2 (C++) \n";
                 sx = new double**[2];sx2 = new double**[2];n = new int**[2];
                 for (int i=0;i<2;i++) {
                     sx[i] = new double*[rt.nscen];
                     sx2[i] = new double*[rt.nscen];
                     n[i] = new int*[rt.nscen];
                     for (int j=0;j<rt.nscen;j++) {
                         sx[i][j] = new double[rt.nparam[j]+rt.nstat];
                         for (int k=0;k<rt.nparam[j]+rt.nstat;k++) sx[i][j][k]=0.0;
                         sx2[i][j] = new double[rt.nparam[j]+rt.nstat];
                         for (int k=0;k<rt.nparam[j]+rt.nstat;k++) sx2[i][j][k]=0.0;
                         n[i][j] = new int[rt.nparam[j]+rt.nstat];
                         for (int k=0;k<rt.nparam[j]+rt.nstat;k++) n[i][j][k]=0;
                     }
                 }
                 paramname[0] = new string*[rt.nscen];
                 for (scen=0;scen<rt.nscen;scen++) paramname[0][scen] = new string[rt.nparam[scen]];
                 for (scen=0;scen<rt.nscen;scen++) {
                     cout <<"\nscenario "<<scen<<"\n"; 
                     for (int j=0;j<rt.nhistparam[scen];j++) cout<<rt.histparam[scen][j].name<<"  ";cout<<"\n";
                     for (int j=0;j<rt.nhistparam[scen];j++) paramname[0][scen][j] = rt.histparam[scen][j].name;
                 }
                 cout<<"Lecture du fichier obtenue avec la version 2 (C++)\n";
                 rt.openfile2();
                 for (int p=0;p<rt.nrec;p++) {
                     bidon=rt.readrecord(&enr);
                     scen=enr.numscen-1;
                     for (int j=0;j<rt.nparam[scen];j++) {
                         n[0][scen][j]++;
                         x=(double)enr.param[j];
                         sx[0][scen][j] += x;
                         sx2[0][scen][j] += x*x;
                     }
                     j0=rt.nparam[scen];
                     for (int j=0;j<rt.nstat;j++) {
                         n[0][scen][j0+j]++;
                         x=(double)enr.stat[j];
                         sx[0][scen][j0+j] += x;
                         sx2[0][scen][j0+j] += x*x;
                     
                     } 
                 }
                 rt.closefile();   cout<<"apres la lecture des "<<rt.nrec<<" enregistrements\n";
                 rename(reftablefilename,linref);
                 rename(winref,reftablefilename);
                 k=readheaders();
                 if (k!=0)  {cout<<"cannot find reftable file\n"; exit(1);}
                 cout<<"Lecture du fichier obtenue avec la version 1 (Delphi)";
                 rt.openfile2();
                 for (int p=0;p<rt.nrec;p++) {
                     bidon=rt.readrecord(&enr);
                     scen=enr.numscen-1;
                     for (int j=0;j<rt.nparam[scen];j++) {
                         n[1][scen][j]++;
                         x=(double)enr.param[j];
                         sx[1][scen][j] += x;
                         sx2[1][scen][j] += x*x;
                     }
                     j0=rt.nparam[scen];
                     for (int j=0;j<rt.nstat;j++) {
                         n[1][scen][j0+j]++;
                         x=(double)enr.stat[j];
                         sx[1][scen][j0+j] += x;
                         sx2[1][scen][j0+j] += x*x;
                     
                     } 
                 }
                 rt.closefile();   cout<<"apres la lecture des "<<rt.nrec<<" enregistrements\n";
                 rename(reftablefilename,winref);
                 paramname[1] = new string*[rt.nscen];
                 for (scen=0;scen<rt.nscen;scen++) paramname[1][scen] = new string[rt.nparam[scen]];
                 for (scen=0;scen<rt.nscen;scen++) {
                     cout <<"\nscenario "<<scen<<"\n"; 
                     for (int j=0;j<rt.nhistparam[scen];j++) cout<<rt.histparam[scen][j].name<<"  ";cout<<"\n";
                     for (int j=0;j<rt.nhistparam[scen];j++) paramname[1][scen][j] = rt.histparam[scen][j].name;
                     for (int j=0;j<rt.nparamut;j++) paramname[1][scen][j+rt.nhistparam[scen]] = rt.mutparam[j].name;
                 }
                 for (scen=0;scen<rt.nscen;scen++) {
                     cout << "\nScenario "<<scen+1<<"\n";
                     for (int j=0;j<rt.nparam[scen];j++) {
                         x0 = sx[0][scen][j]/(double)n[0][scen][j];
                         x2 = sx2[0][scen][j]/(double)n[0][scen][j];
                         //k=0;while (paramname[0][scen][j]!=paramname[1][scen][k]) k++;
                         x1 = sx[1][scen][j]/(double)n[1][scen][j];
                         dx = fabs(x0-x1)/x0;
                         x3 = sx2[1][scen][j]/(double)n[1][scen][j];
                         dx2 = fabs(x2-x3)/x2;
                         if ((dx>=0.1)or(dx2>=0.1))  star="**";
                         else if ((dx>=0.01)or(dx2>=0.01)) star="* ";
                         else star="  ";
                         printf("%10.8f   %10.8f  %s (%s)\n",dx,dx2,star.c_str(),paramname[1][scen][j].c_str());
                         
                     }
                     cout<<"\n";
                     for (int j=0;j<rt.nstat;j++) {
                         x0 = sx[0][scen][j+rt.nparam[scen]]/(double)n[0][scen][j+rt.nparam[scen]];
                         x1 = sx[1][scen][j+rt.nparam[scen]]/(double)n[1][scen][j+rt.nparam[scen]];
                         dx = fabs(x0-x1)/x0;
                         x0 = sx2[0][scen][j+rt.nparam[scen]]/(double)n[0][scen][j+rt.nparam[scen]];
                         x1 = sx2[1][scen][j+rt.nparam[scen]]/(double)n[1][scen][j+rt.nparam[scen]];
                         dx2 = fabs(x0-x1)/x0;
                         if ((dx>=0.1)or(dx2>=0.1))  star="**";
                         else if ((dx>=0.01)or(dx2>=0.01)) star="* ";
                         else star="  ";
                         printf("%10.8f   %10.8f  %s (%s)\n",dx,dx2,star.c_str(),header.statname[j].c_str());
                     }
                 }
                   
  }
	duree=walltime(&debut);
    fprintf(stdout,"durée = %.2f secondes (%.6f)\n",duree,time_file);
     //fprintf(stdout,"durée dans le remplissage de matC = %.2f secondes\n",time_matC);
     //fprintf(stdout,"durée dans call_polytom = %.2f secondes\n",time_call);
     //fprintf(stdout,"durée dans la lecture de la reftable et le tri des enregistrements = %.2f secondes\n",time_readfile);
	return 0;
};
