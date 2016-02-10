#include <string>
#include <vector>
#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <algorithm>
#include <cmath>

#ifdef _MSC_VER
#define _AMD64_
#include "../windirent.h"
#else
#include <dirent.h>
#endif

#include "randomgenerator.hpp"
#include "mesutils.hpp"
#include "history.hpp"
#include "header.hpp"
#include "reftable.hpp"

using namespace std;


extern int nrecneeded;
extern double remtime;
extern clock_t debutr;
extern string progressfilename, path;
extern vector<ScenarioC> scenario;


bool operator<(const enregC& lhs, const enregC& rhs)
{
    return lhs.dist < rhs.dist;
}


void ReftableC::sethistparamname(HeaderC const& header)
{
    cout << "debut de sethistparamname\n";
    int nparamvar = 0, pp;
    this->nparamut = header.nparamut;
    cout << "nparamut=" << header.nparamut << "    nscenarios=" << scenario.size() << "\n";
    if (this->nhistparam != NULL)
    {
        delete [] nhistparam;
        nhistparam = NULL;
    }
    this->nhistparam = new int[scenario.size()];
    this->histparam.resize(scenario.size(), vector<HistParameterC>(0));
    this->histparamlength = scenario.size();
    if (this->mutparam != NULL)
    {
        delete [] mutparam;
        mutparam = NULL;
    }
    if (header.nparamut > 0) this->mutparam = new MutParameterC[header.nparamut];
    cout << "avant la boucle des scenarios  nscenarios =" << scenario.size() << "\n";
    int imax = (int)scenario.size();
    for (int i = 0; i < imax; i++)
    {
        nparamvar = 0;
        for (int p = 0; p < scenario[i].nparam; p++)
            if (not scenario[i].histparam[p].prior.constant)
                nparamvar++;
        cout << "scenario " << i << "   scenario[i].nparam=" << scenario[i].nparam << "  nparamvar=" << nparamvar << "\n";
        this->histparam[i].reserve(nparamvar + 2);
        this->nhistparam[i] = nparamvar;
        pp = -1;//cout<<header.scenario[i].nparam<<"\n";
        for (int p = 0; p < scenario[i].nparam; p++)
            if (not scenario[i].histparam[p].prior.constant)
            {
                pp++;
                this->histparam[i].push_back(scenario[i].histparam[p]);
            }
        if (not this->nparam.empty())
        {
            //cout<<"coucou this->nparam[i] = "<<this->nparam[i]<<"   nparamvar="<<nparamvar<<"   header.nparamut="<<header.nparamut<<"\n";
            {
                if (this->nparam[i] != nparamvar + header.nparamut)
                {
                    cout << "PROBLEME scenario " << i << "  nparam=" << this->nparam[i] << "  nparamvar=" << nparamvar << "   nmutparam=" << nparamut << "\n";
                    exit(1);
                }
            }
        }
        //cout<<"coucou2\n";
        if (header.nparamut > 0)
        {
            //cout<<header.mutparam[0].name<<"\n";
            for (int p = 0; p < this->nparamut; p++) this->mutparam[p] = header.mutparam[p];
        }
    }
    cout << "fin de sethistparamname\n";
}

int ReftableC::readheader(string fname, string flogname, string freftabscen)
{
    //cout <<"debut de readheader\n";
    fstream f0(fname.c_str(), ios::in | ios::out | ios::binary);
    this->filename = fname;
    this->filelog = flogname;
    this->filerefscen = freftabscen;
    cout << "dans rt.readheader this->filerefscen = " << this->filerefscen << "\n";
    if (!f0)
    {
        return 1;
    } //fichier non ouvrable e.g. inexistant
    else
    {
        //cout<<"fichier OK\n";
        f0.seekg(0);
        f0.read((char*)&(this->nrec), sizeof(int));
        f0.read((char*)&(this->nscen), sizeof(int));
        //cout <<"readheader.readheader    nscen = "<<this->nscen<<"   nrec="<<this->nrec<<"\n";
        if (not nrecscen.empty()) nrecscen.clear();
        this->nrecscen = vector<int>(this->nscen);
        for (int i = 0; i < this->nscen; i++)
        {
            f0.read((char*)&(this->nrecscen[i]), sizeof(int));/*cout<<"nrecscen["<<i<<"] = "<<this->nrecscen[i]<<"\n";*/
        }
        if (not this->nparam.empty()) this->nparam.clear();
        this->nparam = vector<int>(nscen);
        for (int i = 0; i < this->nscen; i++)
        {
            f0.read((char*)&(this->nparam[i]), sizeof(int));
            cout << "nparam[" << i << "] = " << this->nparam[i] << "\n";
        }
        f0.read((char*)&(this->nstat), sizeof(int));//cout<<"nstat = "<<this->nstat<<"\n";
        f0.close();
        this->nrec0 = this->nrec;
        return 0; //retour normal
    }
}

int ReftableC::writeheader()
{
    int nb;
    ofstream f1(this->filename.c_str(), ios::out | ios::binary);
    if (!f1.is_open())
    {
        throw std::runtime_error("Unable to open " + this->filename + "\n"); //return 1;
    } //fichier impossible à ouvrir
    nb = this->nrec;
    f1.write((char*)&nb, sizeof(int));
    //cout <<"reftable.writeheader     nscen = "<<this->nscen<<"\n";
    nb = this->nscen;
    f1.write((char*)&nb, sizeof(int));
    for (int i = 0; i < this->nscen; i++)
    {
        nb = this->nrecscen[i];
        f1.write((char*)&nb, sizeof(int));
    }
    for (int i = 0; i < this->nscen; i++)
    {
        nb = this->nparam[i];
        f1.write((char*)&nb, sizeof(int));
        cout << "nparam=" << this->nparam[i] << "\n";
    }
    nb = this->nstat;
    f1.write((char*)&nb, sizeof(int));
    f1.close();
    return 0; //retour normal
}

int ReftableC::readrecord(enregC* enr)
{
    int bidon = 0;
    //cout<<"debut de readrecord\n";
    //try {
    if (this->fifo.eof())
    {
        return 1;
    }
    this->fifo.read((char*)&(enr->numscen), sizeof(int));

    if (enr->numscen > this->nscen)
    {
        cout << "\nThe reftable.bin file is corrupted. numscen=" << enr->numscen << "\n";
        FILE* flog;
        flog = fopen(progressfilename.c_str(), "w");
        fprintf(flog, "%s", "The reftable.bin file is corrupted. Clone the project and simulate a new file.\n");
        fclose(flog);
        stringstream erreur;
        erreur << "The reftable.bin file is corrupted. Clone the project and simulate a new file.\n";
        throw std::runtime_error(erreur.str());
        //exit(1);
    }
    //throw(1);
    //cout<<"numscen = "<<enr->numscen<<"\n";
    for (int i = 0; i < this->nparam[enr->numscen - 1]; i++)
    {
        this->fifo.read((char*)&(enr->param[i]), sizeof(float));/*cout<<enr->param[i]<<"  ";*/
    }
    //throw(2);
    //cout <<"\n";
    for (int i = 0; i < this->nstat; i++)
    {
        this->fifo.read((char*)&(enr->stat[i]), sizeof(float));/*cout<<enr->stat[i]<<"  ";if ((i%10)==9) cout<<"\n";*/
    }
    //throw(3);
    /*}
      catch(int e) {
      cout<<"reftable.bin file is corrupted\n";
      if (e==1) cout<<"Impossible to read the scenario number\n";
      if (e==2) cout<<"Impossible to read one parameter value\n";
      if (e==3) cout<<"Impossible to read one summary statistics value\n";
      bidon=e;
      }*/

    //cout <<"\n";
    //cout<<"fin de readrecord\n";
    //cin>>bidon;
    return bidon;
}

int ReftableC::readparam(vector<float>& param)
{
    int numscen;
    float pa;
    this->fifo.read((char*)&(numscen), sizeof(int));
    //      nm=this->nparam[numscen-1];if(nm<nstat) nm=nstat;
    //      pa = new float[nm];
    for (int i = 0; i < this->nparam[numscen - 1]; i++) this->fifo.read((char*)&(param[i]), sizeof(float));
    for (int i = 0; i < this->nstat; i++) this->fifo.read((char*)&(pa), sizeof(float));
    //      delete []pa;
    return numscen;
}

int ReftableC::writerecords(int nenr, enregC* enr)
{
    int bidon;
    for (int i = 0; i < nenr; i++)
    {
        if (enr[i].message != "OK")
        {
            ofstream f1(this->filelog.c_str(), ios::out);
            f1 << enr[i].message << "\n";
            f1.close();
            return 1;
        }
    }
    //cout<<"avant ouverture de reftable.log nscen="<< this->nscen <<"\n";fflush(stdin);
    ofstream f1(this->filelog.c_str(), ios::out);
    f1 << "OK\n";
    int nb;
    int* nrs = new int[this->nscen];
    for (int i = 0; i < this->nscen; i++) nrs[i] = 0;
    //cout<<"avant le seekp\n";fflush(stdin);
    this->fifo.seekp(0, ios::end);
    //cout<<"apres le seekp\n";fflush(stdin);
    for (int i = 0; i < nenr; i++)
    {
        this->fifo.write((char*)&(enr[i].numscen), sizeof(int));
        nrs[enr[i].numscen - 1]++;
        for (int j = 0; j < this->nparam[enr[i].numscen - 1]; j++) this->fifo.write((char*)&(enr[i].param[j]), sizeof(float));
        for (int j = 0; j < this->nstat; j++) this->fifo.write((char*)&(enr[i].stat[j]), sizeof(float));
    }
    //cout<<"AVANT MODIF      des nombres d'enregistrements\n";fflush(stdin);
    this->fifo.seekp(0, ios::beg);
    this->fifo.read((char*)&(this->nrec), sizeof(int));
    //cout <<"reftable.writerecords     avant lecture nscen = "<<this->nscen<<"\n";
    this->fifo.read((char*)&(bidon), sizeof(int));
    //cout <<"reftable.writerecords     apres lecture nscen = "<<bidon<<"\n";
    for (int i = 0; i < this->nscen; i++) this->fifo.read((char*)&(this->nrecscen[i]), sizeof(int));
    this->fifo.seekp(0, ios::beg);
    //cout<< "position avant écriture du nouveau nrec : "<<this->fifo.tellg()<<"\n";
    nb = this->nrec + nenr;
    this->fifo.write((char*)&nb, sizeof(int));
    //cout<< "position avant écriture de this->nscen: "<<this->fifo.tellg()<<"\n";
    nb = this->nscen;
    this->fifo.write((char*)&nb, sizeof(int));
    //cout<< "position avant écriture des this->nrecscen: "<<this->fifo.tellg()<<"\n";

    for (int i = 0; i < this->nscen; i++)
    {
        nb = this->nrecscen[i] + nrs[i];
        this->fifo.write((char*)&nb, sizeof(int));
    }
    fifo.flush();
    remtime = walltime(debutr) * (nrecneeded - this->nrec - nenr) / (this->nrec + nenr - this->nrec0);
    f1 << this->nrec + nenr << "\n";
    f1 << TimeToStr(remtime) << "\n";
    f1.close();
    ofstream f2(this->filerefscen.c_str(), ios::out);
    for (int i = 0; i < this->nscen; i++)
    {
        nb = this->nrecscen[i] + nrs[i];
        f2 << "scenario " << i + 1 << "   " << nb << " records\n";
    }
    f2.close();
    delete []nrs;
    return 0;
}

int ReftableC::openfile()
{
    cout << "\nOuverture du fichier\n" << this->filename << "\n\n";
    this->fifo.open(this->filename.c_str(), ios::in | ios::out | ios::binary);
    return 0;
}

int ReftableC::openfile2()
{
    this->fifo.open(this->filename.c_str(), ios::in | ios::binary);
    this->fifo.seekg(0);
    this->fifo.read((char*)&(this->nrec), sizeof(int));
    this->fifo.read((char*)&(this->nscen), sizeof(int));
    //cout <<"dans openfile2 nrec = "<<nrec<<"\n";
    if (not nrecscen.empty()) nrecscen.clear();
    this->nrecscen = vector<int>(this->nscen);
    for (int i = 0; i < this->nscen; i++)
    {
        this->fifo.read((char*)&(this->nrecscen[i]), sizeof(int));/*cout<<"nrecscen["<<i<<"] = "<<this->nrecscen[i]<<"\n";*/
    }
    if (not this->nparam.empty()) this->nparam.clear();
    this->nparam = vector<int>(nscen);
    for (int i = 0; i < this->nscen; i++)
    {
        this->fifo.read((char*)&(this->nparam[i]), sizeof(int));/*cout<<"nparam["<<i<<"] = "<<this->nparam[i]<<"\n";*/
    }
    this->fifo.read((char*)&(this->nstat), sizeof(int));//cout<<"nstat = "<<this->nstat<<"\n";
    this->nparamax = 0;
    for (int i = 0; i < this->nscen; i++) if (this->nparam[i] > this->nparamax) this->nparamax = this->nparam[i];
    return 0;
}

int ReftableC::testfile(string reftablefilename, int npart)
{
    bool corrompu = false;
    cout << "\nverification de l'integrite de la table de reference \nfichier" << reftablefilename << "\n";
    this->fifo.open(reftablefilename.c_str(), ios::in | ios::out | ios::binary);
    this->fifo.seekg(0);
    this->fifo.read((char*)&(this->nrec), sizeof(int));
    cout << "nrec = " << this->nrec << "\n";
    if (this->nrec < 1)
    {
        this->fifo.close();
        cout << "fichier reftable.bin vide\n\n";
        return 1;
    }
    this->fifo.read((char*)&(this->nscen), sizeof(int));
    if (not this->nrecscen.empty()) nrecscen.clear();
    this->nrecscen = vector<int>(this->nscen);
    for (int i = 0; i < this->nscen; i++)
    {
        this->fifo.read((char*)&(this->nrecscen[i]), sizeof(int));/*cout<<"nrecscen["<<i<<"] = "<<this->nrecscen[i]<<"\n";*/
    }
    if (not this->nparam.empty()) this->nparam.clear();
    this->nparam = vector<int>(nscen);
    for (int i = 0; i < this->nscen; i++)
    {
        this->fifo.read((char*)&(this->nparam[i]), sizeof(int));/*cout<<"nparam["<<i<<"] = "<<this->nparam[i]<<"\n";*/
    }
    this->fifo.read((char*)&(this->nstat), sizeof(int));//cout<<"nstat = "<<this->nstat<<"\n";

    for (int i = 0; i < this->nscen; i++) this->nrecscen[i] = 0;
    int numscen, k, npmax;
    float x;
    cout << this->nscen << " scenario(s)\n";
    for (k = 0; k < this->nrec; k++)
    {
        if ((k % npart) == 0) cout << k + npart << "\r";
        this->fifo.read((char*)&(numscen), sizeof(int));
        if (numscen <= this->nscen)
        {
            this->nrecscen[numscen - 1]++;
            for (int i = 0; i < this->nparam[numscen - 1]; i++) this->fifo.read((char*)&x, sizeof(float));
            for (int i = 0; i < this->nstat; i++) this->fifo.read((char*)&x, sizeof(float));
        }
        else
        {
            corrompu = true;
            break;
        }
    }
    cout << "\n";
    if (not nrecscen.empty()) this->nrecscen.clear();
    if (not this->nparam.empty()) this->nparam.clear();
    if (not corrompu)
    {
        this->fifo.close();
        cout << "fichier reftable OK\n\n";
        return 0;
    }
    enregC e;
    npmax = 0;
    for (int i = 0; i < this->nscen; i++) if (npmax < this->nparam[i]) npmax = this->nparam[i];
    e.param = vector<float>(npmax);
    e.stat = vector<float>(this->nstat);
    int nb;
    if (k > 0) this->nrec = npart * (k / npart);
    cout << "\n\nfichier corrompu à partir de l'enregistrement " << k << "\n";
    cout << "on va récupérer les " << this->nrec << " premiers enregistrements\n";
    this->fifo.seekg(0);
    this->fifo.read((char*)&(k), sizeof(int));
    this->fifo.read((char*)&(this->nscen), sizeof(int));
    for (int i = 0; i < this->nscen; i++)
    {
        this->fifo.read((char*)&(this->nrecscen[i]), sizeof(int));/*cout<<"nrecscen["<<i<<"] = "<<this->nrecscen[i]<<"\n";*/
    }
    if (not this->nparam.empty()) this->nparam.clear();
    this->nparam = vector<int>(nscen);
    for (int i = 0; i < this->nscen; i++)
    {
        this->fifo.read((char*)&(this->nparam[i]), sizeof(int));/*cout<<"nparam["<<i<<"] = "<<this->nparam[i]<<"\n";*/
    }
    this->fifo.read((char*)&(this->nstat), sizeof(int));//cout<<"nstat = "<<this->nstat<<"\n";
    this->filename0 = reftablefilename;

    cout << this->filename0 << "\n";
    this->filename0.replace(reftablefilename.length() - 1, 1, "s");
    cout << this->filename0 << "\n";

    ofstream f1(this->filename0.c_str(), ios::out | ios::binary);
    if (!f1.is_open())
    {
        cout << "impossible d'ouvrir le fichier\n";
        exit(1);
    } //fichier impossible à ouvrir
    nb = this->nrec;
    f1.write((char*)&nb, sizeof(int));
    //cout <<"reftable.writeheader     nscen = "<<this->nscen<<"\n";
    nb = this->nscen;
    f1.write((char*)&nb, sizeof(int));
    for (int i = 0; i < this->nscen; i++)
    {
        nb = this->nrecscen[i];
        f1.write((char*)&nb, sizeof(int));
    }
    for (int i = 0; i < this->nscen; i++)
    {
        nb = this->nparam[i];
        f1.write((char*)&nb, sizeof(int));
    }
    nb = this->nstat;
    f1.write((char*)&nb, sizeof(int));
    cout << "fin d'ecriture de l'entete   nrec=" << this->nrec << "\n";
    for (int h = 0; h < this->nrec; h++)
    {
        this->fifo.read((char*)&(numscen), sizeof(int));
        for (int i = 0; i < this->nparam[numscen - 1]; i++)
        {
            this->fifo.read((char*)&(e.param[i]), sizeof(float));
        }
        for (int i = 0; i < this->nstat; i++)
        {
            this->fifo.read((char*)&(e.stat[i]), sizeof(float));
        }
        f1.write((char*)&(numscen), sizeof(int));
        for (int i = 0; i < this->nparam[numscen - 1]; i++) f1.write((char*)&(e.param[i]), sizeof(float));
        for (int i = 0; i < this->nstat; i++) f1.write((char*)&(e.stat[i]), sizeof(float));
        if ((((h + 1) % npart) == 0)or (h == this->nrec - 1)) cout << h + 1 << "\r";
    }
    if (not this->nparam.empty()) this->nparam.clear();
    e.param.clear();
    e.stat.clear();
    f1.close();
    this->fifo.close();
    remove(reftablefilename.c_str());

    rename(this->filename0.c_str(), reftablefilename.c_str());
    cout << "\nfin d'ecriture des enregistrements\n";
    return 1;
}


int ReftableC::closefile()
{
    int size;
    this->fifo.seekp(0, ios::end);
    size = this->fifo.tellp();
    this->fifo.close();
    cout << "\nfermeture du fichier " << this->filename << " taille=" << size << "\n\n";
    return 0;
}

void ReftableC::bintotxt()
{
    cout << "Starting bintotxt\n";
    enregC enr;
    int bidon, npm;
    fstream f0(this->filename.c_str(), ios::in | ios::binary);
    if (not f0)
    {
        cout << "no file to translate\n";
        exit(1);
    }
    f0.close();
    string nomfi = path + "reftable.txt";
    FILE* f1;
    f1 = fopen(nomfi.c_str(), "w");
    this->openfile2();
    //bidon=this->readheader(this->filename,this->filelog);
    fprintf(f1, "%d simulated data sets\n", this->nrec);
    printf("%d simulated data sets\n", this->nrec);
    fprintf(f1, "%d scenario(s)\n", this->nscen);
    printf("%d scenario(s)\n", this->nscen);
    for (int i = 0; i < this->nscen; i++)
    {
        fprintf(f1, "scenario %d : %d simulated data sets\n", i + 1, this->nrecscen[i]);
        printf("scenario %d : %d simulated data sets\n", i + 1, this->nrecscen[i]);
    }
    for (int i = 0; i < this->nscen; i++)
    {
        fprintf(f1, "scenario %d : %d parameters\n", i + 1, this->nparam[i]);
        printf("scenario %d : %d parameters\n", i + 1, this->nparam[i]);
    }
    fprintf(f1, "%d summary statistics\n", this->nstat);
    printf("%d summary statistics\n", this->nstat);
    npm = 0;
    for (int i = 0; i < this->nscen; i++) if (npm < this->nparam[i]) npm = this->nparam[i];
    enr.param = vector<float>(npm);
    enr.stat = vector<float>(this->nstat);
    while (not fifo.eof())
    {
        bidon = this->readrecord(&enr);
        if (bidon != 0) cout << "probleme à la lecture du reftable\n";
        fprintf(f1, "%3d    ", enr.numscen);//printf("%3d    ",enr.numscen);
        for (int i = 0; i < this->nparam[enr.numscen - 1]; i++)
        {
            //printf(" %12.4f",enr.param[i]);
            //if (enr.param[i]<1) fprintf(f1,"  %6.4e  ",enr.param[i]);
            //else fprintf(f1," %8.0f ",enr.param[i]);
            fprintf(f1, "  %e  ", enr.param[i]);
        }
        fprintf(f1, "    ");
        //for (int i=0;i<this->nstat;i++) fprintf(f1," %10.4f",enr.stat[i]);
        for (int i = 0; i < this->nstat; i++) fprintf(f1, " %12.12e", enr.stat[i]);
        fprintf(f1, "\n");
    }
    this->closefile();
    enr.param.clear();
    enr.stat.clear();
}

void ReftableC::bintotxt2()
{
    cout << "Starting bintotxt2\n";
    enregC enr;
    int bidon, npm;
    fstream f0(this->filename.c_str(), ios::in | ios::binary);
    if (not f0)
    {
        cout << "no file to translate\n";
        exit(1);
    }
    f0.close();
    string nomfi = path + "reftable.txt";
    FILE* f1;
    f1 = fopen(nomfi.c_str(), "w");
    this->openfile2();
    //bidon=this->readheader(this->filename,this->filelog);
    fprintf(f1, "%d simulated data sets\n", this->nrec);
    printf("%d simulated data sets\n", this->nrec);
    fprintf(f1, "%d scenario(s)\n", this->nscen);
    printf("%d scenario(s)\n", this->nscen);
    for (int i = 0; i < this->nscen; i++)
    {
        fprintf(f1, "scenario %d : %d simulated data sets\n", i + 1, this->nrecscen[i]);
        printf("scenario %d : %d simulated data sets\n", i + 1, this->nrecscen[i]);
    }
    for (int i = 0; i < this->nscen; i++)
    {
        fprintf(f1, "scenario %d : %d parameters\n", i + 1, this->nparam[i]);
        printf("scenario %d : %d parameters\n", i + 1, this->nparam[i]);
    }
    fprintf(f1, "%d summary statistics\n", this->nstat);
    printf("%d summary statistics\n", this->nstat);
    npm = 0;
    for (int i = 0; i < this->nscen; i++) if (npm < this->nparam[i]) npm = this->nparam[i];
    enr.param = vector<float>(npm);
    enr.stat = vector<float>(this->nstat);
    while (not fifo.eof())
    {
        bidon = this->readrecord(&enr);
        if (bidon != 0) cout << "Problem when reading reftable\n";
        fprintf(f1, "%3d    ", enr.numscen);
        printf("%3d    ", enr.numscen);
        for (int i = 0; i < this->nstat; i++) fprintf(f1, " %10.4f", enr.stat[i]);
        fprintf(f1, "    ");
        for (int i = 0; i < this->nparam[enr.numscen - 1]; i++)
        {
            printf(" %12.4f", enr.param[i]);
            if (enr.param[i] < 1) fprintf(f1, "  %6.4f  ", enr.param[i]);
            else fprintf(f1, " %8.0f ", enr.param[i]);
        }
        fprintf(f1, "\n");
    }
    this->closefile();
    enr.param.clear();
    enr.stat.clear();
}

void ReftableC::bintocsv(HeaderC const& header)
{
    enregC enr;
    int bidon, npm;
    ofstream f1, f2;
    string nomfi = path + "reftable.csv";
    f1.open(nomfi.c_str());
    nomfi = path + "scenario.txt";
    f2.open(nomfi.c_str());
    this->openfile2();
    cout << "header.nstat=" << header.nstat << "\n";
    for (int i = 0; i < header.nstat; i++)
    {
        f1 << header.statname[i];
        if (i < header.nstat - 1) f1 << ",";
        else f1 << "\n";
    }
    npm = 0;
    for (int i = 0; i < this->nscen; i++) if (npm < this->nparam[i]) npm = this->nparam[i];
    enr.param = vector<float>(npm);
    enr.stat = vector<float>(this->nstat);
    for (int j = 0; j < this->nrec; j++)
    {
        bidon = this->readrecord(&enr);
        if (bidon != 0) cout << "Problem when reading reftable\n";
        f2 << enr.numscen << "\n";
        for (int i = 0; i < header.nstat; i++)
        {
            f1 << enr.stat[i];
            if (i < header.nstat - 1) f1 << ",";
            else f1 << "\n";
        }
    }
    f1.close();
    f2.close();
    this->closefile();
    enr.param.clear();
    enr.stat.clear();
}

void ReftableC::concat()
{
    DIR* dir = opendir(path.c_str());
    struct dirent* entry;
    cout << "Starting the program concat2\n";
    vector<string> reftabname, pathreftabname;
    //vector <string> reflogname;
    string dname;
    int nf, ns, size, nrecc;
    vector<int> nrecscenc;
    while ((entry = readdir(dir)) != NULL)
    {
        dname = entry->d_name;
        if ((dname.find("reftable_") != string::npos)and (dname.find(".bin") != string::npos))
        {
            pathreftabname.push_back(path + entry->d_name);
            reftabname.push_back(entry->d_name);
        };
    }
    nf = (int)reftabname.size();
    cout << "Looking for files with name structure reftable_$$.bin.\n";
    cout << nf << " files found\n";
    if (nf == 0)
    {
        cout << "No files found with name structure reftable_$$.bin. Program terminated.\n";
        exit(1);
    }
    if (nf == 1)
    {
        cout << "Only one file found with name structure reftable_$$.bin. Concatenation not possible. Program terminated.\n";
        exit(1);
    }
    if (nf > 1)
    {
        fstream f1;
        f1.open(this->filename.c_str(), ios::out | ios::binary);
        cout << "Reading the file " << reftabname[0] << "\n";
        fstream f0(pathreftabname[0].c_str(), ios::in | ios::binary);
        f0.seekp(0, ios::end);
        size = f0.tellp();
        f0.seekp(0, ios::beg);
        char* buffer = new char[size];
        f0.read((char*)&(this->nrec), sizeof(int));
        size -= 4;
        f1.write((char*)&(this->nrec), sizeof(int));
        f0.read((char*)&(this->nscen), sizeof(int));
        size -= 4;
        f1.write((char*)&(this->nscen), sizeof(int));
        this->nrecscen = vector<int>(this->nscen);
        nrecscenc = vector<int>(this->nscen);
        for (int i = 0; i < this->nscen; i++)
        {
            f0.read((char*)&(this->nrecscen[i]), sizeof(int));
            size -= 4;
        }
        for (int i = 0; i < this->nscen; i++) f1.write((char*)&(this->nrecscen[i]), sizeof(int));
        f0.read(buffer, size);
        f1.write(buffer, size);
        f0.close();//remove(reftabname[0].c_str());
        delete [] buffer;
        for (int ifich = 1; ifich < nf; ifich++)
        {
            cout << "Reading the file " << reftabname[ifich] << "\n";
            fstream f0(pathreftabname[ifich].c_str(), ios::in | ios::binary);
            f0.seekp(0, ios::end);
            size = f0.tellp();
            f0.seekp(0, ios::beg);
            char* buffer = new char[size];
            f0.read((char*)&(nrecc), sizeof(int));
            this->nrec += nrecc;
            size -= 4;
            f0.read((char*)&(ns), sizeof(int));
            size -= 4;
            for (int i = 0; i < ns; i++)
            {
                f0.read((char*)&(nrecscenc[i]), sizeof(int));
                size -= 4;
                this->nrecscen[i] += nrecscenc[i];
            }
            for (int i = 0; i < ns; i++)
            {
                f0.read((char*)&(nrecscenc[i]), sizeof(int));
                size -= 4;/*cout<<"scenario "<<i<<"   nparam = "<<nrecscenc[i]<<"\n";*/
            }
            f0.read((char*)&ns, sizeof(int));
            size -= 4;/*cout<<"nstat = "<<ns<<"\n";*/
            f0.read(buffer, size);
            f0.close();//remove(reftabname[ifich].c_str());
            f1.seekp(0, ios::end);
            f1.write(buffer, size);
            f1.seekp(0, ios::beg);
            f1.write((char*)&(this->nrec), sizeof(int));
            f1.seekp(8, ios::beg);
            for (int i = 0; i < this->nscen; i++) f1.write((char*)&(this->nrecscen[i]), sizeof(int));
            ofstream f2(this->filelog.c_str(), ios::out);
            f2 << "OK\n";
            f2 << this->nrec << "\n";
            f2.close();
            delete [] buffer;
        }
        f1.close();
        cout << nf << " Files with name structure reftable_$$.bin have been concatenated. The resulting file reftable.bin includes " << this->nrec << " datasets.\n";
    }
    if (not reftabname.empty()) reftabname.clear();
    if (not nrecscenc.empty()) nrecscenc.clear();
}

void ReftableC::concat2()
{
    DIR* dir = opendir(path.c_str());
    struct dirent* entry;
    cout << "Starting the program concat2\n";
    vector<string> reftabname, pathreftabname;
    //vector <string> reflogname;
    string curname, refname = this->filename;
    int nf, ns, size, nrecc;
    vector<int> nrecscenc;
    cout << "refname=" << refname << "\n";
    int pos = refname.find("RF.bin");
    refname = refname.substr(0, pos + 2) + "_";
    cout << "refname=" << refname << "\n";
    while ((entry = readdir(dir)) != NULL)
    {
        curname = path + entry->d_name;
        if (curname.find(refname) != string::npos)
        {
            pathreftabname.push_back(curname);
            reftabname.push_back(entry->d_name);
        };
    }
    nf = (int)reftabname.size();
    cout << "Looking for files with name structure reftable_$$.bin.\n";
    cout << nf << " files found\n";
    if (nf == 0)
    {
        cout << "No file with name structure reftable_$$.bin have been found. Program terminated.\n";
        exit(1);
    }
    if (nf == 1)
    {
        cout << "Only one file with name structure reftable_$$.bin have been found. No concatenation possible. Program terminated.\n";
        exit(1);
    }
    if (nf > 1)
    {
        fstream f1;
        f1.open(this->filename.c_str(), ios::out | ios::binary);
        cout << "Reading file " << reftabname[0] << "\n";
        fstream f0(pathreftabname[0].c_str(), ios::in | ios::binary);
        f0.seekp(0, ios::end);
        size = f0.tellp();
        f0.seekp(0, ios::beg);
        char* buffer = new char[size];
        f0.read((char*)&(this->nrec), sizeof(int));
        size -= 4;
        f1.write((char*)&(this->nrec), sizeof(int));
        f0.read((char*)&(this->nscen), sizeof(int));
        size -= 4;
        f1.write((char*)&(this->nscen), sizeof(int));
        this->nrecscen = vector<int>(this->nscen);
        nrecscenc = vector<int>(this->nscen);
        for (int i = 0; i < this->nscen; i++)
        {
            f0.read((char*)&(this->nrecscen[i]), sizeof(int));
            size -= 4;
        }
        for (int i = 0; i < this->nscen; i++) f1.write((char*)&(this->nrecscen[i]), sizeof(int));
        f0.read(buffer, size);
        f1.write(buffer, size);
        f0.close();//remove(reftabname[0].c_str());
        delete [] buffer;
        for (int ifich = 1; ifich < nf; ifich++)
        {
            cout << "Reading file " << reftabname[ifich] << "\n";
            fstream f0(pathreftabname[ifich].c_str(), ios::in | ios::binary);
            f0.seekp(0, ios::end);
            size = f0.tellp();
            f0.seekp(0, ios::beg);
            char* buffer = new char[size];
            f0.read((char*)&(nrecc), sizeof(int));
            this->nrec += nrecc;
            size -= 4;
            f0.read((char*)&(ns), sizeof(int));
            size -= 4;
            for (int i = 0; i < ns; i++)
            {
                f0.read((char*)&(nrecscenc[i]), sizeof(int));
                size -= 4;
                this->nrecscen[i] += nrecscenc[i];
            }
            for (int i = 0; i < ns; i++)
            {
                f0.read((char*)&(nrecscenc[i]), sizeof(int));
                size -= 4;/*cout<<"scenario "<<i<<"   nparam = "<<nrecscenc[i]<<"\n";*/
            }
            f0.read((char*)&ns, sizeof(int));
            size -= 4;/*cout<<"nstat = "<<ns<<"\n";*/
            f0.read(buffer, size);
            f0.close();//remove(reftabname[ifich].c_str());
            f1.seekp(0, ios::end);
            f1.write(buffer, size);
            f1.seekp(0, ios::beg);
            f1.write((char*)&(this->nrec), sizeof(int));
            f1.seekp(8, ios::beg);
            for (int i = 0; i < this->nscen; i++) f1.write((char*)&(this->nrecscen[i]), sizeof(int));
            ofstream f2(this->filelog.c_str(), ios::out);
            f2 << "OK\n";
            f2 << this->nrec << "\n";
            f2.close();
            delete [] buffer;
        }
        f1.close();
        cout << nf << " Files with name structure reftable_$$.bin have been concatenated. The resulting file reftableRF.bin includes " << this->nrec << " datasets.\n";
    }
    if (not reftabname.empty()) reftabname.clear();
    if (not nrecscenc.empty()) nrecscenc.clear();
}

/**
 * calcule les variances des statistiques résumées
 * sur les 100000 premiers enregistrements de la table de référence
 */
int ReftableC::cal_varstat()
{
    int nrecutil, iscen, nsOK, bidon, i;
    long double *sx, *sx2, x, an, nr, *min, *max;
    bool scenOK;
    //cout <<"debut de cal_varstat   rt.nrec="<<this->nrec<<"   scenchoisi="<<this->scenchoisi[0]<<"\n";
    enregC enr;
    nrecutil = 100000;
    if (nrecutil > this->nrec) nrecutil = this->nrec;
    nr = 0;
    for (int i = 0; i < nscenchoisi; i++) nr += nrecscen[this->scenchoisi[i] - 1];
    cout << "debut de cal_varstat   nr=" << nr << "   nscenchoisi=" << nscenchoisi;
    if (nrecutil > nr) nrecutil = nr;
    cout << "    nrecutil=" << nrecutil << "\n";
    sx = new long double[this->nstat];
    sx2 = new long double[this->nstat];
    min = new long double[this->nstat];
    max = new long double[this->nstat];
    var_stat = new long double[this->nstat];
    for (int j = 0; j < this->nstat; j++)
    {
        sx[j] = 0.0;
        sx2[j] = 0.0;
        min[j] = 10.0;
        max[j] = -10.0;
    }
    enr.stat = vector<float>(this->nstat);
    this->nparamax = 0;
    for (int i = 0; i < this->nscen; i++) if (this->nparam[i] > this->nparamax) this->nparamax = this->nparam[i];
    enr.param = vector<float>(this->nparamax);
    this->openfile2();
    i = 0;
    while (i < nrecutil)
    {
        //cout<<"avant readrecord\n";
        bidon = this->readrecord(&enr);
        if (bidon != 0)cout << "problème dans la lecture du reftable\n";
        //cout<<"coucou\n";
        scenOK = false;
        iscen = 0;
        while ((not scenOK)and (iscen < this->nscenchoisi))
        {
            scenOK = (enr.numscen == this->scenchoisi[iscen]);
            iscen++;
        }
        //if ((scenOK)and(i==0))for (int j=0;j<this->nstat;j++) cout<<enr.stat[j]<<"   ";cout<<"\n";
        if (scenOK)
        {
            i++;
            for (int j = 0; j < this->nstat; j++)
            {
                x = (long double)enr.stat[j];
                sx[j] += x;
                sx2[j] += x * x;
                if (min[j] > x) min[j] = x;
                if (max[j] < x) max[j] = x;
            }
        }
    }
    //cout<<i<<"   "<<nrecutil<<"\n";
    if (i < nrecutil) nrecutil = i;
    //cout<<i<<"   "<<nrecutil<<"\n";
    this->closefile();
    nsOK = 0;
    an = (long double)nrecutil;
    for (int j = 0; j < this->nstat; j++)
    {
        this->var_stat[j] = (sx2[j] - sx[j] * sx[j] / an) / (an - 1.0);
        if (this->var_stat[j] > 1E-20) nsOK++;
        //printf("sx2[%3d] = %12.8Le   sx[%3d] = %12.8Le\n",j,sx2[j],j,sx[j]);
        //printf("var_stat[%3d] = %12.8Lf   min=%12.8Lf   max=%12.8Lf\n",j,this->var_stat[j],min[j],max[j]);
    }
    delete []sx;
    delete []sx2;
    delete []min;
    delete []max;
    cout << "nstatOK = " << nsOK << "\n";
    //      exit(1);
    return nsOK;
}

/**
 * alloue la mémoire pour enrsel
 */
void ReftableC::alloue_enrsel(int nsel)
{
    int nparamax = 0;
    for (int i = 0; i < this->nscen; i++) if (this->nparam[i] > nparamax) nparamax = this->nparam[i];
    //cout<<"alloue_enrsel nsel="<<nsel<<"   nparamax="<<nparamax<<"\n";fflush(stdout);
    this->enrsel = new enregC[2 * nsel];
    for (int i = 0; i < 2 * nsel; i++)
    {
        this->enrsel[i].param = vector<float>(nparamax);
        this->enrsel[i].stat = vector<float>(this->nstat);
    }
}

/**
 * desalloue la mémoire de enrsel
 */
void ReftableC::desalloue_enrsel(int nsel)
{
    for (int i = 0; i < 2 * nsel; i++)
    {
        this->enrsel[i].param.clear();
        this->enrsel[i].stat.clear();
    }
    delete []this->enrsel;
    this->enrsel = NULL;
}

/**
 * calcule la distance de chaque jeu de données simulé au jeu observé
 * et sélectionne les nsel enregistrements les plus proches (copiés dans enregC *enrsel)
 */
void ReftableC::cal_dist(int nrec, int nsel, float* stat_obs, bool scenarioteste, bool allscenarios)
{
    int nn, nparamax, nrecOK = 0, iscen, bidon, step;
    bool firstloop = true, scenOK;
    long double diff, distmin = 1E100;
    this->nreclus = 0;
    step = nrec / 100;
    if (scenarioteste) cout << "cal_dist scenarioteste = " << this->scenteste << "\n";
    nn = nsel;
    //cout<<"nrec="<<nrec<<"   nscenchoisi="<< this->nscenchoisi<<"   nn="<<nn<<"\n";
    nparamax = 0;
    for (int i = 0; i < this->nscen; i++) if (this->nparam[i] > nparamax) nparamax = this->nparam[i];
    //cout<<"cal_dist nsel="<<nsel<<"   nparamax="<<nparamax<<"   nrec="<<nrec<<"   nreclus="<<this->nreclus<<"   nstat="<<this->nstat<<"   2*nn="<<2*nn<<"\n";
    //cout<<" apres allocation de enrsel\n";
    //for (int j=0;j<this->nstat;j++) cout<<"var_stat["<<j<<"]="<<this->var_stat[j]<<"\n";
    //cout<<"\n";
    //for (int j=0;j<this->nstat;j++) cout<<"stat_obs["<<j<<"]="<<stat_obs[j]<<"\n";
    //cout<<"avant openfile2() nrec="<<nrec<<"   nsel="<<nsel<<"\n";
    this->openfile2();
    cout << "apres openfile2()  firsloop=" << firstloop << "   nn=" << nn << "\n";
    while ((this->nreclus < nrec)and (not fifo.eof()))
    {
        if (firstloop)
        {
            nrecOK = 0;
            firstloop = false;
        }
        else nrecOK = nn;
        while ((nrecOK < 2 * nn)and (this->nreclus < nrec))
        {
            do
            {
                bidon = this->readrecord(&(this->enrsel[nrecOK]));
                if (bidon != 0) cout << "bidon=" << bidon << "\n";
            }
            while (bidon != 0);
            if (not allscenarios)
            {
                scenOK = false;
                if (scenarioteste) scenOK = (this->enrsel[nrecOK].numscen == this->scenteste);
                else
                {
                    iscen = 0;
                    while ((not scenOK)and (iscen < this->nscenchoisi))
                    {
                        scenOK = (this->enrsel[nrecOK].numscen == this->scenchoisi[iscen]);
                        //if (scenOK) cout<<" SCENOK pour this->enrsel["<<nrecOK<<"].numscen = "<<this->enrsel[nrecOK].numscen<<"\n";
                        iscen++;
                    }
                }
            }
            else scenOK = true;
            if (scenOK)
            {
                if (not scenarioteste) this->nreclus++; //cout<<"nreclus = "<<nreclus<<"   nrecOK="<<nrecOK<<"\n";
                this->enrsel[nrecOK].dist = 0.0;
                for (int j = 0; j < this->nstat; j++)
                    if (this->var_stat[j] > 1E-20)
                    {
                        diff = (long double)(this->enrsel[nrecOK].stat[j] - stat_obs[j]);
                        this->enrsel[nrecOK].dist += diff * diff / this->var_stat[j];
                        //cout<<"this->enrsel[nrecOK].dist = "<<this->enrsel[nrecOK].dist<<"\n";
                        //if (nreclus<10) printf("  %12.6f   %12.6f   %12.6Lf   %12.8Lf\n",this->enrsel[nrecOK].stat[j],stat_obs[j],diff*diff,this->enrsel[nrecOK].dist);
                    }
                //if (nreclus<10) cout<<"\n";
                this->enrsel[nrecOK].dist = sqrt(this->enrsel[nrecOK].dist);
                if (distmin > this->enrsel[nrecOK].dist)
                {
                    distmin = this->enrsel[nrecOK].dist;/*cout <<"distmin="<<distmin<<"\n";*/
                }
                nrecOK++;
                if (this->nreclus == nrec) break;
            }
            //cout<<"nrecOK="<<nrecOK<<"\n";
            if (scenarioteste) this->nreclus++;
            //cout<<"nreclus = "<<nreclus<<"\n";
            if ((this->nreclus % step) == 0)
            {
                cout << "\rcal_dist : " << this->nreclus / step << "%";
                fflush(stdout);
            }
        }
        sort(&this->enrsel[0], &this->enrsel[2 * nn]);
    }
    cout << "\rcal_dist : fini   nreclus=" << nreclus << "\n";
    this->closefile();
    int nnuls = 0;
    for (int i = 0; i < 2 * nn; i++) if (this->enrsel[i].dist == 0.0) nnuls++;
    //cout<<nnuls<<" ont une distance nulle (2nn="<<2*nn<<")\n";
    //cout<<"\nnrec_lus = "<<this->nreclus<<"   nrecOK = "<<nrecOK<<"\n";
    //cout<<"    distmin = "<<this->enrsel[0].dist<<"    distmax = "<<this->enrsel[nsel-1].dist<<"\n";
    //printf("    distmin = %12.8Lf   distmax = %12.8Lf\n",this->enrsel[0].dist/(long double)this->nstat,this->enrsel[nsel-1].dist/(long double)this->nstat);
    /*for (int i=0;i<nsel;i++){
      cout<<this->enrsel[i].numscen<<"  ";
      for (int j=0;j<this->nparam[this->enrsel[i].numscen-1];j++) printf("  %6.0f",this->enrsel[i].param[j]);
      for (int j=0;j<this->nstat;j++) printf("  %8.5f",this->enrsel[i].stat[j]);
      cout<<"\n";
      }*/
}

