#include <fstream>
#include <cstring>
#include <iostream>
#include <string>


class ReftableC 
{
public:
    int nrec,*nrecscen,nscen;
    long posnrec;
    char *datapath, *filename, *filelog, *pch;
    int *nparam,nstat,po;
    float *param,*sumstat;
    fstream fifo;
    
    int readheader(char * fname,int nscen,char* datafilename) {
        int nb;
        this->nscen=nscen;
        //cout <<"debut de readheader\n";
        fstream f0(fname,ios::in|ios::out|ios::binary);
        //cout<<"apres fstream\n";
        this->filename = new char[ strlen(fname)+1];
        this->filelog  = new char[ strlen(fname)+1];
        strcpy(this->filename,fname);
        strcpy(this->filelog,fname);
        int p=strcspn(this->filelog,".");
        this->filelog[p]='\0';
        strcat(this->filelog,".log");
        if (!f0) {return 1;}  //fichier non ouvrable e.g. inexistant
        else {
            //cout<<"fichier OK\n";
            f0.seekg(0);
            f0.read((char*)&(this->nrec),sizeof(int));
            //cout <<"nrec = "<<nrec<<"\n";
            this->nrecscen = new int[nscen];
            for (int i=0;i<nscen;i++) {f0.read((char*)&(this->nrecscen[i]),sizeof(int));/*cout<<"nrecscen["<<i<<"] = "<<this->nrecscen[i]<<"\n";*/}
            this->nparam = new int[nscen];
            for (int i=0;i<nscen;i++) {f0.read((char*)&(this->nparam[i]),sizeof(int));/*cout<<"nparam["<<i<<"] = "<<this->nparam[i]<<"\n";*/}
            f0.read((char*)&(this->nstat),sizeof(int));//cout<<"nstat = "<<this->nstat<<"\n";
            f0.close();
            return 0;  //retour normal
        }
    }

    int writeheader() {
        int nb;
        ofstream f1(this->filename,ios::out|ios::binary);
        if (!f1.is_open()) {return 1;}  //fichier impossible à ouvrir
        nb=this->nrec;f1.write((char*)&nb,sizeof(int));
        for (int i=0;i<this->nscen;i++) {nb=this->nrecscen[i];f1.write((char*)&nb,sizeof(int));}
        for (int i=0;i<this->nscen;i++) {nb=this->nparam[i];f1.write((char*)&nb,sizeof(int));}
        nb=this->nstat;f1.write((char*)&nb,sizeof(int));
        f1.close();
        return 0;  //retour normal
    }
    
    enregC readrecord() {
        enregC enr; 
        this->fifo.read((char*)&(enr.numscen),sizeof(int));
        for (int i=0;i<this->nparam[enr.numscen-1];i++) this->fifo.read((char*)&(enr.param[i]),sizeof(float));
        for (int i=0;i<this->nstat;i++) this->fifo.read((char*)&(enr.stat[i]),sizeof(float));
        return enr;
    }
    
    int writerecords(int nenr, enregC *enr) {
        for (int i=0;i<nenr;i++){ 
             if (enr[i].message != "OK") {
                 ofstream f1(this->filelog,ios::out);
                 f1<<enr[i].message<<"\n";
                 f1.close();
                 return 1;
             }
        }
        //cout<<"avant ouverture de reftable.log nscen="<< this->nscen <<"\n";fflush(stdin);
        ofstream f1(this->filelog,ios::out);
        f1<<"OK\n";
        int *nrs,nb;
        nrs = new int[this->nscen];
        for (int i=0;i<this->nscen;i++) nrs[i]=0;
        this->fifo.seekp(0,ios::end);
        //cout<<"apres le seekp\n";fflush(stdin);
        for (int i=0;i<nenr;i++) {
            this->fifo.write((char*)&(enr[i].numscen),sizeof(int));
            nrs[enr[i].numscen-1]++;
            for (int j=0;j<this->nparam[enr[i].numscen-1];j++) this->fifo.write((char*)&(enr[i].param[j]),sizeof(float));
            for (int j=0;j<this->nstat;j++) this->fifo.write((char*)&(enr[i].stat[j]),sizeof(float));
        }
        //cout<<"AVANT MODIF      des nombres d'enregistrements\n";fflush(stdin);
        this->fifo.seekp(0,ios::beg);
        this->fifo.read((char*)&(this->nrec),sizeof(int));
        for (int i=0;i<nscen;i++) this->fifo.read((char*)&(this->nrecscen[i]),sizeof(int));
        this->fifo.seekp(0,ios::beg);
        nb =this->nrec+nenr;this->fifo.write((char*)&nb,sizeof(int));
        for (int i=0;i<nscen;i++) {nb=this->nrecscen[i]+nrs[i];this->fifo.write((char*)&nb,sizeof(int));}
        fifo.flush();
        f1<<this->nrec+nenr<<"\n";
        f1.close();
        return 0;
    }
  
    int openfile() {
        this->fifo.open(this->filename,ios::in|ios::out|ios::binary);
        return 0;
    }
    
    int closefile() {
        this->fifo.close();
        return 0;
    }
  
};