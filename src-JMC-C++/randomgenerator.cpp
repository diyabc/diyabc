/*
 * randomgenerator.cpp
 *
 *  Created on: 23 sept. 2010
 *      Author: cornuet
 */

#include <string>
#include <vector>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <time.h>
#include <math.h>
#include <stdlib.h>

#include "mesutils.hpp"
#include "randomgenerator.hpp"

extern"C" {
#include "../dcmt0.6.1/include/dc.h"
}

using namespace std;

double sqr(double x) {
    return x * x;
}

#if defined (_MSC_VER) && !defined (__INTEL_COMPILER)
extern __declspec(thread) mt_struct* r;
#else
extern mt_struct* r;
#pragma omp threadprivate(r)
#endif

extern string path;

void MwcGen::randinit(unsigned int indice, unsigned int indice2) {
    //              this->mult = multipli[indice % 10000];
    //        srand(indice2);
    //        this->mwcx = rand() % MODULO;
    //              srand(indice2*indice);
    //        this->carry= rand() % MODULO;
    //        //std::cout << "mwcx = " << this->mwcx << "   carry = " << this->carry << "\n";
}

void MwcGen::ecris() {
    //              std::cout<<"mwcx ="<<this->mwcx<<"    mult="<<this->mult<<"   carry="<<this->carry<<"\n";
    //
}

double MwcGen::random() {
    //              unsigned long int y ;
    //              double r;
    //              y = (this->mult)*(this->mwcx) + (this->carry);
    //              this->mwcx = y % MODULO;
    //              this->carry = y / MODULO;
    //              r = (double) this->mwcx /(double)MODULO;
    //              return r;
    return genrand_mt(r) / 4294967296.0 ;
}

int MwcGen::rand0(int max) {
    double r;
    r = floor(max * this->random());
    return (int) r;
}

int MwcGen::rand1(int max) {
    double r;
    r = floor(max * this->random()) + 1.0;
    return (int) r;
}

int* MwcGen::randperm(int n) {
    int *index, j;
    index = new int[n];
    index[0] = 0;
    for (int i = 1; i < n; i++) {
        j = this->rand0(i + 1);
        index[i] = index[j];
        index[j] = i;
    }
    return index;
}

int MwcGen::poisson(double lambda) {
    if (lambda <= 0) {
        return 0;
    }
    int i = 0;
    double p = 0.0;
    double ra;
    while (p < lambda) {
        i++;
        do {
            ra = this->random();
        } while (ra == 0.0);
        p -= log(ra);
    }
    return (i - 1);
}

int MwcGen::binom(int n, double p) { // version light valable pour n<=30
    int k = 0;
    if (n < 1) return 0;
    for (int i = 0; i < n; i++) if (this->random() < p) k++;
    return k;
}

double MwcGen::gunif(double min, double max) {
    double u;
    u = min + (max - min) * this->random();
    return u;
}

double MwcGen::glogunif(double min, double max) {
    double u, logmin, logmax;
    logmin = log(min);
    logmax = log(max);
    u = logmin + (logmax - logmin) * this->random();
    //cout <<logmin<<"   "<<logmax<<"   "<<u<<"   "<<exp(u)<<"\n";
    return exp(u);
}

double MwcGen::gnormal() {
    const double TWOPI = 6.283185307; //8.0*atan(1.0);
    double u1, u2, c, d;
    do {
        u1 = this->random();
    } while (u1 < 0.0000001);
    u2 = this->random();
    c = sqrt(-2.0 * log(u1));
    d = TWOPI * u2;
    if (this->random() < 0.5) return c * cos(d);
    else return c * sin(d);
}

double MwcGen::gnorm(double mean, double sd) {
    return mean + sd * this->gnormal();
}

double MwcGen::glognorm(double mean, double sd) {
    double u, logmean, logsd;
    logmean = log(mean);
    logsd = log(sd);
    u = gnorm(logmean, logsd);
    return exp(u);
}

double MwcGen::ggamma(double shape, double rate) {
    const double log4 = 1.3862943611198906188344642429164, sg_magic = 2.5040773967762740733732583523869, e = 2.7182818284590452353602874713527;
    double alpha, beta, ainv, bbb, ccc, u1, u2, v, x, z, r, res = 0.0;
    bool fin = false;;
    //cout<<"debut de ggamma\n";
    alpha = shape;
    beta = 1.0 / rate;
    //cout<<"alpha = "<<alpha<<"   beta = "<<beta<<"   rate = "<<rate<<"\n";
    if (alpha > 1.0) {
        ainv = sqrt(2.0 * alpha - 1.0);
        bbb = alpha - log4;
        ccc = alpha + ainv;
        while (not fin) {
            do {
                u1 = this->random();
            } while ((u1 < 0.0000001)or (u1 > 0.9999999));
            u2 = 1.0 - this->random();
            v = log(u1 / (1.0 - u1)) / ainv;
            x = alpha * exp(v);
            z = u1 * u1 * u2;
            r = bbb + ccc * v - x;
            fin = ((r + sg_magic - 4.5 * z >= 0.0)or (r >= log(z)));
        }
        res = x * beta;
    } else {
        if (alpha == 1.0) {
            do {
                u1 = this->random();
            } while (u1 < 0.0000001);
            res = -log(u1) * beta;
        }
        if (alpha < 1.0) {
            while (not fin) {
                u2 = this->random();
                bbb = (e + alpha) / e;
                v = bbb * u2;
                if (v <= 1.0) x = exp((1.0 / alpha) * log(v));
                else x = -log((bbb - v) / alpha);
                u1 = this->random();
                if (v > 1.0) fin = (u1 <= exp((alpha - 1.0) * log(x)));
                else fin = (u1 <= exp(x));
            }
            res = x * beta;
        }
    }
    return res;
}

double MwcGen::ggamma2(double mean, double sd) {
    double shape, rate;
    if (sd <= 0.0) return -1.0;
    shape = sqr(mean / sd);
    rate = mean / sqr(sd);
    return this->ggamma(shape, rate);
}

double MwcGen::ggamma3(double mean, double shape) {
    double sd;
    if (shape == 0.0) return 1.0;
    sd = mean / sqrt(shape);
    return this->ggamma2(mean, sd);
}

void MwcGen::resample(int nw, int n, vector<int>& index) {
    bool* deja;
    deja = new bool[nw];
    for (int i = 0; i < nw; i++) deja[i] = false;
    for (int i = 0; i < n; i++) {
        do {
            index[i] = this->rand0(nw);
        } while (deja[index[i]]);
        deja[index[i]] = true;
    }
    delete []deja;
}

void MwcGen::samplewith(int nw, int n, vector<int>& index) {
    for (int i = 0; i < n; i++) index[i] = this->rand0(nw);
}


void write_mt_struct(ofstream& fout, mt_struct* mt) {
    fout.write((char*)&(mt->aaa), sizeof(uint32_t));
    // FIXME: sizeof(int) peut changer d'une plateforme à l'autre !!!
    fout.write((char*)&(mt->mm), sizeof(int));
    fout.write((char*)&(mt->nn), sizeof(int));
    fout.write((char*)&(mt->rr), sizeof(int));
    fout.write((char*)&(mt->ww), sizeof(int));

    fout.write((char*)&(mt->wmask), sizeof(uint32_t));
    fout.write((char*)&(mt->umask), sizeof(uint32_t));
    fout.write((char*)&(mt->lmask), sizeof(uint32_t));

    fout.write((char*)&(mt->shift0), sizeof(int));
    fout.write((char*)&(mt->shift1), sizeof(int));
    fout.write((char*)&(mt->shiftB), sizeof(int));
    fout.write((char*)&(mt->shiftC), sizeof(int));

    fout.write((char*)&(mt->maskB), sizeof(uint32_t));
    fout.write((char*)&(mt->maskC), sizeof(uint32_t));

    fout.write((char*)&(mt->i), sizeof(int));

    fout.write((char*)mt->state, mt->nn * sizeof(uint32_t));
}

void read_mt_struct(ifstream& fichier, mt_struct* mt) {
    fichier.read((char*)&(mt->aaa), sizeof(uint32_t));
    // FIXME: sizeof(int) peut changer d'une plateforme à l'autre !!!
    fichier.read((char*)&(mt->mm), sizeof(int));
    fichier.read((char*)&(mt->nn), sizeof(int));
    fichier.read((char*)&(mt->rr), sizeof(int));
    fichier.read((char*)&(mt->ww), sizeof(int));

    fichier.read((char*)&(mt->wmask), sizeof(uint32_t));
    fichier.read((char*)&(mt->umask), sizeof(uint32_t));
    fichier.read((char*)&(mt->lmask), sizeof(uint32_t));

    fichier.read((char*)&(mt->shift0), sizeof(int));
    fichier.read((char*)&(mt->shift1), sizeof(int));
    fichier.read((char*)&(mt->shiftB), sizeof(int));
    fichier.read((char*)&(mt->shiftC), sizeof(int));

    fichier.read((char*)&(mt->maskB), sizeof(uint32_t));
    fichier.read((char*)&(mt->maskC), sizeof(uint32_t));

    fichier.read((char*)&(mt->i), sizeof(int));

    // dimensionner mt->state, puis le lire
    mt->state = (uint32_t*)malloc(mt->nn * sizeof(uint32_t));
    fichier.read((char*)mt->state, mt->nn * sizeof(uint32_t));
}

int saveRNG(mt_struct** mts, int sizeofmts, string filename) {
    if (mts == NULL) {
        cout << "Pas réussi à sauver le RNG car pointe sur NULL" << endl;
        return 1;
    }
    ofstream fout(filename.c_str(), ios::out | ios::binary);
    if (!fout.is_open()) {
        cout << "Pas réussi à ouvrir le fichier de sauvegarde" << endl;
        return 1;
    }
    fout.write((char*)&sizeofmts, sizeof(int));
    for (int i = 0; i < sizeofmts; ++i) {
        write_mt_struct(fout, mts[i]);
    }
    fout.close();
    return 0;
}

mt_struct** loadRNG(int& sizeofmts, string filename) {
    ifstream fichier(filename.c_str(), ios::in | ios::binary);
    fichier.seekg(0);

    fichier.read((char*)&sizeofmts, sizeof(int));
    // FIXME: comparer le nombre de mts avec le nombre de threads
    // erreur si pas assez !!!
    mt_struct** mts = (mt_struct**) malloc(sizeofmts * sizeof(mt_struct*));
    for (int i = 0; i < sizeofmts; ++i) {
        mts[i] = (mt_struct*)malloc(sizeof(mt_struct));
        read_mt_struct(fichier, mts[i]);
    }
    fichier.close();
    return mts;
}


// converti l'entier number en une string de 4 caractères.
string convertInt4(int number) {
    stringstream ss;//create a stringstream
    ss << number;//add number to the stream
    if ((number < 0) or (number > 9999)) {
        throw std::range_error("You are trying to convert " + ss.str() +
                               " into a string of length 4. Not possible.");
    }
    string ans = ss.str();
    while (ans.size() < 4) {
        ans = "0" + ans;
    }
    return ans;//return a string with the contents of the stream
}


void doinitRNG(string rngpar) {
    string progressfilename = path + "RNG_progress.txt";
    vector<string> ss;
    size_t ns;
    int number_of_threads = 16, number_of_computers = 1;
    time_t seed = time(NULL);
    bool force = false;
    int countRNGs;

    cout << "Beginning of doinitRNG  options : " << rngpar << "\n";
    splitwords(rngpar, ";", ss);
    ns = ss.size();
    for (int i = 0; i < ns; i++) { //cout<<ss[i]<<"\n";
        string s0 = ss[i].substr(0, 2);
        string s1 = ss[i].substr(2);
        if (s0 == "t:") {
            number_of_threads = atoi(s1.c_str());
            cout << "Maximal number of threads per computer is set to " << number_of_threads << endl;
        } else if (s0 == "c:") {
            number_of_computers = atoi(s1.c_str());
            cout << "Maximal number of computers in the cluster is set to " << number_of_computers << endl;
        } else if (s0 == "s:") {
            seed = atoi(s1.c_str());
            cout << "First seed is set to " << seed << endl;
        } else if (s0 == "f:") {
            force = true;
            cout << "Force creation of new RNG's and overwrite existing RNG's states.\n";
        } // fin if
    } // fin for
    cout << "End of doinitRNG options\n";

    // Quitter si RNG_state_0000.bin existe deja
    string firstRNGfile = path + "RNG_state_0000.bin";
    ifstream test_file(firstRNGfile.c_str(), ios::in);
    if ((test_file.is_open()) && (force == false)) {
        cout << "Some files saving the RNG states already exist." << endl
             << "Use f flag if you want to overwrite it." << endl;
        test_file.close();
        exit(1);
    }
    test_file.close();

    // Création des états des RNGs indépendants
    // Il y a autant de RNG que num_nodes * num_threads
    cout << "Beginning of RNGs initializations\n";
    mt_struct** mtss0 = get_mt_parameters_st(32, 607, 0,
                                             (number_of_computers * number_of_threads) - 1, 4172, &countRNGs);
    // Fixons les graines
    for (int i = 0; i < countRNGs; ++i) sgenrand_mt(i + (int)seed, mtss0[i]);
    cout << "End of RNGs initializations\n";


    // Sauvegarde par blocs de num_threads (un fichier par noeud du cluster)
    for (int k = 0; k < number_of_computers; ++k) {
        string RNG_file;
        try {
            RNG_file = path + string("RNG_state_") + convertInt4(k) + string(".bin");
        } catch (std::exception& e) {
            cout << e.what() << endl;
            throw;
        };
        saveRNG(mtss0 + k * (number_of_threads), number_of_threads, RNG_file);
        cout << k << "-th batch of RNGs' states are saved into " << RNG_file << endl;
    }

    // Terminer
    free_mt_struct_array(mtss0, countRNGs);
}

