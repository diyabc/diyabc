/*
 * initRNG.cpp
 *
 *  Created on: 12 déc. 2011
 *      Author: ppudlo
 */
#include <string>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include "randomgenerator.hpp"
#ifdef _MSC_VER
#include "../wingetopt/src/getopt.h"
#include <iso646.h>
#else
#include <unistd.h>
#endif
#include <time.h>

using namespace std;

// r est inutile ici, mais on en a besoin pour compatibilité avec randomgenerator
mt_struct* r;
#pragma omp threadprivate(r)


mt_struct **mtss;
int countRNG;
int num_threads;
int num_nodes;
string path;
bool multithread;

int main(int argc, char *argv[]){

    bool flagp = false, flags = false, force = false;
    string soptarg;
    int num_nodes = 1;
    int seed = 0;
    int optchar;
    int num_threads = 1;

    // récupération des paramètres de la ligne de commande
    while((optchar = getopt(argc,argv,"hp:c:s:t:f")) !=-1) {
        if (optarg!=NULL) soptarg = string(optarg);
        switch (optchar) {
        case 'h' :
            cout << "Build independent parallel RNGs" << endl << endl;
            cout << "USAGE :\n";
            cout << "-h <print this help>\n";
            cout << "-p <directory of header.txt>\n";
            cout << "-c <number of nodes in the cluster>\n";
            cout << "-s <seed for the random generator>\n";
            cout << "-t <number of threads (per node of the cluster)>\n";
            cout << "-f <force creation of new files if exist>\n";
            return 0;
        case 'p' :
            path = soptarg;
            flagp = true;
            break;
        case 't' :
            num_threads = atoi(optarg);
            break;
        case 'c' :
            num_nodes = atoi(optarg);
            break;
        case 's' :
            seed = atoi(optarg);
            flags = true;
            break;
        case 'f' :
            force = true;
            break;
        }
    }

    // Quitter si pas de répertoire de projet
    if (not flagp) {
        cout << "option -p is compulsory\n";
        exit(1);
    }

    // Quitter si RNG_state_0000.bin existe deja
    string firstRNGfile = path + "RNG_state_0000.bin";
    ifstream test_file(firstRNGfile.c_str(), ios::in);
    if((test_file.is_open()) && (force == false)){
        cout << "Some files saving the RNG states already exist." << endl
             << "Use -f if you want to overwrite it." << endl;
        return 0;
    }
    test_file.close();

    // Choisir une seed si pas donné par la ligne de commande
    if(not flags){
        seed = time(NULL);
    }

    // Création des états des RNGs indépendants
    // Il y a autant de RNG que num_nodes * num_threads
    mtss = get_mt_parameters_st(32, 607, 0, (num_nodes*num_threads)-1, 4172, &countRNG);
    // Fixons les graines
    for (int i = 0; i < countRNG; ++i)
        sgenrand_mt(i + seed, mtss[i]);

    // Test des exceptions
    try{
        string test = convertInt4(-2);
    } catch(std::exception &e) {
        cout << e.what() << endl;
    }

    // Sauvegarde par blocs de num_threads (un fichier par noeud du cluster)
    for(int k = 0; k < num_nodes; ++k){
        string RNG_file;
        try{
            RNG_file = path + string("RNG_state_") + convertInt4(k) + string(".bin");
        } catch(std::exception &e) {
            cout << e.what() << endl;
        };
        saveRNG(mtss + k * (num_threads), num_threads, RNG_file);
    }

    // Terminer
    free_mt_struct_array(mtss, countRNG);

    return 0;
}
