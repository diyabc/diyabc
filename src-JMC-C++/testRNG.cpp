/*
 * testRNG.cpp
 *
 *  Created on: 10 déc. 2011
 *      Author: ppudlo
 */

#include <math.h>
#include <cstdlib>
#include <iostream>
#include <omp.h>
using namespace std;
#include "randomgenerator\.hpp"


mt_struct* r;
#pragma omp threadprivate(r)
mt_struct **mtss;
mt_struct **mtss2;
int countRNG;

int nb_threads;

void initRNG (int seed)
{
	cout << "Debut initialisation RNG." << endl;
	int NB_THREADS;
#pragma omp parallel
	{
		NB_THREADS = omp_get_num_threads();
	}
	// mt_struct **mtss; global variable
	mtss = get_mt_parameters_st(32, 521, 0, NB_THREADS-1, 4172, &countRNG);
	for (int idxThread = 0; idxThread < countRNG; idxThread++)
		sgenrand_mt(idxThread+seed, mtss[idxThread]);

#pragma omp parallel
	{
		r = mtss[omp_get_thread_num()];
	}
	cout << "Fin initialisation RNG." << endl;
}
void freeRNG (void)
{
	free_mt_struct_array(mtss, countRNG);
}

int main(){
#pragma omp parallel
	{
		nb_threads = omp_get_num_threads();
	}
    int seed = 1234;
    mtss = NULL;
    initRNG(seed);
    int imax = 1000, imax2 = 100;
    cout << "Nombre de threads = " << nb_threads << endl;
    uint32_t resultat1[nb_threads];
    uint32_t resultat2[nb_threads];
#pragma omp parallel
    {
    // Avant sauvegarde
    	for(int i = 0; i < imax; ++i)
    		genrand_mt(r);
    }
    // sauvegarde
    saveRNG(mtss, countRNG, string("RNG_state.bin"));
    // Recharge
    //freeRNG();
    //cout << "freeRNG OK" << endl;
    mtss2 = loadRNG(countRNG, string("RNG_state.bin"));
#pragma omp parallel
    {
    // Apres sauvegarde
    	for(int i = 0; i < imax2; ++i)
    		genrand_mt(r);
    	resultat1[omp_get_thread_num()] = genrand_mt(r);
    }
#pragma omp parallel
    {
    // Apres recharge
    	r = mtss2[omp_get_thread_num()];
    	for(int i = 0; i < imax2; ++i)
    		genrand_mt(r);
    	resultat2[omp_get_thread_num()] = genrand_mt(r);

    }

    // compare après sauvegarde et après recharge, doit être égal.
    for(int i = 0; i < nb_threads; ++i){
    	cout << i << endl;
    	cout << resultat1[i] << endl;
    	cout << resultat2[i] << endl << endl;
    }

	return 0;
}
