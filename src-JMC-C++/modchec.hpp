/*
 * modchec.h
 *
 *  Created on: 13 déc. 2011
 *      Author: ppudlo
 */

#pragma once

#include <string>
#include <vector>

//bool resetstats(std::string s);
std::string pseudoprior2(long double x);
int detphistarOK(int nsel, long double** phistar);
void call_loc(int npart, int nrec, int nsel, long double** ss, float* stat_obs);
void call_acp(int nr, int ns, int nstat, int* numscen, long double** ssref,
              long double** ssphistar, float* stat_obs);
void domodchec(std::string opt, int seed);

int compatphistar(std::vector<int>& numcompat);
