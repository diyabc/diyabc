#include <string>
#include <vector>
#include <cstdlib>
#include <stdexcept>
#include <iostream>
#include <sstream>
#include <iomanip>

#include <string.h>
#include <iso646.h>
#include <sys/stat.h>

#ifdef _MSC_VER
#include "../wingetopt/src/getopt.h"
#else
#include <unistd.h>
#endif

#ifdef _OPENMP
#include <omp.h>
#else
typedef int omp_int_t;

inline omp_int_t omp_get_thread_num() {
	return 0;
}

inline omp_int_t omp_get_max_threads() {
	return 1;
}

inline omp_int_t omp_get_num_threads() {
	return 1;
}

inline void omp_set_num_threads(int) {}
#endif

extern "C" {
#include "../dcmt0.6.1/include/dc.h"
}

#include "acploc.hpp"
#include "header.hpp"
#include "particleset.hpp"
#include "reftable.hpp"
#include "estimparam.hpp"
#include "comparscen.hpp"
#include "bias.hpp"
#include "conf.hpp"
#include "simfile.hpp"
#include "modchec.hpp"
#include "randomgenerator.hpp"
#include "data.hpp"
#include "history.hpp"
#include "randforest.hpp"
#include "mesutils.hpp"


extern "C" void __libc_freeres(void);

#define NSTAT 47

DataC dataobs, datasim;
vector<LocusGroupC> groupe;
vector<ScenarioC> scenario;


string* stat_type;
int* stat_num;

ofstream fprog;
ofstream fpar;

void initstat_typenum() {
	string stat_type0[NSTAT] = {"PID","NAL","HET","VAR","MGW","N2P","H2P","V2P","FST","LIK","DAS","DM2","AML","NHA","NSS","MPD","VPD","DTA","PSS","MNS","VNS","NH2","NS2","MP2","MPB","HST","SML","HP0","HM1","HV1","HMO","NP0","NM1","NV1","NMO","FP0","FM1","FV1","FMO","AP0","AM1","AV1","AMO","RP0","RM1","RV1","RMO"};
	int stat_num0[NSTAT] = {0 , 1 , 2 , 3 , 4 , 5 , 6 , 7 , 8 , 9 , 10 , 11 , 12 , -1 , -2 , -3 , -4 , -5 , -6 , -7 , -8 , -9 , -10 , -11 , -12 , -13 , -14 , 21 , 22 , 23 , 24 , 25 , 26 , 27 , 28 , 29 , 30 , 31 , 32 , 33 , 34 , 35 , 36 , 37 , 38 , 39 , 40};
	/*  Numérotation des stat
	 * 	1 : nal			-1 : nha			 21 : moyenne des genic diversities
	 *  2 : het			-2 : nss             22 : variance des genic diversities
	 *  3 : var			-3 : mpd			 23 : premier quartile des genic diversities
	 *  4 : MGW			-4 : vpd			 24 : troisième quartile des genic diversities
	 *  5 : Fst			-5 : dta			 25 : moyenne des distances de Nei
	 *  6 : lik			-6 : pss			 26 : variance des distances de Nei
	 *  7 : dm2			-7 : mns			 27 : premier quartile des distances de Nei
	 *  8 : n2P			-8 : vns			 28 : troisième quartile des distances de Nei
	 *  9 : h2P			-9 : nh2			 29 : moyenne des distances Fst
	 * 10 : v2P		   -10 : ns2			 30 : variance des distances Fst
	 * 11 : das        -11 : mp2			 31 : premier quartile des distances Fst
	 * 12 : Aml        -12 : mpb			 32 : troisième quartile des distances Fst
	 *      		   -13 : fst			 33 : moyenne des estimations d'admixture
	 * 				   -14 : aml			 34 : variance des estimations d'admixture
	 * 				            			 35 : premier quartile des estimations d'admixture
	 * 										 36 : troisième quartile des estimations d'admixture
	 *										 37 : moyenne des estimations de F3
	 *										 38 : variance des estimations de F3
	 *										 39 : premier quartile des estimations de F3
	 *										 40 : dernier quartile des estimations de F3 
	 * 
	 */
	stat_type = new string[NSTAT];
	for (int i = 0; i < NSTAT; ++i)
		stat_type[i] = stat_type0[i];

	stat_num = new int[NSTAT];
	for (int i = 0; i < NSTAT; ++i)
		stat_num[i] = stat_num0[i];
}


/* Début: pour le nouveau générateur de nombre aléatoires */

mt_struct* r;
#pragma omp threadprivate(r)
mt_struct** mtss;
int countRNG;

void initRNG(int seed) {
	cout << "Debut initialisation RNG." << endl;
	int NB_THREADS;
#pragma omp parallel
	{
		NB_THREADS = omp_get_num_threads();
	}
	// mt_struct **mtss; global variable
	mtss = get_mt_parameters_st(32, 521, 0, NB_THREADS - 1, 4172, &countRNG);
	for (int idxThread = 0; idxThread < countRNG; idxThread++)
		sgenrand_mt(idxThread + seed, mtss[idxThread]);

#pragma omp parallel
	{
		r = mtss[omp_get_thread_num()];
	}
	cout << "Fin initialisation RNG." << endl;
}

void freeRNG(void) {
	free_mt_struct_array(mtss, countRNG);
}

/* Fin: pour le nouveau générateur de nombre aléatoires */


ReftableC rt;
HeaderC header;
ParticleSetC ps;
struct stat stFileInfo;
enregC *enreg, *enregOK;

//char *headerfilename, *reftablefilename,*datafilename,*statobsfilename, *reftablelogfilename,*path,*ident,*stopfilename, *progressfilename;
//char *;
string headerfilename, headersimfilename, reftablefilename, datafilename, statobsfilename, reftablelogfilename, path, ident, stopfilename, progressfilename, scsufilename;
string reftabscen, paramfilename, statfilename;
bool multithread = false, randomforest = false;
int nrecneeded, nenr = 100, nenrOK, *neOK, *netot;
int debuglevel = 0;
int num_threads = 0;
string sremtime, scurfile;
double duree, debutf, dureef, time_file = 0.0, time_reftable = 0.0, remtime;
clock_t debut, debutr;
ofstream fd;

bool RNG_must_be_saved;

/**
 * lecture du fichier header.txt, calcul des stat_obs et lecture de l'entête de reftable.bin
 */

int readheaders() {
	int k;
	string message;
	if (debuglevel == 1) cout << "lecture des entêtes\n";
	k = header.readHeader(headerfilename);
	if (debuglevel == 1) cout << "apres header.readHeader k=" << k << "\n";
	if (k != 0) {
		stringstream erreur;
		erreur << "Erreur dans readheaders().\n" << header.message;
		throw std::runtime_error(erreur.str()); // exit(1)
	}
	message = header.calstatobs(statobsfilename);
	if (debuglevel == 1) cout << "apres header.calstatobs\n";
	datafilename = header.datafilename;
	if (debuglevel == 1) cout << "datafile name : " << header.datafilename << "\n";
	k = rt.testfile(reftablefilename, nenr);
	if (k == 0) {
		k = rt.readheader(reftablefilename, reftablelogfilename, reftabscen);
		if (debuglevel == 1) cout << "apres rt.readheader k=" << k << "   rt.nparam[0]=" << rt.nparam[0] << "\n";
		rt.sethistparamname(header);
		if (debuglevel == 1) cout << "sethistparamname" << "\n";
	}
	return k;
}

/**
 * lecture du fichier headersim.txt
 */
int readheadersim() {
	int k;
	if (debuglevel == 1) cout << "avant header.readHeadersim    headersimfilename=" << headersimfilename << "\n";
	k = header.readHeadersim(headersimfilename);
	return k;
}


/**
 * 
 */
void analyseRNG(string& modpar) {
	string RNGfilename = modpar;

	ifstream fichier(RNGfilename.c_str(), ios::in | ios::binary);
	if (!fichier.is_open()) {
		stringstream erreur;
		erreur << "File " << RNGfilename << " does not exist.\n"
				<< "I cannot analyse it.\n";
		throw runtime_error(erreur.str());

	}
	fichier.seekg(0);
	int sizeoffile;
	fichier.read((char*)&sizeoffile, sizeof(int));
	fichier.close();

	string name = modpar.substr(0, modpar.find_last_of(".")) + "_cores.txt";
	ofstream outfile(name.c_str());
	outfile << sizeoffile << endl;
	outfile.close();
}


/**
 * programme principal : lecture et interprétation des options et lancement des calculs choisis
 */

int main(int argc, char* argv[]) {
#ifdef ___VALGRIND___
    atexit(__libc_freeres);
#endif
	string RNG_filename;
	bool exception_caught = false;
	bool erreur_scenario = false;
	try {
		cout << "debut\n";
		initstat_typenum();
		RNG_must_be_saved = false;
		bool firsttime;
		int k, seed;
		int optchar;
		int computer_identity = 0; // should be 0 if diyabc runs on a single computer
		char action = 'a';
		bool flagp = false, flagi = false, flags = false, simOK, stoprun = false;
		string message, soptarg, estpar, comppar, confpar, acplpar, biaspar, modpar, rngpar, randforpar;

		debut = clock();
		srand(time(NULL));
		seed = rand() % 1000;
		while ((optchar = getopt(argc, argv, "i:p:z:r:e:s:b:c:qkf:g:d:hmqj:a:t:n:w:xyl:o:R:F:Q")) != -1) {
			if (optarg != NULL) soptarg = string(optarg);
			switch (optchar) {

				case 'h': cout << "USAGE :\n";
					cout << "-p <directory of header.txt>\n";
					cout << "-r <number of required data sets in the reftable>\n";
					cout << "-i <name given to the analysis\n";
					cout << "-g <minimum number of particles simulated in a single bunch (default=100)>\n";
					cout << "-m <multithreaded version of the program\n";
					cout << "-q to merge all reftable_$j.bin \n";
					cout << "-Q to merge all reftableRF_$j.bin \n";
					cout << "-s <seed for the random generator (deprecated!!!)>\n"; // TODO: à changer en identité dans le cluster
					cout << "-t <required number of threads>\n";
					cout << "-w <computer's number if in a cluster (0 by default) >\n";
					cout << "          (each computer in the cluster is numbered between 0 and the maximal number of computers in the cluster.)\n";
					cout << "-l for producing a csv reftable file (reftable.csv) and a scenario file (scenario.txt) used by program rf\n";
					cout << "-x for translating a reftable.bin in reftable.txt\n";
					cout << "-y for translating a reftable.bin in reftable.txt writing summary statistics before parameters\n";
					cout << "\n-z <path/RNGfilename.bin> write the number of streams of the RNG into path/RNGfilename_cores.txt\n";
					cout << "\n-n for INITIALIZATION OF THE PARALLEL RNG'S (with parameters as a string including the following options separated par a semi-colon)\n";
					cout << "           t:<maximal number of the threads (per computers if cluster, 16 by default)>\n";
					cout << "           c:<maximal number of computers in the cluster (1 by default)>\n";
					cout << "           s:<seed of the first RNG (1 by default)>\n";
					cout << "           f:<forcing creation of new RNG's and overriding the old ones>\n";

					cout << "\n-e for ABC PARAMETER ESTIMATION (with parameters as a string including the following options separated par a semi-colon)\n";
					cout << "           s:<chosen scenario[s separated by a comma]>\n";
					cout << "           n:<number of simulated datasets taken from reftable>\n";
					cout << "           m:<number of simulated datasets used for the local regression>\n";
					cout << "           t:<number of the transformation (1,2,3 or 4)>\n";
					cout << "           p:<o for original, c for composite, s for scaled, oc,os for both, ocs for all>\n";

					cout << "\n-c for ABC COMPUTATION OF POSTERIOR PROBALITY OF SCENARIOS (idem)\n";
					cout << "           s:<chosen scenarios separated by a comma>\n";
					cout << "           n:<number of simulated datasets taken from reftable>\n";
					cout << "           d:<number of simulated datasets used in the direct approach>\n";
					cout << "           l:<number of simulated datasets used in the logistic regression>\n";
					cout << "           m:<number of requested logistic regressions>\n";

					cout << "\n-F to run a comparison of scenarios based on random forests(idem)\n";
					cout << "           s:<chosen scenario[s separated by a comma]>\n";
					cout << "           n:<number of simulated datasets taken from reftable>\n";
					cout << "           t:<number of trees (default=500)>\n";
					cout << "           b:<size of the bootstrap sample per tree>\n";
					cout << "           k:<number of summary statistics to draw at each node (defaut=sqrt(number of stat))\n";
					cout << "           d:<if present, add linear discriminant scores to summary statistics\n";
					cout << "           o:<number of most effective summary statistics to be shown (default=30case)\n";

					cout << "\n-b for BIAS/PRECISION COMPUTATIONS (idem)\n";
					cout << "           s:<chosen scenario<\n";
					cout << "           n:<number of simulated datasets taken from reftable>\n";
					cout << "           m:<number of simulated datasets used for the local regression>\n";
					cout << "           t:<number of the transformation (1,2,3 or 4)>\n";
					cout << "           p:<o for original, c for composite, s for scaled, oc,os for both, ocs for all>\n";
					cout << "           d:<number of requested test data sets>\n";
					cout << "           h:<histparameter values/priors (see below)>\n";
					cout << "                histparameter values (separated by a space): <parameter name>=<parameter value>\n";
					cout << "                histparameter priors (separated by a space): <parameter name>=<parameter prior as in header.txt>\n";
					cout << "           u:<mutparameter values/priors for successive groups (see below)> groups are named G1, G2 and separated by a star : G1*G2-...\n";
					cout << "                mutparameter values/priors of a given as a set of 6 values/priors : Gx(vx1,vx2,vx3,vx4,vx5,vx6) with :\n";
					cout << "                vx1=<mean mutation rate/prior for group x>    vx2=<shape value/locus mutation rate prior for group x>\n";
					cout << "                vx3=<mean P value/mean P prior for group x>   vx4=<shape value/locus P prior for group x>\n";
					cout << "                vx5 and vx6 correspond to sni mutation rate.\n";
					cout << "                For a DNA sequence group, replace P and sni by k1 and k2 respectively\n";
					cout << "           po to draw parameters for the pseudo-observed data sets from parameter posterior distributions\n";

					cout << "\n-f for CONFIDENCE IN SCENARIO CHOICE COMPUTATIONS (idem)\n";
					cout << "           s:<chosen scenarios among which to choose>\n";
					cout << "           r:<chosen scenario to be tested<\n";
					cout << "           n:<number of simulated datasets taken from reftable>\n";
					cout << "           a:<number of simulated datasets taken from reftable>\n";
					cout << "           t:<number of requested test data sets>\n";
					cout << "           d:<number of simulated datasets used in the direct approach>\n";
					cout << "           l:<number of simulated datasets used in the logistic regression>\n";
					cout << "           m:<number of requested logistic regressions>\n";
					cout << "           h:<histparameter values/priors (as in bias/precision)>\n";
					cout << "           u:<mutparameter values/priors for successive groups (as in bias/precision)\n";
					cout << "           f:<0 if logistic regression on SS, 1 if logistic regression on FDA components>\n";
					cout << "           po to draw parameters for the pseudo-observed data sets from parameter posterior distributions\n";
					cout << "           z:<number of simulated datasets used for the local regression>\n";
					cout << "           b:<number of simulated datasets used for the prior predictive error>\n";
					cout << "           c:<number of simulated datasets used for the posterior predictive error>\n";


					cout << "\n-d for ABC PRIOR/SCENARIO CHECKING (idem)\n";
					cout << "           a:<p for PCA, l for locate observed, pl for both>\n";

					cout << "\n-j for ABC MODEL CHECKING (idem)\n";
					cout << "           s:<chosen scenario[s separated by a comma]>\n";
					cout << "           n:<number of simulated datasets taken from reftable>\n";
					cout << "           m:<number of simulated datasets used for the local regression>\n";
					cout << "           q:<number of datasets simulated from posterior>\n";
					cout << "           t:<number of the transformation (1,2,3 or 4)>\n";
					cout << "           v:<list of summary stat names separated by a comma (if empty keep those of reftable)>\n";

					cout << "\n-k to SIMULATE DATA FILES\n";
					cout << "\n-o to simulate summary statistics from a text file containing scenario and parameter values\n";
					cout << "\n-R <number of required data sets in a reftable file (reftableRF.bin) with all possible summary statistics> (also produces a corresponding statobsRF.txt)\n";
					action = 'h';
					break;


				case 'a':
					debuglevel = atoi(optarg);
					break;

				case 'w':
					computer_identity = atoi(optarg);
					break;

				case 'i':
					ident = soptarg;
					flagi = true;
					break;

				case 'o':
					paramfilename = soptarg;
					action = 'o';
					break;

				case 'p':
					headerfilename = soptarg + "header.txt";
					headersimfilename = soptarg + "headersim.txt";
					reftablefilename = soptarg + "reftable.bin";
					reftablelogfilename = soptarg + "reftable.log";
					statobsfilename = soptarg + "statobs.txt";
					stopfilename = soptarg + ".stop";
					reftabscen = soptarg + "reftabscen.txt";
					path = soptarg;
					flagp = true;
					if (stat(stopfilename.c_str(), &stFileInfo) == 0) remove(stopfilename.c_str());
					break;

				case 's':
					seed = atoi(optarg);
					flags = true;
					break;

				case 'r':
					nrecneeded = atoi(optarg);
					action = 'r';
					break;

				case 'g':
					nenr = atoi(optarg);
					break;

				case 'm':
					multithread = true;
					break;

				case 'c':
					comppar = soptarg;
					action = 'c';
					break;

				case 'e':
					estpar = soptarg;
					action = 'e';
					break;

				case 'b':
					biaspar = soptarg;
					action = 'b';
					break;

				case 'f':
					confpar = soptarg;
					action = 'f';
					break;

				case 'd':
					acplpar = soptarg;
					action = 'd';
					cout << "option -d " << optarg << "      " << soptarg << "\n";
					break;

				case 'R':
					nrecneeded = atoi(optarg);
					action = 'r';
					randomforest = true;
					reftablefilename = path + "reftableRF.bin";
					statobsfilename = path + "statobsRF.txt";
					cout << "simulating data sets with all summary statistics\n";
					break;

				case 'F':
					randforpar = soptarg;
					randomforest = true;
					reftablefilename = path + "reftableRF.bin";
					statobsfilename = path + "statobsRF.txt";
					action = 'F';
					break;

				case 'k':
					action = 'k';
					break;

				case 'l':
					nrecneeded = atoi(optarg);
					action = 'r';
					randomforest = true;
					cout << "simulating for random forest\n";
					break;

				case 'q':
					action = 'q';
					break;

				case 'Q':
					action = 'Q';
					break;

				case 'j':
					modpar = soptarg;
					action = 'j';
					break;

				case 't':
					num_threads = atoi(optarg);
					multithread = true;
					break;

				case 'n':
					rngpar = soptarg;
					action = 'n';
					break;

				case 'x':
					action = 'x';
					break;
				case 'y':
					action = 'y';
					break;
				case 'z':
					rngpar = soptarg;
					action = 'z';
					cout << "RNGpar : " << rngpar << endl;
					flagp = true;
					break;
			}
		}
		if (not flagp) {
			cout << "option -p is compulsory\n";
			exit(1);
		}
		if (not flagi) {
			if (action == 'e') ident = strdup("estim1");
			if (action == 'c') ident = strdup("compar1");
			if (action == 'b') ident = strdup("bias1");
			if (action == 'f') ident = strdup("conf1");
			if (action == 'd') ident = strdup("pcaloc1");
			if (action == 'j') ident = strdup("modchec1");
			if (action == 'o') ident = strdup("statfile.txt");
			if (action == 'F') ident = strdup("rf1");
		}
		if (not flags) seed = time(NULL); // TODO: remove this
		if (num_threads > 0) omp_set_num_threads(num_threads);
		//cout<<"apres les flags action="<<action<<"\n";

		/* Debut: pour le nouveau RNG      */
		if ((action != 'n') and (action != 'h') and (action != 'a') and (action != 'q') and (action != 'Q') and (action != 'z')) {
			// Je dois lire l'état courant des RNG
			mtss = NULL;

			RNG_filename = path + string("RNG_state_") + convertInt4(computer_identity) + string(".bin");
			ifstream test_file(RNG_filename.c_str(), ios::in);
			//cout<<"avant la lecture du ficher "<<RNG_filename<<"\n";
			if (!test_file.is_open()) {
				stringstream erreur;
				erreur << "File " << RNG_filename << " does not exist.\n"
						<< "Use option -n to create it before doing anything else.\n";
				throw runtime_error(erreur.str());
			}

			// lit le fichier RNG_state
			mtss = loadRNG(countRNG, RNG_filename);
			/*
			 *    	 #pragma omp parallel
			 *    	 {
			 *    		 num_threads = omp_get_num_threads();
		}
		*/
			//cout<<"Dans general countRNG="<<countRNG<<"   num_threads="<<num_threads<<"\n";
			if (countRNG < num_threads) {
				stringstream erreur;
				erreur << "I do not have enough states into " << RNG_filename;
				erreur << " to run " << num_threads << " threads, but only " << countRNG << endl;
				erreur << "Reduce the number of threads, or create new RNGs' states." << endl;
				throw std::runtime_error(erreur.str());
			}
#pragma omp parallel
			{
				r = mtss[omp_get_thread_num()];//cout<<omp_get_thread_num()<<"\n";
			}
			cout << "I have read RNGs' states from file " << RNG_filename << endl;
			RNG_must_be_saved = true;
			//initRNG(seed);
		}
		else { // Je n'ai pas besoin de RNGs
			mtss = NULL;
		}
		/* Fin: pour le nouveau RNG      */

		switch (action) {
			case 'r':
				cout << "debut de l'action r\n";
				k = readheaders();
				cout << "apres readheader   k=" << k << "\n";
				cout << dataobs.title << "\n nloc = " << dataobs.nloc << "   nsample = " << dataobs.nsample << "   ";
				fflush(stdin);
				cout << "k=" << k << "\n";
				if (k == 1) {
					cout << "general k==1\n";
					rt.datapath = datafilename;
					rt.nscen = header.nscenarios;
					rt.nrec = 0;
					int scenariosize = (int)scenario.size();
					rt.nrecscen = vector<int>(scenariosize);
					for (int i = 0; i < scenariosize; i++) rt.nrecscen[i] = 0;
					rt.nparam = vector<int>(scenariosize);
					for (int i = 0; i < scenariosize; i++) rt.nparam[i] = scenario[i].nparamvar;
					rt.nstat = header.nstat;
					cout << "general avant rt.writeheader\n";
					rt.filename = reftablefilename;
					rt.writeheader();
					rt.sethistparamname(header);
				}
				else if (k == 2) {
					throw std::runtime_error("cannot create reftable file\n");
					//cout<<"cannot create reftable file\n"; exit(1);
				}
				cout << "DEBUT  nrecneeded=" << nrecneeded << "   rt.nrec=" << rt.nrec << "    rt.nstat=" << rt.nstat << "   nscenarios=" << scenario.size() << "\n";
				if (nrecneeded > rt.nrec) {
					rt.openfile();
					enreg = new enregC[nenr];
					for (int p = 0; p < nenr; p++) {
						enreg[p].stat = vector<float>(header.nstat);
						//cout<<"enreg.param = new float["<<header.nparamtot+3*header.ngroupes <<"]\n";
						enreg[p].param = vector<float>(header.nparamtot + 3 * header.ngroupes);
						enreg[p].numscen = 1;
					}
					cout << "nparammax=" << header.nparamtot + 3 * header.ngroupes << "\n";
					firsttime = true;
					stoprun = false;
					debutr = clock();
					if (not header.drawuntil) {
						neOK = new int[header.nscenarios];
						netot = new int[header.nscenarios];
						for (int p = 0; p < header.nscenarios; p++) {
							neOK[p] = 0;
							netot[p] = 0;
						}
						scsufilename = path + "scenariosuccess.txt";
					}
					while ((not stoprun)and (nrecneeded > rt.nrec)) {
						//cout<<"avant dosimultabref rt.nrec="<<rt.nrec<<"    nenr="<<nenr<<"   nscenarios="<<scenario.size()<<"\n";
						ps.dosimultabref(nenr, false, multithread, firsttime, 0, seed, 0);
						//cout<<"retour de dosimultabref header.drawuntil="<<header.drawuntil<<"\n";

						if (header.drawuntil) {
							simOK = true;
							for (int i = 0; i < nenr; i++)
								if (enreg[i].message != "OK") {
									simOK = false;
									message = enreg[i].message;
								}
							//cout <<"simOK="<<simOK<<"\n";
							if (simOK) {
								//cout<<"simOK=true   nenr="<<nenr<<"\n";
								//debutf=walltime(&clock_zero);
								rt.writerecords(nenr, enreg);
								saveRNG(mtss, countRNG, RNG_filename); // Ajout Pierre le 30/11/2012
								//dureef=walltime(&debutf);time_file += dureef;
								rt.nrec += nenr;
								cout << rt.nrec;
								//if (firsttime) writecourant();
								//cout<<"à la place de writecourant\n";
								ofstream f1(reftablelogfilename.c_str(), ios::out);
								f1 << "OK\n" << rt.nrec << "\n" << TimeToStr(remtime) << "\n";
								f1.close();
								if (((rt.nrec % 1000) == 0)and (rt.nrec < nrecneeded))cout << "   (" << TimeToStr(remtime) << ")""\n";
								else cout << "\n";
								stoprun = (stat(stopfilename.c_str(), &stFileInfo) == 0);
								if (stoprun) remove(stopfilename.c_str());
							}
							else {
								fprog.open(progressfilename.c_str());
								fprog << message << "\n";
								fprog.close();
								cout << "\n\n\n" << message << "\n";
								stoprun = true;
								erreur_scenario = true;
							}
						}
						else {
							nenrOK = 0;
							for (int i = 0; i < nenr; i++) {
								if (enreg[i].message == "OK") {
									nenrOK++;
									neOK[enreg[i].numscen - 1]++;
								}
								netot[enreg[i].numscen - 1]++;
							}
							fd.open(scsufilename.c_str());
							fd << "Numbers of parameter combinations complying with all conditions\n";
							fd << "scenario      tested      successful\n";
							for (int p = 0; p < header.nscenarios; p++)
								fd << setiosflags(ios::fixed) << setw(5) << p + 1 << setw(15) << netot[p] << setw(15) << neOK[p] << "\n";
							fd.close();
							cout << "nenrOK=" << nenrOK << "\n";
							if (nenrOK > 0) {
								enregOK = new enregC[nenrOK];
								for (int p = 0; p < nenrOK; p++) {
									enregOK[p].stat = vector<float>(header.nstat);
									enregOK[p].param = vector<float>(header.nparamtot + 3 * header.ngroupes);
									enregOK[p].numscen = 1;
								}
								nenrOK = 0;
								for (int i = 0; i < nenr; i++) {
									if (enreg[i].message == "OK") {
										for (int p = 0; p < header.nstat; p++) enregOK[nenrOK].stat[p] = enreg[i].stat[p];
										for (int p = 0; p < header.nparamtot + 3 * header.ngroupes; p++) enregOK[nenrOK].param[p] = enreg[i].param[p];
										enregOK[nenrOK].numscen = enreg[i].numscen;
										enregOK[nenrOK].message = "OK";
										nenrOK++;
									}
								}
								rt.writerecords(nenrOK, enregOK);
								cout << "apres rt.writerecords\n";
								saveRNG(mtss, countRNG, RNG_filename); // Ajout Pierre le 30/11/2012
								rt.nrec += nenrOK;
								ofstream f1(reftablelogfilename.c_str(), ios::out);
								f1 << "OK\n" << rt.nrec << "\n" << TimeToStr(remtime) << "\n";
								f1.close();
								if (((rt.nrec % 1000) == 0)and (rt.nrec < nrecneeded))cout << "   (" << TimeToStr(remtime) << ")""\n";
								else cout << "\n";
								stoprun = (stat(stopfilename.c_str(), &stFileInfo) == 0);
								if (stoprun) remove(stopfilename.c_str());
								for (int i = 0; i < nenrOK; i++) {
									enregOK[i].param.clear();
									enregOK[i].stat.clear();
								}
								delete [] enregOK;
							}
						}
						if (firsttime) firsttime = false;
						//if (stoprun) cout<<"STOPRUN=TRUE\n";
					}
					//cout<<"fin du while\n";
					for (int i = 0; i < nenr; i++) {
						if (not enreg[i].param.empty()) enreg[i].param.clear();
						if (not enreg[i].stat.empty()) enreg[i].stat.clear();
					}
					//cout<<"avant delete [] enreg\n";
					delete [] enreg;
					//cout<<"apres delete [] enreg\n";
					//ps.libere(nenr);
					rt.closefile();
					if (randomforest) {
						cout << "ecriture du fichier csv  rt.nstat = " << rt.nstat << "\n";
						rt.bintocsv(header);
					}
					//cout<<"apres rt.closefile\n";
					if (nrecneeded == rt.nrec) {
						ofstream f1(reftablelogfilename.c_str(), ios::out);
						f1 << "END\n" << rt.nrec << "\n";
						f1.close();
					}
					//cout<<"avant header.libere\n";
					//header.libere();
					//cout<<"apres header.libere\n";
					//exit(1);
				}
				else {
					ofstream f1(reftablelogfilename.c_str(), ios::out);
					f1 << "END\n\n";
					f1.close();
				}
				//cout<<"avant le break\n";
				break;
			case 'n':
				doinitRNG(rngpar);
				break;

			case 'e': k = readheaders();
				if (k == 1) {
					cout << "no file reftable.bin in the current directory\n";
					exit(1);
				}
				doestim(estpar);
				break;

			case 'c': k = readheaders();
				if (k == 1) {
					cout << "no file reftable.bin in the current directory\n";
					exit(1);
				}
				docompscen(comppar);
				break;

			case 'b': k = readheaders();
				if (k == 1) {
					cout << "no file reftable.bin in the current directory\n";
					exit(1);
				}

				dobias(biaspar, seed);
				break;

			case 'f': k = readheaders();
				if (k == 1) {
					cout << "no file reftable.bin in the current directory\n";
					exit(1);
				}
				doconf(confpar, seed);
				break;

			case 'F': k = readheaders();
				if (k == 1) {
					cout << "no file reftableRF.bin in the current directory\n";
					exit(1);
				}
				dorandfor(randforpar, seed);
				break;

			case 'k': k = readheadersim();
				dosimfile(seed);
				break;

			case 'd': k = readheaders();
				if (k == 1) {
					cout << "no file reftable.bin in the current directory\n";
					exit(1);
				}
				cout << "avant doacpl soptarg=" << soptarg << "\n";
				doacpl(acplpar);
				break;

			case 'j': k = readheaders();
				if (k == 1) {
					cout << "no file reftable.bin in the current directory\n";
					exit(1);
				}
				nenr = 10;
				domodchec(modpar, seed);
				break;

			case 'q': //header.readHeader(headerfilename);
				k = rt.readheader(reftablefilename, reftablelogfilename, reftabscen);
				rt.concat();
				break;

			case 'Q': //header.readHeader(headerfilename);
				randomforest = true;
				reftablefilename = path + "reftableRF.bin";
				statobsfilename = path + "statobsRF.txt";
				k = rt.readheader(reftablefilename, reftablelogfilename, reftabscen);
				rt.concat2();
				break;
			case 'x':
				k = readheaders();
				if (k == 1) {
					cout << "no file reftable.bin in the current directory\n";
					exit(1);
				}
				k = rt.readheader(reftablefilename, reftablelogfilename, reftabscen);
				rt.bintotxt();
				break;
			case 'y':
				k = readheaders();
				if (k == 1) {
					cout << "no file reftable.bin in the current directory\n";
					exit(1);
				}
				k = rt.readheader(reftablefilename, reftablelogfilename, reftabscen);
				rt.bintotxt2();
				break;
			case 'z':
				analyseRNG(rngpar);
				break;

			case 'o': k = header.readHeader(headerfilename);
				if (debuglevel == 1) cout << "apres header.readHeader k=" << k << "\n";
				if (k != 0) {
					stringstream erreur;
					erreur << "Erreur dans readheaders().\n" << header.message;
					throw std::runtime_error(erreur.str()); // exit(1)
				}
				if (k == 0) {
					rt.sethistparamname(header);
					if (debuglevel == 1) cout << "sethistparamname" << "\n";
				}
				paramfilename = path + paramfilename;
				statfilename = path + ident;
				dosimstat(seed);
				break;

		}
		//cout<<"avant les delete []\n";
		if (rt.mutparam != NULL) {
			delete [] rt.mutparam;
			rt.mutparam = NULL;
		}
		//cout<<"1\n";
		if (not rt.nparam.empty()) rt.nparam.clear();
		//cout<<"2\n";
		if (not rt.nrecscen.empty()) rt.nrecscen.clear();
		//cout<<"3\n";
		/*if (action !='k'){
		 *		for (int iscen=0;iscen<header.nscenarios;iscen++) {
		 *			if (rt.histparam[iscen] !=NULL) delete [] rt.histparam[iscen];
	} 
	}*/
		//cout<<"4\n";
		//if (rt.histparam !=NULL) {delete [] rt.histparam; rt.histparam=NULL;}
		//cout<<"5\n";
		delete [] stat_type;
		//cout<<"6\n";
		delete [] stat_num;
		//cout<<"apres les delete []\n";
		/* Debut: pour le nouveau RNG      */
		// sauvegarde des RNGs' states
		if (RNG_must_be_saved) {
			saveRNG(mtss, countRNG, RNG_filename);
			if (debuglevel < 0)cout << "I have saved current RNGs' states into " << RNG_filename << endl;
			RNG_must_be_saved = false;
		}
		/* Fin: pour le nouveau RNG      */
		duree = walltime(debut);
		if (debuglevel < 0) cout << "durée =" << TimeToStr(duree) << "\n";
		//int aaa;
		//cin>>aaa;
		//fprintf(stdout,"durée = %.2f secondes (%.6f)\n",duree,time_file);
		//fprintf(stdout,"durée dans le remplissage de matC = %.2f secondes\n",time_matC);
		//fprintf(stdout,"durée dans call_polytom = %.2f secondes\n",time_call);
		//fprintf(stdout,"durée dans la lecture de la reftable et le tri des enregistrements = %.2f secondes\n",time_readfile);
	}
	catch (exception& e) {
		cout << "\n\n!!! I have caught an exception which is: \n";
		cout << e.what() << endl << endl;
		exception_caught = true;

	} catch (...) {
		cout << "\n\n!!! I have caught an unknown exception. \n\n";
		exception_caught = true;
	}

	if (RNG_must_be_saved) {
		cout << "I will try to save current RNGs' states.\n";
		saveRNG(mtss, countRNG, RNG_filename);
		cout << "I have saved current RNGs' states into " << RNG_filename << endl;
		RNG_must_be_saved = false;
	}
	freeRNG();
	if ((exception_caught)or (erreur_scenario)) return 1;
	return 0;
};
