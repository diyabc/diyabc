#define BOOST_TEST_MODULE f3reichPool
#include <boost/test/included/unit_test.hpp>
#include <boost/filesystem.hpp>

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

inline omp_int_t omp_get_thread_num()
{
  return 0;
}

inline omp_int_t omp_get_max_threads()
{
  return 1;
}

inline omp_int_t omp_get_num_threads()
{
  return 1;
}

inline void omp_set_num_threads(int) {}
#endif

extern "C" {
#include "../dcmt/include/dc.h"
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

string *stat_type;
int *stat_num;

ofstream fprog;
ofstream fpar;

void initstat_typenum()
{
  string stat_type0[NSTAT] = {"PID", "NAL", "HET", "VAR", "MGW", "N2P", "H2P", "V2P", "FST", "LIK", "DAS", "DM2", "AML", "NHA", "NSS", "MPD", "VPD", "DTA", "PSS", "MNS", "VNS", "NH2", "NS2", "MP2", "MPB", "HST", "SML", "HP0", "HM1", "HV1", "HMO", "NP0", "NM1", "NV1", "NMO", "FP0", "FM1", "FV1", "FMO", "AP0", "AM1", "AV1", "AMO", "RP0", "RM1", "RV1", "RMO"};
  int stat_num0[NSTAT] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, -1, -2, -3, -4, -5, -6, -7, -8, -9, -10, -11, -12, -13, -14, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40};
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

mt_struct *r;
#pragma omp threadprivate(r)
mt_struct **mtss;
int countRNG;
atomic<int> numloop {0};
atomic<int> rejectedbymrc {0};

void freeRNG(void)
{
  free_mt_struct_array(mtss, countRNG);
}

/* Fin: pour le nouveau générateur de nombre aléatoires */

ReftableC rt;
HeaderC header;
ParticleSetC ps;
struct stat stFileInfo;
vector<enregC> enreg, enregOK;

//char *headerfilename, *reftablefilename,*datafilename,*statobsfilename, *reftablelogfilename,*path,*ident,*stopfilename, *progressfilename;
//char *;
string headerfilename, headersimfilename, reftablefilename, datafilename, statobsfilename, reftablelogfilename, path, ident, stopfilename, progressfilename, scsufilename;
string reftabscen, paramfilename, statfilename;
bool multithread = false, randomforest = false;
int nrecneeded, nenr = 100, nenrOK, *neOK, *netot;
int debuglevel = 1;
int num_threads = 0;
string sremtime, scurfile;
double duree, debutf, dureef, time_file = 0.0, time_reftable = 0.0, remtime;
clock_t debut, debutr;
ofstream fd;

bool RNG_must_be_saved;

/**
 * lecture du fichier header.txt, calcul des stat_obs et lecture de l'entête de reftable.bin
 */

int readheaders()
{
  int k;
  string message;
  if (debuglevel == 1)
    cout << "lecture des entêtes\n";
  k = header.readHeader(headerfilename);
  if (debuglevel == 1)
    cout << "apres header.readHeader k=" << k << "\n";
  if (k != 0)
  {
    stringstream erreur;
    erreur << "Erreur dans readheaders().\n"
           << header.message;
    throw std::runtime_error(erreur.str()); // exit(1)
  }
  message = header.calstatobs(statobsfilename);
  if (debuglevel == 1)
    cout << "apres header.calstatobs\n";
  datafilename = header.datafilename;
  if (debuglevel == 1)
    cout << "datafile name : " << header.datafilename << "\n";
  k = rt.testfile(reftablefilename, nenr);
  if (k == 0)
  {
    k = rt.readheader(reftablefilename, reftablelogfilename, reftabscen);
    if (debuglevel == 1)
      cout << "apres rt.readheader k=" << k << "   rt.nparam[0]=" << rt.nparam[0] << "\n";
    rt.sethistparamname(header);
    if (debuglevel == 1)
      cout << "sethistparamname"
           << "\n";
  }
  return k;
}

/**
 * lecture du fichier headersim.txt
 */
int readheadersim()
{
  int k;
  if (debuglevel == 1)
    cout << "avant header.readHeadersim    headersimfilename=" << headersimfilename << "\n";
  k = header.readHeadersim(headersimfilename);
  return k;
}

/**
 * 
 */
void analyseRNG(string &modpar)
{
  string RNGfilename = modpar;

  ifstream fichier(RNGfilename.c_str(), ios::in | ios::binary);
  if (!fichier.is_open())
  {
    stringstream erreur;
    erreur << "File " << RNGfilename << " does not exist.\n"
           << "I cannot analyse it.\n";
    throw runtime_error(erreur.str());
  }
  fichier.seekg(0);
  int sizeoffile;
  fichier.read((char *)&sizeoffile, sizeof(int));
  fichier.close();

  string name = modpar.substr(0, modpar.find_last_of(".")) + "_cores.txt";
  ofstream outfile(name.c_str());
  outfile << sizeoffile << endl;
  outfile.close();
}

BOOST_AUTO_TEST_CASE(f3reich_pool_test)
/* Compare with void free_test_function() */
{
  headerfilename = "header.txt";
  string RNG_filename;
  bool exception_caught = false;
  bool erreur_scenario = false;
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

  headerfilename = "header.txt";
  headersimfilename = "headersim.txt";
  reftablefilename = "reftable.bin";
  reftablelogfilename = "reftable.log";
  statobsfilename = "statobs.txt";
  stopfilename = ".stop";
  reftabscen = "reftabscen.txt";
  path = "./";
  nrecneeded = 100;
  multithread = true;
  if (num_threads > 0)
    omp_set_num_threads(num_threads);

  mtss = NULL;

  RNG_filename = path + string("RNG_state_") + convertInt4(computer_identity) + string(".bin");
  ifstream test_file(RNG_filename.c_str(), ios::in);
  //cout<<"avant la lecture du ficher "<<RNG_filename<<"\n";
  if (!test_file.is_open())
  {
    stringstream erreur;
    erreur << "File " << RNG_filename << " does not exist.\n"
           << "Use option -n to create it before doing anything else.\n";
    throw runtime_error(erreur.str());
  }

  // lit le fichier RNG_state
  mtss = loadRNG(countRNG, RNG_filename);
  if (countRNG < num_threads)
  {
    stringstream erreur;
    erreur << "I do not have enough states into " << RNG_filename;
    erreur << " to run " << num_threads << " threads, but only " << countRNG << endl;
    erreur << "Reduce the number of threads, or create new RNGs' states." << endl;
    throw std::runtime_error(erreur.str());
  }
#pragma omp parallel
  {
    r = mtss[omp_get_thread_num()]; //cout<<omp_get_thread_num()<<"\n";
  }
  cout << "I have read RNGs' states from file " << RNG_filename << endl;
  RNG_must_be_saved = true;

  boost::filesystem::path full_path(boost::filesystem::current_path());
  std::cout << "Current path is : " << full_path << std::endl;
  int res = readheaders();

  k = readheaders();
  cout << "apres readheader   k=" << k << "\n";
  cout << dataobs.title << "\n nloc = " << dataobs.nloc << "   nsample = " << dataobs.nsample << "   ";
  fflush(stdin);
  cout << "k=" << k << "\n";
  if (k == 1)
  {
    cout << "general k==1\n";
    rt.datapath = datafilename;
    rt.nscen = header.nscenarios;
    rt.nrec = 0;
    int scenariosize = (int)scenario.size();
    rt.nrecscen = vector<int>(scenariosize);
    for (int i = 0; i < scenariosize; i++)
      rt.nrecscen[i] = 0;
    rt.nparam = vector<int>(scenariosize);
    for (int i = 0; i < scenariosize; i++)
      rt.nparam[i] = scenario[i].nparamvar;
    rt.nstat = header.nstat;
    cout << "general avant rt.writeheader\n";
    rt.filename = reftablefilename;
    rt.writeheader();
    rt.sethistparamname(header);
  }
  else if (k == 2)
  {
    throw std::runtime_error("cannot create reftable file\n");
    //cout<<"cannot create reftable file\n"; exit(1);
  }
  cout << "DEBUT  nrecneeded=" << nrecneeded << "   rt.nrec=" << rt.nrec << "    rt.nstat=" << rt.nstat << "   nscenarios=" << scenario.size() << "\n";
  if (nrecneeded > rt.nrec)
  {
    cout << "je rentre là" << endl;
    rt.openfile();
    enreg = vector<enregC>(nenr);
    for (int p = 0; p < nenr; p++)
    {
      enreg[p].stat = vector<float>(header.nstat);
      //cout<<"enreg.param = new float["<<header.nparamtot+3*header.ngroupes <<"]\n";
      enreg[p].param = vector<float>(header.nparamtot + 3 * header.ngroupes);
      enreg[p].numscen = 1;
    }
    cout << "nparammax=" << header.nparamtot + 3 * header.ngroupes << "\n";
    firsttime = true;
    stoprun = false;
    debutr = clock();
    if (not header.drawuntil)
    {
      neOK = new int[header.nscenarios];
      netot = new int[header.nscenarios];
      for (int p = 0; p < header.nscenarios; p++)
      {
        neOK[p] = 0;
        netot[p] = 0;
      }
      scsufilename = path + "scenariosuccess.txt";
    }
    while ((not stoprun) and (nrecneeded > rt.nrec))
    {
      //cout<<"avant dosimultabref rt.nrec="<<rt.nrec<<"    nenr="<<nenr<<"   nscenarios="<<scenario.size()<<"\n";
      ps.dosimultabref(nenr, false, multithread, firsttime, 0, seed, 0);
      //cout<<"retour de dosimultabref header.drawuntil="<<header.drawuntil<<"\n";

      if (header.drawuntil)
      {
        simOK = true;
        for (int i = 0; i < nenr; i++)
          if (enreg[i].message != "OK")
          {
            simOK = false;
            message = enreg[i].message;
          }
        //cout <<"simOK="<<simOK<<"\n";
        if (simOK)
        {
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
          f1 << "OK\n"
             << rt.nrec << "\n"
             << TimeToStr(remtime) << "\n";
          f1.close();
          if (((rt.nrec % 1000) == 0) and (rt.nrec < nrecneeded))
            cout << "   (" << TimeToStr(remtime) << ")"
                                                    "\n";
          else
            cout << "\n";
          stoprun = (stat(stopfilename.c_str(), &stFileInfo) == 0);
          if (stoprun)
            remove(stopfilename.c_str());
        }
        else
        {
          fprog.open(progressfilename.c_str());
          fprog << message << "\n";
          fprog.close();
          cout << "\n\n\n"
               << message << "\n";
          stoprun = true;
          erreur_scenario = true;
        }
      }
      else
      {
        nenrOK = 0;
        for (int i = 0; i < nenr; i++)
        {
          if (enreg[i].message == "OK")
          {
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
        if (nenrOK > 0)
        {
          enregOK = vector<enregC>(nenrOK);
          for (int p = 0; p < nenrOK; p++)
          {
            enregOK[p].stat = vector<float>(header.nstat);
            enregOK[p].param = vector<float>(header.nparamtot + 3 * header.ngroupes);
            enregOK[p].numscen = 1;
          }
          nenrOK = 0;
          for (int i = 0; i < nenr; i++)
          {
            if (enreg[i].message == "OK")
            {
              for (int p = 0; p < header.nstat; p++)
                enregOK[nenrOK].stat[p] = enreg[i].stat[p];
              for (int p = 0; p < header.nparamtot + 3 * header.ngroupes; p++)
                enregOK[nenrOK].param[p] = enreg[i].param[p];
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
          f1 << "OK\n"
             << rt.nrec << "\n"
             << TimeToStr(remtime) << "\n";
          f1.close();
          if (((rt.nrec % 1000) == 0) and (rt.nrec < nrecneeded))
            cout << "   (" << TimeToStr(remtime) << ")"
                                                    "\n";
          else
            cout << "\n";
          stoprun = (stat(stopfilename.c_str(), &stFileInfo) == 0);
          if (stoprun)
            remove(stopfilename.c_str());
          enregOK.clear();
        }
      }
      if (firsttime)
        firsttime = false;
      //if (stoprun) cout<<"STOPRUN=TRUE\n";
    }
    cout << "Nombre de simu avant mrc : " << numloop << endl;
		cout << "Nombre de rejets par mrc : " << rejectedbymrc << endl;
    //cout<<"fin du while\n";
    //cout<<"avant delete [] enreg\n";
    enreg.clear();
    //cout<<"apres delete [] enreg\n";
    //ps.libere(nenr);
    rt.closefile();
    //cout<<"apres rt.closefile\n";
    if (nrecneeded == rt.nrec)
    {
      ofstream f1(reftablelogfilename.c_str(), ios::out);
      f1 << "END\n"
         << rt.nrec << "\n";
      f1.close();
    }
  }
  vector<float> f3reichpool_vals = {0.58700000,
                                    0.62800000,
                                    0.59500000,
                                    0.21477114,
                                    0.22705896,
                                    0.20285954,
                                    0.02720213,
                                    0.02759929,
                                    0.02882971,
                                    0.08870048,
                                    0.08446593,
                                    0.08215811,
                                    0.28200000,
                                    0.25100000,
                                    0.60300000,
                                    0.22319429,
                                    0.20970107,
                                    0.04614485,
                                    0.07134203,
                                    0.06915003,
                                    0.00386333,
                                    0.16025350,
                                    0.15706610,
                                    0.01831950,
                                    0.23700000,
                                    0.20400000,
                                    0.48900000,
                                    0.13069657,
                                    0.12145902,
                                    0.01506385,
                                    0.06529012,
                                    0.06152789,
                                    0.00188179,
                                    0.09972148,
                                    0.09668138,
                                    0.00769763,
                                    0.42386831,
                                    0.59692308,
                                    0.61295419,
                                    0.98545512,
                                    0.70685981,
                                    0.64222233,
                                    0.00906400,
                                    0.12304726,
                                    0.13125621,
                                    0.56775192,
                                    0.28491888,
                                    0.24856946,
                                    0.40600000,
                                    0.80900000,
                                    0.80900000,
                                    0.14325643,
                                    0.03442794,
                                    0.02646005,
                                    0.06093279,
                                    0.00368489,
                                    0.00263723,
                                    0.08504927,
                                    0.00322839,
                                    0.00093312};

  double ecart = 0.0;
  for(int i = 0; i < rt.nstat; i++) ecart += abs(header.stat_obs[i] - f3reichpool_vals[i])/(header.stat_obs[i] + f3reichpool_vals[i])/(double) 2;
  ecart /= (double) rt.nstat;
  cout << "ecart moyen " << std::setprecision(std::numeric_limits<double>::digits10 + 1) << ecart << endl;
  BOOST_TEST( ecart < 0.01 );
}
