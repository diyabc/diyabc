#define BOOST_TEST_MODULE f3reich
#include <boost/test/included/unit_test.hpp>
#include <boost/filesystem.hpp>

#include <string>
#include <vector>


#include "reftable.hpp"
#include "header.hpp"
#include "particleset.hpp"


#include "testvalues.hpp"

#define NSTAT 47

DataC dataobs, datasim;
vector<LocusGroupC> groupe;
vector<ScenarioC> scenario;

string *stat_type;
int *stat_num;

ofstream fprog;

namespace bf = boost::filesystem;

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

int debuglevel = 1;


/* Fin: pour le nouveau générateur de nombre aléatoires */

ReftableC rt;
HeaderC header;
ParticleSetC ps;
struct stat stFileInfo;
vector<enregC> enreg, enregOK;

//char *headerfilename, *reftablefilename,*datafilename,*statobsfilename, *reftablelogfilename,*path,*ident,*stopfilename, *progressfilename;
//char *;
string headersimfilename, reftablefilename, datafilename, statobsfilename, reftablelogfilename, path, ident, stopfilename, progressfilename, scsufilename;
string reftabscen, paramfilename, statfilename;
bool multithread = false, randomforest = false;
int nrecneeded, nenr = 100, nenrOK, *neOK, *netot;
int num_threads = 0;
string sremtime, scurfile;
double duree, debutf, dureef, time_file = 0.0, time_reftable = 0.0, remtime;
clock_t debut, debutr;
ofstream fd;

bool RNG_must_be_saved;

/**
 * lecture du fichier header.txt, calcul des stat_obs et lecture de l'entête de reftable.bin
 */

string headerfilename = "header.txt";

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
  return k;
}

float gap(float a, float b)
{
  return abs(a - b)/(a + b)/2.0;
}

float gap_vec(vector<float> a, vector<float> b)
{
  int n = header.nstat;
  float res = 0.0;
  for (int i = 0; i < n; i++) res += gap(a[i],b[i]);
  return res/(float) n;
}

BOOST_AUTO_TEST_CASE(f3reich_pool_test)
{
  initstat_typenum();
  chdir("pool");
  cout << "debut poolseq\n";
  int res = readheaders();
  float ecart = gap_vec(header.stat_obs,f3reichpool_vals);
  cout << "ecart moyen " << std::setprecision(std::numeric_limits<float>::digits10 + 1) << ecart << endl;
  BOOST_TEST( ecart < 0.01 );

}

BOOST_AUTO_TEST_CASE(f3reich_snp_test)
{
  initstat_typenum();
  chdir("../snp");
  bf::path full_path( bf::current_path() );
  cout << "current directory " << full_path << endl;
  cout << "debut snp\n";
  header = HeaderC();
  int res = readheaders();
  float ecart = gap_vec(header.stat_obs,f3reichsnp_vals);
  cout << "ecart moyen " << std::setprecision(std::numeric_limits<float>::digits10 + 1) << ecart << endl;
  BOOST_TEST( ecart < 0.01 );

}
