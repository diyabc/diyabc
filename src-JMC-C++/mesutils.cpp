#include <string>
#include <vector>
#include <iostream>
#include <sstream>
#include <algorithm>
#include <stdexcept>
#include <limits>
#include <cstring>
#include <time.h>
#include <math.h>


#include "matrices.hpp"
#include "randomgenerator.hpp"
#include "mesutils.hpp"


using namespace std;

#define PI 3.141592653589793


/**
 * calcule le temps écoulé depuis t0 (récupéré sur internet)
 */
double walltime(clock_t& t0) {
	clock_t t;
	t = clock() - t0;
	return (static_cast<float>(t)) / CLOCKS_PER_SEC;
}

/**
 * renvoit la valeur entière contenue dans la num-ième sous-chaîne de s
 * fonctionne également si la valeur entière est encadrée par des parenthèses ou des crochets
 */
int getwordint(string s, int num) {
	s.append(" ");
	while (s.find(" ") == 0) s = s.substr(1);
	int i = 1, j;
	while ((i < num)and (s.length() > 1)) {
		j = s.find(" ");
		s = s.substr(j);
		while (s.find(" ") == 0) s = s.substr(1);
		i++;
	}
	s = s.substr(0, s.find(" "));
	if ((s.find("(") == 0)or (s.find("[") == 0)) s = s.substr(1, s.length() - 2);
	return atoi(s.c_str());
}

/**
 * renvoit la valeur flottante contenue dans la num-ième sous-chaîne de s
 * fonctionne également si la valeur flottante est encadrée par des parenthèses ou des crochets
 */
double getwordfloat(string s, int num) {
	s.append(" ");
	while (s.find(" ") == 0) s = s.substr(1);
	int i = 1, j;
	while ((i < num)and (s.length() > 1)) {
		j = s.find(" ");
		s = s.substr(j);
		while (s.find(" ") == 0) s = s.substr(1);
		i++;
	}
	s = s.substr(0, s.find(" "));
	if ((s.find("(") == 0)or (s.find("[") == 0)) s = s.substr(1, s.length() - 2);
	return atof(s.c_str());
}

/**
 * découpe la chaîne s en sous-chaînes séparées par le séparateur sep
 * le nombre de sous chaînes est donné par *k
 */
string* splitwords(string s, string sep, int* k) {
	int j = 0, j0;
	while (s.find(sep) == 0) s = s.substr(1);
	*k = 0;
	s.append(sep);
	//cout<<"s="<<s<<"\n";
	string *sb, s0, s1;
	sb = NULL;
	s1 = string();
	for (int i = 0; i < (int)s.length(); i++) {
		s0 = s.substr(i, 1);
		if (s0 == sep) {
			j++;
			if (j == 1) {
				s1.append(s0);
				if (j == 1) (*k)++;
				//cout <<" j=1  k="<<*k<<"\n";
			}
		}
		else {
			s1.append(s0);
			j = 0;
		}
	}
	sb = new string[*k];
	for (int i = 0; i < *k; i++) {
		j0 = s1.find(sep);
		sb[i] = s1.substr(0, j0);
		s1 = s1.substr(j0 + 1, s.length());

	}
	//cout <<"k="<<*k<<"\n";
	return sb;
}

void splitwords(string s, string sep, vector<string>& resultat) {
	int j = 0, j0;
	while (s.find(sep) == 0) s = s.substr(1);
	int k = 0;
	s.append(sep);
	//cout<<"s="<<s<<"\n";
	vector<string> sb(0);
	string s0, s1;
	s1 = string();
	for (int i = 0; i < (int)s.length(); i++) {
		s0 = s.substr(i, 1);
		if (s0 == sep) {
			j++;
			if (j == 1) {
				s1.append(s0);
				if (j == 1) (k)++;
				//cout <<" j=1  k="<<*k<<"\n";
			}
		}
		else {
			s1.append(s0);
			j = 0;
		}
	}
	sb.resize(k);
	for (int i = 0; i < k; i++) {
		j0 = s1.find(sep);
		sb[i] = s1.substr(0, j0);
		s1 = s1.substr(j0 + 1, s.length());
	}
	//cout <<"k="<<*k<<"\n";
	resultat.swap(sb);
}

/**
 * change les tabulations en espaces
 */

string purgetab(string s) {
	string ss = "", s0;
	for (int i = 0; i < (int)s.length(); i++) {
		s0 = s.substr(i, 1);
		if (s0 == "\t") s0 = " ";
		ss += s0;
	}
	//	if (ss.at(ss.length()-1)=='\r')cout<<"HELLO\n";
	return ss;
}

/**
 * découpe la chaîne s en sous-chaînes séparées par le séparateur sep
 * le découpage s'arrête quand le nombre de sous-chaînes atteint m
 */
string* splitwordsR(string s, string sep, int m, int* k) {
	int j = 0, j0;
	while (s.find(sep) == 0) s = s.substr(1);
	*k = 0;
	string *sb, s0, s1;
	if (s.length() == 0) {
		sb = new string[1];
		sb[0] = "";
		return sb;
	}
	s.append(sep);
	s1 = string();
	for (int i = 0; i < (int)s.length(); i++) {
		s0 = s.substr(i, 1);
		if (s0 == sep) {
			j++;
			if (j == 1) {
				s1.append(s0);
				if (j == 1) (*k)++;
				//cout <<" j=1  k="<<*k<<"\n";
			}
		}
		else {
			s1.append(s0);
			j = 0;
		}
		if ((*k) == m) break;
	}
	sb = new string[*k];
	for (int i = 0; i < *k; i++) {
		j0 = s1.find(sep);
		sb[i] = s1.substr(0, j0);
		s1 = s1.substr(j0 + 1, s.length());

	}
	//cout <<"k="<<*k<<"\n";
	return sb;
}

void splitwordsR(string s, string sep, int m, vector<string>& resultat) {
	int j = 0, j0;
	while (s.find(sep) == 0) s = s.substr(1);
	int k = 0;
	string s0, s1;
	if (s.length() == 0) {
		resultat.resize(0);
		return;
	}
	s.append(sep);
	s1 = string();
	for (int i = 0; i < (int)s.length(); i++) {
		s0 = s.substr(i, 1);
		if (s0 == sep) {
			j++;
			if (j == 1) {
				s1.append(s0);
				if (j == 1) k++;
				//cout <<" j=1  k="<<*k<<"\n";
			}
		}
		else {
			s1.append(s0);
			j = 0;
		}
		if (k == m) break;
	}
	vector<string> sb(k);
	for (int i = 0; i < k; i++) {
		j0 = s1.find(sep);
		sb[i] = s1.substr(0, j0);
		s1 = s1.substr(j0 + 1, s.length());

	}
	resultat.swap(sb);
	//cout <<"k="<<*k<<"\n";
}

/**
 * Centre la chaîne s dans une chaîne de longeur totale k
 */
string centre(string const& s, int k) {
	int l = s.length();
	string sb = s;
	int av, ap;
	av = (k - l) / 2;
	ap = k - l - av;
	if ((av <= 0)and (ap <= 0)) return s;
	if (av <= 0) return s + " ";
	sb = s;
	for (int i = 0; i < av; i++) sb = " " + sb;
	if (ap == 0) return sb;
	for (int i = 0; i < ap; i++) sb = sb + " ";
	return sb;
}

/**
 * Renvoie la position de la chaîne st2 dans la chaîne st1
 */
int strpos(char* st1, char* st2) {
	int i, p = 0, n1, n2;
	char* s;
	n1 = strlen(st1);
	n2 = strlen(st2);
	if (n2 < n1) {
		s = new char[n2 + 1];
		s[strlen(st2)] = '\0';
		for (p = 0; p < n1 - n2; p++) {
			for (i = 0; i < n2; i++) s[i] = st1[p + i];
			if (strcmp(s, st2) == 0) break;
		}
		delete[]s;
	}
	return p;

}

/**
 * Calcule le nombre de décimales pour l'affichage d'un nombre compris entre mini et maxi
 */
int ndecimales(double mini, double maxi) {
	double p;
	int k;
	if (mini >= 1) return 0;
	if (maxi < mini) cout << "dans ndecimales maxi<mini \n";
	p = mini / 100;
	k = 0;
	while (abs((long int)(p - round(p))) > 1E-12) {
		k++;
		p = 10 * p;
	}
	return k;
}

/**
 * transforme un nombre entier en string
 */
string IntToString(int number) {
	std::ostringstream oss;
	oss << number;
	return oss.str();
}

/**
 * transforme un short int en string
 */
string ShortIntToString(short int number) {
	std::ostringstream oss;
	oss << number;
	return oss.str();
}

/**
 *  converti l'entier number en une string de 3 caractères.
 */
string IntToString3(int number) {
	stringstream ss;//create a stringstream
	ss << number;//add number to the stream
	if ((number < 0) or (number > 999)) {
		throw std::range_error("You are trying to convert " + ss.str() +
			" into a string of length 3. Not possible.");
	}
	string ans = ss.str();
	while (ans.size() < 3) {
		ans = "0" + ans;
	}
	return ans;//return a string with the contents of the stream
}

/**
 * transforme un nombre flottant en string
 */
string FloatToString(float number) {
	std::ostringstream oss;
	oss << number;
	return oss.str();
}


/**
 * transforme un nombre <double> en string
 */
string DoubleToString(double number) {
	std::ostringstream oss;
	oss << number;
	return oss.str();
}

/**
 * transforme un nombre <long double> en string
 */
string LongDoubleToString(long double number) {
	std::ostringstream oss;
	oss << number;
	return oss.str();
}

/**
 * transforme un string en mettant en majuscules les lettres minuscules
 */
string majuscules(string s) {
	string s2;
	char* c;
	c = new char[s.length() + 1];
	for (int i = 0; i < (int)s.length(); i++) {
		c[i] = s[i];
		c[i] = toupper(c[i]);
	}
	s2 = string(c, s.length());
	delete [] c;
	return s2;
}


/**
 * transforme un string en char*
 */
/*char* str2char(string x)
{
  return x.c_str();
}
*/

/**
 * élimine les espaces et les tabulations bordant aux extrémités d'un string
 */
string trim(string s) {
	while ((s.length() > 0)and (s.substr(0, 1) == " ")) s = s.substr(1, s.length() - 1);
	while ((s.length() > 0)and (s.substr(s.length() - 1, 1) == " ")) s = s.substr(0, s.length() - 1);
	while ((s.length() > 0)and (s.substr(0, 1) == "	")) s = s.substr(1, s.length() - 1);
	while ((s.length() > 0)and (s.substr(s.length() - 1, 1) == "	")) s = s.substr(0, s.length() - 1);
	return s;
}

/**
 * arrondit un double à l'entier le plus proche
 */
int arrondi(double a) {
	return (int)(a + 0.5);
}

/**
 * Formate en string une durée exprimée en double
 */
string TimeToStr(double ti) {
	string stime = "";
	//int m_s;
	int sec, min, hou, day;
	day = int(ti) / 86400;
	ti = ti - (double)day * 86400;
	hou = int(ti) / 3600;
	ti = ti - (double)hou * 3600;
	min = int(ti) / 60;
	ti = ti - (double)min * 60;
	sec = int(ti);
	ti = ti - (double)sec;
	//m_s = int(1000*ti);
	if (day > 0) {
		stime = IntToString(day) + " day";
		if (day > 1) stime += "s ";
		else stime += " ";
	}
	if (hou > 0) stime += IntToString(hou) + " h ";
	if (min > 0) stime += IntToString(min) + " min ";
	if ((day > 0)or (hou > 0)) return stime;
	if (sec > 0) stime += IntToString(sec) + " s ";
	if (min > 0) return stime;
	//if (m_s>0) stime += IntToString(m_s)+" ms ";
	return stime;
}

/**
 * Calcule la moyenne d'un vecteur de doubles
 */
double cal_moy(int n, double* x) {
	double sx = 0;
	for (int i = 0; i < n; i++) sx += x[i];
	if (n > 0) return sx / (double)n;
	else return 0.0;
}

/**
 * Calcule la moyenne d'un vecteur de long doubles
 */
long double cal_moyL(int n, long double* x) {
	long double sx = 0;
	for (int i = 0; i < n; i++) sx += x[i];
	if (n > 0) return sx / (long double)n;
	else return 0.0;
}

/**
 * Calcule le mode d'un vecteur de doubles
 */
double cal_mode(int n, double* x) {
	int l0 = 1, l2 = n, dl, lmin;
	double min;
	dl = (l2 - l0) / 2;
	while ((x[l0 - 1] != x[l2 - 1])and (l0 != l2)) {
		min = x[l0 + dl - 1] - x[l0 - 1];
		lmin = l0;
		for (int l = l0; l < l2 - dl; l++) {
			if (x[l + dl - 1] - x[l - 1] < min) {
				min = x[l + dl - 1] - x[l - 1];
				lmin = l;
			}
		}
		l0 = lmin;
		l2 = lmin + dl;
		dl = (l2 - l0) / 2;
	}
	return x[l0 - 1];
}

/**
 * Calcule le mode d'un vecteur de long doubles
 */
long double cal_modeL(int n, long double* x) {
	int l0 = 1, l2 = n, dl, lmin;
	long double min;
	dl = (l2 - l0) / 2;
	while ((x[l0 - 1] != x[l2 - 1])and (l0 != l2)) {
		min = x[l0 + dl - 1] - x[l0 - 1];
		lmin = l0;
		for (int l = l0; l < l2 - dl; l++) {
			if (x[l + dl - 1] - x[l - 1] < min) {
				min = x[l + dl - 1] - x[l - 1];
				lmin = l;
			}
		}
		l0 = lmin;
		l2 = lmin + dl;
		dl = (l2 - l0) / 2;
		//cout<<"l0="<<l0<<"    l2="<<l2<<"\n";
	}
	return x[l0 - 1];
}

/**
 * Calcule la médiane d'un vecteur de doubles
 */
double cal_med(int n, double* x) {
	sort(&x[0], &x[n]);
	if ((n % 2) == 0) return 0.5 * (x[n / 2 - 1] + x[n / 2]);
	else return x[n / 2];
}

/**
 * Calcule la médiane d'un vecteur de long doubles
 */
long double cal_medL(int n, long double* x) {
	sort(&x[0], &x[n]);
	if ((n % 2) == 0) return 0.5 * (x[n / 2 - 1] + x[n / 2]);
	else return x[n / 2];
}

/**
 * Calcule l'écart-type d'un vecteur de doubles
 */
double cal_sd(int n, double* x) {
	double sx, sx2, a = x[0];
	bool ident = true;
	int i = 1;
	while ((ident) and (i < n)) {
		ident = (x[i] == a);
		i++;
	}
	if (ident) return 0.0;
	sx = 0.0;
	sx2 = 0.0;
	for (int i = 0; i < n; i++) {
		sx += x[i];
		sx2 += x[i] * x[i];
	}
	if (n < 2) return 0.0;
	else return sqrt((sx2 - (sx / (double)n * sx)) / (double)(n - 1));
}

/**
 * Calcule l'écart-type d'un vecteur de long doubles
 */
long double cal_sdL(int n, long double* x) {
	long double sx, sx2, a = x[0];
	bool ident = true;
	int i = 1;
	while ((ident) and (i < n)) {
		ident = (x[i] == a);
		i++;
	}
	if (ident) return 0.0;
	sx = 0.0;
	sx2 = 0.0;
	for (int i = 0; i < n; i++) {
		sx += x[i];
		sx2 += x[i] * x[i];
	}
	if (n < 2) return 0.0;
	else return sqrt((sx2 - (sx / (long double)n * sx)) / (long double)(n - 1));
}

/**
 * Calcule la variance d'un vecteur de long doubles
 */
long double cal_varL(int n, long double* x) {
	long double sx, sx2, a = x[0];
	bool ident = true;
	int i = 1;
	while ((ident) and (i < n)) {
		ident = (x[i] == a);
		i++;
	}
	if (ident) return 0.0;
	sx = 0.0;
	sx2 = 0.0;
	for (int i = 0; i < n; i++) {
		sx += x[i];
		sx2 += x[i] * x[i];
	}
	if (n < 2) return 0.0;
	else return (sx2 - (sx / (long double)n * sx)) / (long double)(n - 1);
}

/**
 * retourne le logarithme de la densité d'une loi normale de moyenne m et d'écart-type s
 */
double lnormal_dens(double x, double m, double s) {
	return -0.5 * ((x - m) / s) * ((x - m) / s) - log(s) - 0.5 * log(2 * PI);
}

/**
 * fonction utilisée par la fonction pnorm5
 */
double do_del(double temp, double x) { //calcule cum avec ccum=1.0-cum
	double xsq, del;
	xsq = floor(x * 16.0) / 16.0;
	del = (x - xsq) * (x + xsq);
	return exp(-0.5 * xsq * xsq) * exp(-0.5 * del) * temp;
}

/**
 *
 */
double pnorm5(double x, double mu, double sigma) {
	double a[5] = {2.2352520354606839287,161.02823106855587881,1067.6894854603709582,18154.981253343561249,0.065682337918207449113};
	double b[4] = {47.20258190468824187,976.09855173777669322,10260.932208618978205,45507.789335026729956};
	double c[9] = {0.39894151208813466764,8.8831497943883759412,93.506656132177855979,597.27027639480026226,2494.5375852903726711,6848.1904505362823326,11602.651437647350124,9842.7148383839780218,1.0765576773720192317e-8};
	double d[8] = {22.266688044328115691,235.38790178262499861,1519.377599407554805,6485.558298266760755,18615.571640885098091,34900.952721145977266,38912.003286093271411,19685.429676859990727};
	double p[6] = {0.21589853405795699,0.1274011611602473639,0.022235277870649807,0.001421619193227893466,2.9112874951168792e-5,0.02307344176494017303};
	double q[5] = {1.28426009614491121,0.468238212480865118,0.0659881378689285515,0.00378239633202758244,7.29751555083966205e-5};
	//    double std::numeric_limits<double>::epsilon()=2.2204460492503131E-16;
	double sqrt2PI = 2.506628274631000502415765284811;
	double pp, xden, xnum, temp, eps, xsq, y, cum;

	if (sigma <= 0.0) {
		if (x < mu) return 0.0;
		else return 1.0;
	}
	pp = (x - mu) / sigma;
	if (pp < 0) y = -pp;
	else y = pp;
	eps = 0.5 * std::numeric_limits<double>::epsilon();
	if (y <= 0.67448975) {
		if (y > eps) {
			xsq = pp * pp;
			xnum = a[4] * xsq;
			xden = xsq;
			for (int i = 0; i < 3; i++) {
				xnum = (xnum + a[i]) * xsq;
				xden = (xden + b[i]) * xsq;
			}
		}
		else {
			xnum = 0.0;
			xden = 0.0;
		}
		temp = pp * (xnum + a[3]) / (xden + b[3]);
		cum = 0.5 + temp;
	}
	else {
		if (y <= sqrt(32.0)) {
			xnum = c[8] * y;
			xden = y;
			for (int i = 0; i < 7; i++) {
				xnum = (xnum + c[i]) * y;
				xden = (xden + d[i]) * y;
			}
			temp = (xnum + c[7]) / (xden + d[7]);
			cum = do_del(temp, y);
			if (pp > 0.0) cum = 1.0 - cum;
		}
		else {
			xsq = 1.0 / (pp * pp);
			xnum = p[5] * xsq;
			xden = xsq;
			for (int i = 0; i < 4; i++) {
				xnum = (xnum + p[i]) * xsq;
				xden = (xden + q[i]) * xsq;
			}
			temp = xsq * (xnum + p[4]) / (xden + q[4]);
			temp = (1.0 / sqrt2PI - temp) / y;
			cum = do_del(temp, pp);
			if (pp > 0.0) cum = 1.0 - cum;
		}
	}
	return cum;
}


rescov covarianceponderee(int nl, int nc, long double** A, long double* w) {
	rescov mcv;
	mcv.n = nc;
	mcv.mu = new long double[nc];
	mcv.cov = new long double*[nc];
	for (int i = 0; i < nc; i++) mcv.cov[i] = new long double[nc];
	long double *wn, sw2, sw;
	//sx = new double[nc];
	wn = new long double[nl];
	sw = 0.0;
	for (int k = 0; k < nl; k++) sw += w[k];
	for (int k = 0; k < nl; k++) wn[k] = w[k] / sw;
	sw2 = 0.0;
	for (int k = 0; k < nl; k++) sw2 += wn[k] * wn[k];
	for (int i = 0; i < nc; i++) {
		mcv.mu[i] = 0.0;
		for (int k = 0; k < nl; k++) mcv.mu[i] += wn[k] * A[k][i];
		//mu[i]=sx[i];
		for (int j = 0; j <= i; j++) {
			mcv.cov[i][j] = 0.0;
			for (int k = 0; k < nl; k++) mcv.cov[i][j] += wn[k] * A[k][i] * A[k][j];
			mcv.cov[j][i] = mcv.cov[i][j];
		}
	}
	for (int i = 0; i < nc; i++) {
		for (int j = 0; j < nc; j++) mcv.cov[i][j] = (mcv.cov[i][j] - mcv.mu[i] * mcv.mu[j]) / (1.0 - sw2);
	}
	delete []wn;
	return mcv;
}


resAFD AFD(int nl, int nc, int* pop, long double* omega, long double** X, long double prop) {
	resAFD res;
	rescov mcv;
	res.proportion = prop;
	vector<int> nk, numpop; //npop est remplace par numpop.size() ou nk.size()
	bool trouve;
	long double **matC, **matCT, **matTI, **matM, *valprop, **mk, piv, sl, co, somega, *w, *wk;
	int snk, j, err;
	//cout<<"debut AFD  nl="<<nl<<"\n";
	if (not numpop.empty()) numpop.clear();
	numpop.push_back(pop[0]);
	if (not nk.empty()) nk.clear();
	nk.push_back(1);
	for (int i = 1; i < nl; i++) {
		trouve = false;
		for (j = 0; j < (int)numpop.size(); j++) {
			trouve = (numpop[j] == pop[i]);
			if (trouve) break;
		}
		if (trouve) {
			nk[j]++;
		}
		else {
			numpop.push_back(pop[i]);
			nk.push_back(1);
		}
	}
	snk = 0;
	for (size_t i = 0; i < nk.size(); i++) snk += nk[i];
	//cout<<"apres calcul snk\n";
	if (nl != snk) {
		cout << "dans AFD nl(" << nl << ") !=snk(" << snk << ")\n";
		exit(1);
	}
	matC = new long double*[nc];
	for (int i = 0; i < nc; i++) matC[i] = new long double[numpop.size()];
	matCT = new long double*[numpop.size()];
	for (size_t i = 0; i < numpop.size(); i++) matCT[i] = new long double[nc];
	//matM = new *double[nc];for (int i=0;i<nc;i++) matM[i] = new double[nc];
	valprop = new long double[nc];
	res.vectprop = new long double*[nc];
	for (int i = 0; i < nc; i++) res.vectprop[i] = new long double[nc];
	//anl=1.0/(double)nl;
	somega = 0.0;
	for (int i = 0; i < nl; i++) somega += omega[i];
	//cout<<"somega="<<somega<<"\n";
	w = new long double[nl];
	for (int i = 0; i < nl; i++) w[i] = omega[i] / somega;
	//cout<<"avant covarianceponderee\n";
	mcv = covarianceponderee(nl, nc, X, w);
	//cout<<"apres covarianceponderee\n";
	res.moy = new long double[mcv.n];
	for (j = 0; j < mcv.n; j++) res.moy[j] = mcv.mu[j];
	mk = new long double*[numpop.size()];
	for (size_t i = 0; i < numpop.size(); i++) mk[i] = new long double[nc];
	wk = new long double[numpop.size()];
	//cout<<"numpop.size="<<numpop.size()<<"\n";
	for (size_t k = 0; k < numpop.size(); k++) {
		wk[k] = 0.0;
		for (int i = 0; i < nl; i++) {
			if (pop[i] == numpop[k]) {
				wk[k] += w[i];
			}
		}
		//cout<<"wk["<<k<<"]="<<wk[k]<<"\n";
		for (j = 0; j < nc; j++) {
			mk[k][j] = 0.0;
			for (int i = 0; i < nl; i++) if (pop[i] == numpop[k]) mk[k][j] += w[i] * X[i][j] / wk[k];
		}
	}
	//cout<<"apres la premiere boucle\n";
	for (size_t k = 0; k < numpop.size(); k++) {
		//cout<<"k="<<k<<"\n";
		//cout<<"   w[k]="<<wk[k]<<"\n";
		co = sqrt(wk[k]);
		for (j = 0; j < nc; j++) {
			matC[j][k] = co * (mk[k][j] - res.moy[j]);
			matCT[k][j] = matC[j][k];
			//cout<<"k="<<k<<"   j="<<j<<"   co="<<co<<"   mk[k][j]="<<mk[k][j]<<"   res.moy[j]="<<res.moy[j]<<"\n";
		}
	}
	matTI = new long double*[mcv.n];
	for (int i = 0; i < mcv.n; i++) matTI[i] = new long double[mcv.n];
	//cout<<"avant invM\n";
	err = inverse_Tik(mcv.n, mcv.cov, matTI);
	if (err > 0) cout << "err=" << err << " dans inverse_Tik\n";
	matM = prodML(numpop.size(), nc, numpop.size(), prodML(numpop.size(), nc, nc, matCT, matTI), matC);
	//cout<<"avant jacobi\n";
	jacobiL(numpop.size(), matM, valprop, res.vectprop);
	//cout<<"apres jaconi\n";
	res.vectprop = prodML(nc, numpop.size(), numpop.size(), prodML(nc, nc, numpop.size(), matTI, matC), res.vectprop);
	//tri des valeurs et vecteurs propres par ordre décroissant des valeurs propres
	for (int i = 0; i < (int)numpop.size() - 1; i++) {
		for (j = i + 1; j < (int)numpop.size(); j++)
			if (valprop[i] < valprop[j]) {
				piv = valprop[i];
				valprop[i] = valprop[j];
				valprop[j] = piv;
				for (int k = 0; k < nc; k++) {
					piv = res.vectprop[k][i];
					res.vectprop[k][i] = res.vectprop[k][j];
					res.vectprop[k][j] = piv;
				}
			}
	}
	res.slambda = 0.0;
	for (size_t i = 0; i < numpop.size(); i++) res.slambda += valprop[i];
	res.nlambda = 1;
	sl = valprop[0];
	while ((res.nlambda < (int)numpop.size())and (sl / res.slambda < prop)) {
		res.slambda += valprop[res.nlambda];
		res.nlambda++;
	}
	res.lambda = new long double[res.nlambda];
	for (int i = 0; i < res.nlambda; i++) res.lambda[i] = valprop[i];
	res.princomp = new long double*[nl];
	for (int i = 0; i < nl; i++) res.princomp[i] = new long double[res.nlambda];
	for (int i = 0; i < nl; i++) {
		for (j = 0; j < res.nlambda; j++) {
			res.princomp[i][j] = 0.0;
			for (int k = 0; k < nc; k++) res.princomp[i][j] += (X[i][k] - res.moy[k]) * res.vectprop[k][j];
		}
	}
	for (j = 0; j < mcv.n; j++) delete []mcv.cov[j];
	delete []mcv.cov;
	delete []mcv.mu;
	for (int i = 0; i < nc; i++) delete []matC[i];
	delete []matC;
	for (size_t i = 0; i < numpop.size(); i++) delete []matCT[i];
	delete []matCT;
	for (int i = 0; i < nc; i++) delete []matTI[i];
	delete []matTI;
	for (size_t i = 0; i < numpop.size(); i++) delete []matM[i];
	delete []matM;
	for (size_t i = 0; i < numpop.size(); i++) delete []mk[i];
	delete []mk;
	delete []valprop;
	delete []w;
	delete []wk;
	return res;
}

int PGCD(int a, int b) {
	while (1) {
		a = a % b;
		if (a == 0) return b;
		b = b % a;
		if (b == 0) return a;
	}
}

int PPCM(int a, int b) {
	return a * b / PGCD(a, b);
}

bool compC::operator()(const matC& lhs, const matC& rhs) const {
	return lhs.v < rhs.v;
}

long double DCVM(int n, int m, long double* x, long double* y) {
	matC* A;
	A = new matC[n + m];
	for (int i = 0; i < n; i++) {
		A[i].v = x[i];
		A[i].num = 0;
		A[i].ind = i;
	}
	for (int i = 0; i < m; i++) {
		A[n + i].v = y[i];
		A[n + i].num = 1;
		A[n + i].ind = i;
	}
	sort(&A[0], &A[n + m], compC());
	long double L, a, b, h = 0.0, sh = 0.0, d;
	L = (long double)PPCM(n, m);
	a = L / (long double)n;
	b = L / (long double)m;
	for (int i = 0; i < m + n; i++) {
		if (A[i].num == 0) h += 1;
		else h -= 1;
		sh += sqr(h);
	}
	a = (long double)n;
	b = (long double)m;
	d = (a * b) / sqr(a + b) / sqr(L) * sh;
	return d;
}

void combrank2(int n, int m, long double* x, long double* y, long double* rangx, long double* rangy) {
	matC* A;
	int k, i0, i1;
	//	long double rmin,rmax;
	A = new matC[n + m];
	for (int i = 0; i < n; i++) {
		A[i].v = x[i];
		A[i].num = 0;
		A[i].ind = i;
	}
	for (int i = 0; i < m; i++) {
		A[n + i].v = y[i];
		A[n + i].num = 1;
		A[n + i].ind = i;
	}
	sort(&A[0], &A[n + m], compC());
	i0 = 0;
	i1 = 0;
	for (int i = 0; i < n; i++) rangx[i] = 0.0;
	for (int i = 0; i < m; i++) rangy[i] = 0.0;
	while (i1 < n + m) {
		while ((i1 < n + m - 1)and (A[i1 + 1].v == A[i0].v)) i1++;
		for (int i = i0; i <= i1; i++) {
			k = A[i].ind;
			if (A[i].num == 0) rangx[k] = 1.0 + 0.5 * (long double)(i0 + i1);
			else rangy[k] = 1.0 + 0.5 * (long double)(i0 + i1);
		}
		i0 = i1 + 1;
		i1 = i0;
	}
	delete []A;
}


void ecritTable(int** table, int m, int n) {
	cout << endl << endl;
	for (int i = 0; i < m; ++i) {
		for (int j = 0; j < n; ++j) {
			cout << table[i][j] << " ";
		}
		cout << endl;
	}
}
