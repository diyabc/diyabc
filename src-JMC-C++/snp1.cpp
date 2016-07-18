
#include <fstream>
#include <iostream>
#include <string.h>
#include <string>
#include <sstream>
#include <unistd.h>
#include <sys/time.h>
#include <sys/stat.h>
#include <cstdlib>
#include <stdexcept>
#include <vector>

#include "mesutils.h"


using namespace std;

/**
 * renvoit la valeur entière contenue dans la num-ième sous-chaîne de s
 * fonctionne également si la valeur entière est encadrée par des parenthèses ou des crochets
 */
int getwordint(string s,int num){
	s.append(" ");
	while (s.find(" ")== 0) s=s.substr(1);
	int i=1,j;
	while ((i<num)and(s.length()>1)) {
		j=s.find(" ");
		s=s.substr(j);
		while (s.find(" ")== 0) s=s.substr(1);
		i++;
	}
	s=s.substr(0,s.find(" "));
	if ((s.find("(")==0)or(s.find("[")==0)) s=s.substr(1,s.length()-2);
	return atoi(s.c_str());
}

/**
 * renvoit la valeur flottante contenue dans la num-ième sous-chaîne de s
 * fonctionne également si la valeur flottante est encadrée par des parenthèses ou des crochets
 */
double getwordfloat(string s,int num){
	s.append(" ");
	while (s.find(" ")== 0) s=s.substr(1);
	int i=1,j;
	while ((i<num)and(s.length()>1)) {
		j=s.find(" ");
		s=s.substr(j);
		while (s.find(" ")== 0) s=s.substr(1);
		i++;
	}
	s=s.substr(0,s.find(" "));
	if ((s.find("(")==0)or(s.find("[")==0)) s=s.substr(1,s.length()-2);
	return atof(s.c_str());
}

/**
 * découpe la chaîne s en sous-chaînes séparées par le séparateur sep
 * le nombre de sous chaînes est donné par *k
 */
string* splitwords(string s,string sep,int *k){
	int j=0,j0;
	while (s.find(sep)== 0) s=s.substr(1);
	*k=0;
	s.append(sep);
	//cout<<"s="<<s<<"\n";
	string *sb,s0,s1;
	sb=NULL;
	s1=string();
	for (int i = 0; i < (int)s.length(); i++) {
		s0=s.substr(i,1);
		if (s0==sep){
			j++;
			if (j==1) {
				s1.append(s0);
				if (j==1) (*k)++;
				//cout <<" j=1  k="<<*k<<"\n";
			}
		} else {
			s1.append(s0);
			j=0;
		}
	}
	sb = new string[*k];
	for (int i=0;i<*k;i++) {
		j0=s1.find(sep);
		sb[i]=s1.substr(0,j0);
		s1=s1.substr(j0+1,s.length());

	}
	//cout <<"k="<<*k<<"\n";
	return sb;
}
void splitwords(string s, string sep, vector<string>& resultat) {
	int j=0,j0;
	while (s.find(sep) == 0) s=s.substr(1);
	int k=0;
	s.append(sep);
	//cout<<"s="<<s<<"\n";
	vector<string> sb(0);
	string s0,s1;
	s1=string();
	for (int i = 0; i < (int)s.length(); i++) {
		s0=s.substr(i,1);
		if (s0==sep){
			j++;
			if (j==1) {
				s1.append(s0);
				if (j==1) (k)++;
				//cout <<" j=1  k="<<*k<<"\n";
			}
		} else {
			s1.append(s0);
			j=0;
		}
	}
	sb.resize(k);
	for (int i=0;i<k;i++) {
		j0=s1.find(sep);
		sb[i]=s1.substr(0,j0);
		s1=s1.substr(j0+1,s.length());
	}
	//cout <<"k="<<*k<<"\n";
	resultat.swap(sb);
}
/**
 * change les tabulations en espaces
 */

string purgetab(string s) {
	string ss="",s0;
	for (int i=0; i < (int)s.length(); i++) {
		s0=s.substr(i,1);
		if (s0=="\t") s0=" ";
		ss +=s0;
	}
//	if (ss.at(ss.length()-1)=='\r')cout<<"HELLO\n";
	return ss;
}

/**
 * découpe la chaîne s en sous-chaînes séparées par le séparateur sep
 * le découpage s'arrête quand le nombre de sous-chaînes atteint m
 */
string* splitwordsR(string s,string sep,int m,int *k){
	int j=0,j0;
	while (s.find(sep)== 0) s=s.substr(1);
	*k=0;
	string *sb,s0,s1;
	if (s.length()==0) {sb = new string[1];sb[0]="";return sb;}
	s.append(sep);
	s1=string();
	for (int i=0; i < (int)s.length(); i++) {
		s0=s.substr(i,1);
		if (s0==sep){
			j++;
			if (j==1) {
				s1.append(s0);
				if (j==1) (*k)++;
				//cout <<" j=1  k="<<*k<<"\n";
			}
		} else {
			s1.append(s0);
			j=0;
		}
		if ((*k)==m) break;
	}
	sb = new string[*k];
	for (int i=0;i<*k;i++) {
		j0=s1.find(sep);
		sb[i]=s1.substr(0,j0);
		s1=s1.substr(j0+1,s.length());

	}
	//cout <<"k="<<*k<<"\n";
	return sb;
}

void splitwordsR(string s,string sep,int m, vector<string>& resultat){
	int j=0,j0;
	while (s.find(sep)== 0) s=s.substr(1);
	int k=0;
	string s0,s1;
	if (s.length()==0) {resultat.resize(0); return;}
	s.append(sep);
	s1=string();
	for (int i=0; i < (int)s.length(); i++) {
		s0=s.substr(i,1);
		if (s0==sep){
			j++;
			if (j==1) {
				s1.append(s0);
				if (j==1) k++;
				//cout <<" j=1  k="<<*k<<"\n";
			}
		} else {
			s1.append(s0);
			j=0;
		}
		if (k==m) break;
	}
	vector<string> sb(k);
	for (int i=0;i<k;i++) {
		j0=s1.find(sep);
		sb[i]=s1.substr(0,j0);
		s1=s1.substr(j0+1,s.length());

	}
	resultat.swap(sb);
	//cout <<"k="<<*k<<"\n";
}

int lifichier(string filename,string outfile){
    ifstream file0(filename.c_str(), ios::in);
    if (not file0) {cout<<"file0=NULL\n";return -1;}
	string ligne;
	getline(file0,ligne);
	//splitwords(ligne," ", ss); nss = ss.size();
	int ngenes,nloc,n0,nind,k1,k2,loc;
	ngenes=getwordint(ligne,2);nind=ngenes/2;
	nloc=getwordint(ligne,3);
	vector <int> hist;
	hist.resize(ngenes+1);
	for (int i=0;i<ngenes+1;i++) hist[i]=0;
	vector<vector<int> > locus;
    locus = vector<vector<int> >(nind);
    for (int i=0;i<nind;i++) {
        locus[i] = vector<int>(nloc);
    }	
	cout<<"samplesize="<<ngenes<<"   nloc="<<nloc<<"\n";
	getline(file0,ligne);
	for (int loc=0;loc<nloc;loc++) {
		for (int i=0;i<5;i++) getline(file0,ligne);
		n0=0;
		for (int i=0;i<nind;i++) {
			getline(file0,ligne);
			k1=getwordint(ligne,1);
			if (k1==0) n0++;
			getline(file0,ligne);
			k2=getwordint(ligne,1);
			if (k2==0) n0++;
			locus[i][loc]=k1+k2;
		}
		hist[n0]++;
	}
	file0.close();
	for (int i=0;i<ngenes+1;i++) cout<<i<<"      ";cout<<"\n";
	for (int i=0;i<ngenes+1;i++) cout<<hist[i]<<"    ";cout<<"\n";
	ofstream file1(outfile.c_str(), ios::out);
	file1<<"<NM=1.0NF>\n";
	file1<<"IND   SEX   POP   ";
	string sind="P1_";
	for (loc=0;loc<nloc;loc++) file1<<" A";file1<<"\n";
	for (int i=0;i<nind;i++) {
		file1<<sind<<i<<"   F     P1    ";
		for (loc=0;loc<nloc;loc++) file1<<" "<<locus[i][loc];
		file1<<"\n";
	}
	file1.close();
}

int main(int argc, char *argv[]){
	if (argc<3) {cout<<"Usage snp1 <infile> <outfile>"<<"\n"<<argc<<"n";exit(1);}
	string nomfiin,nomfiout;
	int err;
	nomfiin=argv[1];
	cout<<"fichier en entrée: "<<nomfiin<<"\n";
	nomfiout=argv[2];
	cout<<"fichier en sortie: "<<nomfiout<<"\n";
	err=lifichier(nomfiin,nomfiout);
};