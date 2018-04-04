/*
 * sumstat.cpp
 *
 *  Created on: 8 déc. 2011
 *      Author: ppudlo
 */

#include <string>
#include <vector>
#include <iostream>
#include <utility>
#include <algorithm>
#include <numeric>
#include <cmath>
#include <cstdlib>
#include <math.h>

#include "randomgenerator.hpp"
#include "mesutils.hpp"
#include "data.hpp"
#include "history.hpp"
#include "particuleC.hpp"
#include "statdefs.hpp"

extern int debuglevel;
extern DataC dataobs;
extern vector<statn> stats;
/***********************************************************************************************************************/

int ParticleC::samplesize(int loc, int sample) {
	int cat = this->locuslist[loc].type % 5;
	int n = dataobs.ssize[cat][sample];
	if (this->locuslist[loc].type < 10) {
		//cout<<"samplesize["<<loc<<"]["<<sample<<"]   nmisshap="<<dataobs.nmisshap<<"\n";
		if (dataobs.nmisshap > 0) {
			for (int i = 0; i < dataobs.nmisshap; i++) {
				if ((dataobs.misshap[i].locus == loc)and (dataobs.misshap[i].sample == sample)) n--;
			}
		}
	}
	else {
		n = this->locuslist[loc].samplesize[sample];
		/*
	  if (dataobs.nmisssnp>0) {
		  for (int i=0;i<dataobs.nmisssnp;i++){
			  if ((dataobs.misssnp[i].locus == loc)and(dataobs.misssnp[i].sample == sample)) n--;
		  }
	  }*/
	}
	return n;
}

void ParticleC::calfreq(int gr) {
	//int n=0;
	//cout <<"debut de calfreq nloc = "<<this->nloc<<"  groupe = "<<gr <<"   this->nsample="<<this->nsample<<"\n";
	int loc, iloc, cat, sasize;
	for (iloc = 0; iloc < this->grouplist[gr].nloc; iloc++) {
		loc = this->grouplist[gr].loc[iloc];
		cat = this->locuslist[loc].type % 5;
		//cout <<"\nLocus "<< loc <<"   kmin="<<this->locuslist[loc].kmin<<"   kmax="<<this->locuslist[loc].kmax<<"\n";
		this->locuslist[loc].nal = this->locuslist[loc].kmax - this->locuslist[loc].kmin + 1;
		this->locuslist[loc].freq = vector<vector<long double>>(this->nsample);
		//		this->locuslist[loc].freq.resize(this->nsample);
		//cout <<"nal ="<<this->locuslist[loc].nal<<"\n";
		for (int samp = 0; samp < this->nsample; samp++) {
			this->locuslist[loc].freq[samp] = vector<long double>(this->locuslist[loc].nal);
			for (int i = 0; i < this->locuslist[loc].nal; i++) this->locuslist[loc].freq[samp][i] = 0.0;
			//cout << dataobs.ssize[cat][samp] <<"\n";
			for (int i = 0; i < dataobs.ssize[cat][samp]; i++) {
				//cout<<"i="<<i<<"\n";
				//cout<<"this->locuslist[loc].haplomic[samp][i]="<<this->locuslist[loc].haplomic[samp][i]<<"\n";
				//exit(1);
				if (this->locuslist[loc].haplomic[samp][i] != MICMISSING) {
					//cout <<"  "<<this->locuslist[loc].haplomic[samp][i]<<"\n";
					this->locuslist[loc].freq[samp][this->locuslist[loc].haplomic[samp][i] - this->locuslist[loc].kmin] += 1.0;
					//n++;
				}
			}
			//cout <<" fini\n";
			for (int k = 0; k < this->locuslist[loc].kmax - this->locuslist[loc].kmin + 1; k++) {
				sasize = samplesize(loc, samp);
				if (sasize > 0) this->locuslist[loc].freq[samp][k] /= sasize;

			}
			//cout << "loc "<<loc <<"  sa "<<samp;
			//for (int k=0;k<this->locuslist[loc].nal;k++) {if (this->locuslist[loc].freq[samp][k]>0.0) cout << "   "<<this->locuslist[loc].freq[samp][k]<<"("<<k+this->locuslist[loc].kmin<<")";}
			//cout << "\n";
		}
	}
	//cout <<"fin de calfreq \n";//exit(1);
}

void ParticleC::calfreqsnp(int gr) {
	int loc, iloc, sasize;
	//int cat;
	short int g0 = 0;
	//long double h0,h1,h2,m0=0.0,m1=0.0,m2=0.0;
	//cout << "calfreqsnp début   nloc="<<this->grouplist[gr].nloc<<"\n";
	for (iloc = 0; iloc < this->grouplist[gr].nloc; iloc++) {
		loc = this->grouplist[gr].loc[iloc];
		//cout<<"locus "<<loc<<"\n";
		if (this->locuslist[loc].weight > 0.0) {
			//cout<<"this->locuslist[loc].type="<<this->locuslist[loc].type<<"\n";
			//cat=this->locuslist[loc].type%5;
			//cout<<"cat = "<<cat<<"   nsample="<<this->nsample<<"\n";
			this->locuslist[loc].freq = vector<vector<long double>>(this->nsample);
			for (int samp = 0; samp < this->nsample; samp++) {
				//cout<<"sample "<<samp<<"\n";
				sasize = this->locuslist[loc].samplesize[samp];
				this->locuslist[loc].freq[samp] = vector<long double>(2);
				this->locuslist[loc].freq[samp][0] = 0.0;
				this->locuslist[loc].freq[samp][1] = 0.0;
				//cout<<"avant la boucle i\n";
				//cout<<"dataobs.ssize[cat][samp]="<<dataobs.ssize[cat][samp]<<"\n";
				for (int i = 0; i < sasize; i++) {
					//if ((this->locuslist[loc].haplosnp[samp][i] == g0)/*and(this->dat[cat][samp][i])*/) this->locuslist[loc].freq[samp][0] +=1.0;
					if (this->locuslist[loc].haplosnp[samp][i] == g0) this->locuslist[loc].freq[samp][0] += 1.0;
				}
				//cout<<"apres la boucle i\n";
				this->locuslist[loc].freq[samp][0] /= sasize;
				this->locuslist[loc].freq[samp][1] = 1.0 - this->locuslist[loc].freq[samp][0];
				//if (loc<10) printf("loc %3d    sa %d     samplesize %3d   freq[0]=%10.6Lf\n",loc,samp,this->locuslist[loc].samplesize[samp],this->locuslist[loc].freq[samp][0]);
			}
		}
	}
}

void ParticleC::calfreqpool(int gr) {
	int loc, iloc, cat;
	double x1, x;
	for (iloc = 0; iloc < this->grouplist[gr].nloc; iloc++) {
		loc = this->grouplist[gr].loc[iloc];
		//cout<<"locus "<<loc<<"\n";
		if (this->locuslist[loc].weight > 0.0) {
			//cout<<"this->locuslist[loc].type="<<this->locuslist[loc].type<<"\n";
			cat = this->locuslist[loc].type % 5;
			//cout<<"cat = "<<cat<<"   nsample="<<this->nsample<<"\n";
			this->locuslist[loc].freq = vector<vector<long double>>(this->nsample);
			for (int samp = 0; samp < this->nsample; samp++) {
				this->locuslist[loc].freq[samp] = vector<long double>(2);
				x1 = (double)this->locuslist[loc].nreads1[samp];
				x = (double)this->locuslist[loc].nreads[samp];
				this->locuslist[loc].freq[samp][0] = x1 / x;
				this->locuslist[loc].freq[samp][1] = 1.0 - this->locuslist[loc].freq[samp][0];
				/*/-------------------------------------------------A SUPPRIMER
									short int g0=0;
									long double freq0=0.0;
									int sasize=this->locuslist[loc].samplesize[samp];
									for (int i=0;i<sasize;i++){
										if (this->locuslist[loc].haplosnp[samp][i] == g0)  freq0 +=1.0;
									}
									freq0 /=sasize;
									double perr;
									if (this->locuslist[loc].freq[samp][0]>0.0) perr = (this->locuslist[loc].freq[samp][0]-freq0)/this->locuslist[loc].freq[samp][0];
									else perr = (this->locuslist[loc].freq[samp][0]-freq0); 
				//-------------------------------------------------*/
				if ((debuglevel == 18)and (loc < 5)) {
					//printf("loc %3d    sa %d     samplesize %3d   freq[0]=%10.6Lf (%10.6Lf)  perr=%10f  %d\n",loc,samp,this->locuslist[loc].samplesize[samp],this->locuslist[loc].freq[samp][0],freq0,perr,(int)x1);
					printf("loc %3d    sa %d     samplesize %3d   freq[0]=%10.6Lf \n", loc, samp, this->locuslist[loc].samplesize[samp], this->locuslist[loc].freq[samp][0]);
				}
			}
		}
	}
}

void ParticleC::liberefreq(int gr) {
	int loc, iloc;
	for (iloc = 0; iloc < this->grouplist[gr].nloc; iloc++) {
		loc = this->grouplist[gr].loc[iloc];
		if (this->locuslist[loc].groupe > 0) {
			if (this->locuslist[loc].type < 5) {
				this->locuslist[loc].freq.clear();
				this->locuslist[loc].haplomic.clear();
			}
			else if ((this->locuslist[loc].type > 9)and (this->locuslist[loc].weight > 0.0)) {
				this->locuslist[loc].freq.clear();
				this->locuslist[loc].haplosnp.clear();
			}
		}
	}
}

void ParticleC::liberefreqsnp(int gr) {
	int loc, iloc;
	for (iloc = 0; iloc < this->grouplist[gr].nloc; iloc++) {
		loc = this->grouplist[gr].loc[iloc];
		if (this->locuslist[loc].groupe > 0) {
			if (this->locuslist[loc].type < 5) {
				this->locuslist[loc].freq.clear();
				this->locuslist[loc].haplomic.clear();
			}
		}
	}
}

void ParticleC::libere_freq(int gr) {
	int kloc, iloc;
	for (iloc = 0; iloc < this->grouplist[gr].nloc; iloc++) {
		kloc = this->grouplist[gr].loc[iloc];
		//cout<<"libere_freq 1  this-<nsample="<<this->nsample<<"\n";
		this->locuslist[kloc].freq.clear();
		//cout<<"libere_freq 2\n";
		this->locuslist[kloc].haplomic.clear();
	}
}

void ParticleC::liberednavar(int gr) {
	int kloc, iloc;
	for (iloc = 0; iloc < this->grouplist[gr].nloc; iloc++) {
		kloc = this->grouplist[gr].loc[iloc];
		this->locuslist[kloc].haplodnavar.clear();
		this->locuslist[kloc].haplodna.clear();
	}
}

/*
  void ParticleC::cal_snhet(int gr,int numsnp) {
  long double *het;
  int iloc,kloc,nl=0;
  int sample=this->grouplist[gr].sumstatsnp[numsnp].samp-1;
  het = new long double[this->grouplist[gr].nloc];
  for (iloc=0;iloc<this->grouplist[gr].nloc;iloc++){
  kloc=this->grouplist[gr].loc[iloc];
  if (this->locuslist[kloc].samplesize[sample]>1){
  het[nl]=1.0;
  for (int k=0;k<2;k++) het[nl] -= sqr(this->locuslist[kloc].freq[sample][k]);
  het[nl] *= ((long double)this->locuslist[kloc].samplesize[sample]/((long double)this->locuslist[kloc].samplesize[sample]-1));
  nl++;
  }
  }
  this->grouplist[gr].sumstatsnp[numsnp].n = nl;
  this->grouplist[gr].sumstatsnp[numsnp].x =new long double[nl];
  for (int loc=0;loc<nl;loc++) this->grouplist[gr].sumstatsnp[numsnp].x[loc] = het[loc];
  this->grouplist[gr].sumstatsnp[numsnp].defined=true;
  delete []het;
  }
*/

////////////// Stats SNP
void ParticleC::cal_snhet(int gr, int numsnp) {
	//cout<<"debut de cal_snhet   numsnp="<<numsnp<<"\n";
	long double het;
	int iloc, kloc, sasize;
	long double x = 0, w;
	auto&& samp = this->grouplist[gr].sumstatsnp[numsnp].samp.get();
	int sample = samp[0];
	//ofstream f10("/home/cornuet/workspace/diyabc/src-JMC-C++/hetero.txt",ios::out);
	this->grouplist[gr].sumstatsnp[numsnp].n = this->grouplist[gr].nloc;
	for (iloc = 0; iloc < this->grouplist[gr].nloc; iloc++) {
		kloc = this->grouplist[gr].loc[iloc];
		//cout<<"kloc="<<kloc<<"    weight="<<this->locuslist[kloc].weight<<"\n";
		w = this->locuslist[kloc].weight;
		if (w > 0.0) {
			//cout<<"avant sasize\n";
			sasize = samplesize(kloc, sample);
			//cout<<"sasize="<<sasize<<"\n";
			if (sasize > 1) {
				het = 1.0;
				for (int k = 0; k < 2; k++) het -= sqr(this->locuslist[kloc].freq[sample][k]);
				het *= ((long double)sasize / (long double)(sasize - 1));
				x = het;
				//cout<<"het["<<iloc<<"] = "<<het<<"\n";
				//f10<<het<<"\n";
			}
			this->cal_snpstatRedacc(gr, numsnp, x, w);			
		}
	}
	//f10.close();
	//cout<<"fin de cal_snhet\n";
	//exit(3);
}


void ParticleC::cal_snnei(int gr, int numsnp) {
	//cout<<"debut de cal_snnei\n";
	long double x = 0, w;
	long double nei, fi, fj, gi, gj;
	int iloc, loc, n1, n2;
	auto&& samp = this->grouplist[gr].sumstatsnp[numsnp].samp.get();
	int sample = samp[0];
	int sample1 = samp[1];
	//cout<<"numsnp="<<numsnp<<"sample = "<<sample<<"    sample1 = "<<sample1<<"\n";
	this->grouplist[gr].sumstatsnp[numsnp].n = this->grouplist[gr].nloc;
	for (iloc = 0; iloc < this->grouplist[gr].nloc; iloc++) {
		loc = this->grouplist[gr].loc[iloc];
		w = this->locuslist[loc].weight;
		if (w > 0.0) {
			n1 = samplesize(loc, sample);
			n2 = samplesize(loc, sample1);
			if ((n1 > 0)and (n2 > 0)) {
				fi = this->locuslist[loc].freq[sample][0];
				fj = this->locuslist[loc].freq[sample1][0];
				gi = this->locuslist[loc].freq[sample][1];
				gj = this->locuslist[loc].freq[sample1][1];
				nei = 1.0 - (fi * fj + gi * gj) / sqrt(fi * fi + gi * gi) / sqrt(fj * fj + gj * gj);
				x = nei;
					//if (iloc<10) cout<<"numsnp="<<numsnp <<"    nei["<<iloc<<"]="<<nei<<"\n";
				}
			this->cal_snpstatRedacc(gr, numsnp, x, w);						
		}
	}
}
void ParticleC::cal_snfl(int gr, int numsnp, int npop) {
	int iloc,loc, npopr;
	long double w;
	if (npop == 0) {
		npopr = this->nsample;
	} else npopr = min(npop,this->nsample);
	auto& statC = grouplist[gr].sumstatsnp[numsnp];
	vector<int> samples;
	if (npopr == nsample) {
		samples.resize(npopr);
		iota(samples.begin(),samples.end(),0);
	} else {
		samples = 	grouplist[gr].sumstatsnp[numsnp].samp.get();
		samples.resize(npopr);
	}
	grouplist[gr].sumstatsnp[numsnp].n = grouplist[gr].nloc;
	for (iloc = 0; iloc < grouplist[gr].nloc; iloc++) {
		loc = grouplist[gr].loc[iloc];
		w = locuslist[loc].weight;
		if (w > 0.0) {
			auto freq0 = locuslist[loc].freq[samples[0]][0];
			bool fixed = freq0 == 0.0 || freq0 == 1.0;
			if (fixed) {
				if (samples.size() > 1) {
					for (auto i = 1; i < samples.size(); i++) 
						fixed = fixed && locuslist[loc].freq[samples[i]][0] == freq0;
				}
			}
			statC.sw++;
			if (fixed) statC.sw0++;
		}
	}

}
void ParticleC::cal_snfstd(int gr, int numsnp, int npop) {
	long double numt = 0, dent = 0, w;
	long double num = 0, den = 0, x = 0;
	int iloc,loc, npopr;
	if (npop == 0) {
		npopr = this->nsample;
	} else npopr = min(npop,this->nsample);
	auto& statC = grouplist[gr].sumstatsnp[numsnp];
	vector<int> samples;
	if (npopr == nsample) {
		samples.resize(npopr);
		iota(samples.begin(),samples.end(),0);
	} else {
		samples = 	grouplist[gr].sumstatsnp[numsnp].samp.get();
		samples.resize(npopr);
	}
	grouplist[gr].sumstatsnp[numsnp].n = grouplist[gr].nloc;
	for (iloc = 0; iloc < grouplist[gr].nloc; iloc++) {
		loc = grouplist[gr].loc[iloc];
		w = locuslist[loc].weight;
		if (w > 0.0) {
			if (grouplist[gr].type == 3) {// POOLSEQ
				long double R_1 = 0.0, R_2 = 0.0;
				long double C_1 = 0.0, C_1_star = 0.0, n_c;
				long double MSI, MSP;
				long double pi1 = 0.0, pi2 = 0.0;
				long double SSI = 0.0, SSP = 0.0;
				for (auto r: samples)
				{	
					auto c = static_cast<long double>(samplesize(loc, r));
					auto n = static_cast<long double>(locuslist[loc].nreads[r]);
					auto x1 = locuslist[loc].nreads1[r];
					auto x2 = n - locuslist[loc].nreads1[r];
					if (c > 0.0) {
						C_1 += n/c + (c -1)/c;
						C_1_star += n * ( n/c + (c-1)/c);
						R_1 += n;
						R_2 += n * n;
						if (n > 0.0)
						{
							pi1 += x1;
							pi2 += x2;
							SSI += x1 - (x1 * x1) / n + x2 - (x2 * x2) / n;
						}
					}
				}
				pi1 /= R_1;
				pi2 /= R_1;
				C_1_star /= R_1;
				for (auto r: samples)
				{
					auto n = static_cast<long double>(locuslist[loc].nreads[r]);
					auto x1 = locuslist[loc].nreads1[r];
					auto x2 = n - locuslist[loc].nreads1[r];
					if (n > 0.0) {
						auto r1 = x1 / n - pi1;
						auto r2 = x2 / n - pi2;
						SSP += n * (r1 * r1 + r2 * r2);
					}
				}
				n_c = (R_1 - R_2/R_1)/(C_1 - C_1_star);
				MSI = SSI /(R_1 - C_1);
				MSP = SSP /(C_1 - C_1_star);
				num = (MSP - MSI);
				den = (MSP + (n_c - 1.0)*MSI);
			} else {  //SNP
				long double S_1 = 0, S_2 = 0;
				long double n_c, MSI, MSP;
				long double SSI = 0.0, SSP = 0.0;
				vector<long double> p(2);
				vector<long double> pi_hat { 0.0, 0.0 };
				int n0 = 0;
				for (auto r : samples)
				{
					auto n = static_cast<long double>(samplesize(loc, r));
					if (n > 0)
					{
						S_1 += n;
						S_2 += n * n;
						for(auto k = 0; k < 2; k++) {
							p[k] = locuslist[loc].freq[r][k];
							pi_hat[k] += n * p[k];
							SSI += n * p[k] * (1 - p[k]);
						}
					} else n0++;
				}
				for(auto k = 0; k < 2; k++) pi_hat[k] /= S_1;
				for (auto r : samples)
				{
					auto n = static_cast<long double>(samplesize(loc, r));
					if (n > 0)
					{
						for(auto k = 0; k < 2; k++) {
							p[k] = locuslist[loc].freq[r][k];
							auto r = p[k] - pi_hat[k];
							SSP += n * (r * r);
						}						
					}
				}
				auto n_d = static_cast<long double>(npopr - n0);
				n_c = (S_1 - S_2 / S_1) / (n_d - 1.0);
				MSI = SSI / (S_1 - n_d);
				MSP = SSP / (n_d - 1.0);
				num = (MSP - MSI);
				den = (MSP + (n_c - 1.0) * MSI);

				// if (iloc < 10 && samples[0] == 0 && samples[1] == 1 && samples.size() == 2) 
				// 	cout << "[FDC] locus " << iloc << " : " << (MSP - MSI)/(MSP + (n_c - 1.0) * MSI) << " " << (MSP - MSI) << " " << (MSP + (n_c - 1.0) * MSI) << endl;
			}
			if (abs(den) > 0.0) x = num/den;
			cal_snpstatRedacc(gr, numsnp, x , w);
			numt += num;
			dent += den;
		}
	}
	statC.mx = numt/dent;
}

void ParticleC::cal_snfst(int gr, int numsnp) {
	//cout<<"debut de cal_snfst\n";
	long double x = 0, w;
	long double fst, fi, fj, n1, n2;
	long double sniA, sniAA, sni, sni2, s2A, s1l, s3l, nc, MSI, MSP, s2I, s2P;
	int iloc, loc;
	auto&& samp = this->grouplist[gr].sumstatsnp[numsnp].samp.get();
	int sample = samp[0];
	int sample1 = samp[1];
	this->grouplist[gr].sumstatsnp[numsnp].n = this->grouplist[gr].nloc;
	for (iloc = 0; iloc < this->grouplist[gr].nloc; iloc++) {
		loc = this->grouplist[gr].loc[iloc];
		w = this->locuslist[loc].weight;
		if (w > 0.0) {
			n1 = (long double)samplesize(loc, sample);
			n2 = (long double)samplesize(loc, sample1);
			if ((n1 > 0)and (n2 > 0)) {
				s1l = 0.0;
				s3l = 0.0;
				for (int all = 0; all < 2; all++) {
					fi = this->locuslist[loc].freq[sample][all];
					fj = this->locuslist[loc].freq[sample1][all];
					sniAA = fi * n1 + fj * n2;
					sniA = 2.0 * sniAA;
					sni = n1 + n2;
					sni2 = n1 * n1 + n2 * n2;
					s2A = 2.0 * (n1 * fi * fi + n2 * fj * fj);
					nc = (sni - sni2 / sni);
					MSI = (0.5 * sniA + sniAA - s2A) / (sni - 2.0);
					MSP = (s2A - 0.5 * sniA * sniA / sni);
					s2I = 0.5 * MSI;
					s2P = (MSP - MSI) / 2.0 / nc;
					s1l += s2P;
					s3l += s2P + s2I;
				}
				fst = 0.0;
				if ((s3l > 0.0)and (s1l > 0.0)) fst = s1l / s3l;
				if (fst < 0.0) fst = 0.0;
				x = fst;
				if (iloc < 10 && sample == 0 && sample1 == 1) 
					cout << "[JMC] locus " << iloc << " : " << fst << " " << s1l << " " << s3l << endl;

				//printf("s1l=%10.5Lf   s3l=%10.5Lf   fst[%3d] = %10.5Lf\n", s1l,s3l,nl,fst);
			}
			this->cal_snpstatRedacc(gr, numsnp, x, w);			
		}
	}
}


void ParticleC::cal_snaml(int gr, int numsnp) {
	//cout<<"debut de cal_snaml\n";
	long double x = 0, w;
	long double aml, f1, f2, f3;
	int iloc, loc;
	auto&& samp = this->grouplist[gr].sumstatsnp[numsnp].samp.get();
	int sample = samp[0];
	int sample1 = samp[1];
	int sample2 = samp[2];
	this->grouplist[gr].sumstatsnp[numsnp].n = this->grouplist[gr].nloc;
	for (iloc = 0; iloc < this->grouplist[gr].nloc; iloc++) {
		loc = this->grouplist[gr].loc[iloc];
		w = this->locuslist[loc].weight;
		if (w > 0.0) {
			if ((samplesize(loc, sample1) > 0)and (samplesize(loc, sample2) > 0)) {
				f1 = this->locuslist[loc].freq[sample1][0];
				f2 = this->locuslist[loc].freq[sample2][0];
				f3 = this->locuslist[loc].freq[sample][0];
				if (f1 != f2) {
					aml = (f3 - f2) / (f1 - f2);
					if ((aml < 0.0)or (aml > 1.0)) {
//						w = 0.0;
						if (aml < 0) x = 0.0;
						else x = 1.0;
					}
					else x = aml;
				} else { x = 0.5; w = 0.0; }
			}
			this->cal_snpstatRedacc(gr, numsnp, x, w);			
		}
		//if (this->grouplist[gr].sumstatsnp[numsnp].w[iloc]>0.0) cout<<this->grouplist[gr].sumstatsnp[numsnp].x[iloc]<<"\n";
	}
	//cout<<"cal_snaml  sw="<<this->grouplist[gr].sumstatsnp[numsnp].sw<<"\n";
}

void ParticleC::cal_snf3r(int gr, int numsnp)
{
	long double x = 0, w;
	long double f1, f2, f3;
	int iloc, loc;
	auto&& samp = this->grouplist[gr].sumstatsnp[numsnp].samp.get();
	int sample = samp[0];
	int sample1 = samp[1];
	int sample2 = samp[2];
	this->grouplist[gr].sumstatsnp[numsnp].n = this->grouplist[gr].nloc;
	for (iloc = 0; iloc < this->grouplist[gr].nloc; iloc++)
	{
		loc = this->grouplist[gr].loc[iloc];
		w = this->locuslist[loc].weight;
		if (w > 0.0)
		{
			if ((samplesize(loc, sample1) > 0) and (samplesize(loc, sample2) > 0))
			{
				int np = samplesize(loc, sample);
				if (this->grouplist[gr].type == 3)
				{ //Poolseq
					double a1p = this->locuslist[loc].nreads1[sample];
					double c1p = this->locuslist[loc].nreads[sample];
					double a2p = this->locuslist[loc].nreads1[sample1];
					double c2p = this->locuslist[loc].nreads[sample1];
					double a3p = this->locuslist[loc].nreads1[sample2];
					double c3p = this->locuslist[loc].nreads[sample2];
					//						double alpha = np * a1p * ( a1p - c1p ) / ( np - 1 ) / c1p / (c1p - 1);
					double alpha = ((np * a1p * (a1p - 1)) / (c1p * (c1p - 1)) - (a1p / c1p)) / (np - 1);
					double betaBC = (a2p * a3p) / (c2p * c3p);
					double betaAB = (a1p * a2p) / (c1p * c2p);
					double betaAC = (a1p * a3p) / (c1p * c3p);
					x = alpha + betaBC - betaAB - betaAC;
				}
				else
				{ // Genepop/SNP
					f1 = this->locuslist[loc].freq[sample][0];
					f2 = this->locuslist[loc].freq[sample1][0];
					f3 = this->locuslist[loc].freq[sample2][0];
					//double alpha = f1 * (1 - f1) * np / (np - 1);
					// version JMC + MAT corrected
					double alpha = f1 * (1 - f1) * 1 / (np - 1);
					x = (f1 - f2) * (f1 - f3) - alpha;
				}
			}
			this->cal_snpstatRedacc(gr, numsnp, x, w);						
		}
	}
}

void ParticleC::cal_snf4r(int gr, int numsnp)
{
	long double x = 0, w;
	long double a,b,c,d;
	int iloc, loc;
	auto&& samp = this->grouplist[gr].sumstatsnp[numsnp].samp.get();
	int sample = samp[0];
	int sample1 = samp[1];
	int sample2 = samp[2];
	int sample3 = samp[3];
	this->grouplist[gr].sumstatsnp[numsnp].n = this->grouplist[gr].nloc;
	for (iloc = 0; iloc < this->grouplist[gr].nloc; iloc++)
	{
		loc = this->grouplist[gr].loc[iloc];
		w = this->locuslist[loc].weight;
		if (w > 0.0)
		{
			if ((samplesize(loc, sample1) > 0) and (samplesize(loc, sample2) > 0) and (samplesize(loc, sample3) > 0))
			{
				int np = samplesize(loc, sample);
				if (this->grouplist[gr].type == 3)
				{ //Poolseq
					double a1p = this->locuslist[loc].nreads1[sample];
					double c1p = this->locuslist[loc].nreads[sample];
					double a2p = this->locuslist[loc].nreads1[sample1];
					double c2p = this->locuslist[loc].nreads[sample1];
					double a3p = this->locuslist[loc].nreads1[sample2];
					double c3p = this->locuslist[loc].nreads[sample2];
					double a4p = this->locuslist[loc].nreads1[sample3];
					double c4p = this->locuslist[loc].nreads[sample3];
					a = (a1p / c1p);
					b = (a2p / c2p);
					c = (a3p / c3p);
					d = (a4p / c4p);
				}
				else
				{
					a = this->locuslist[loc].freq[sample][0];
					b = this->locuslist[loc].freq[sample1][0];
					c = this->locuslist[loc].freq[sample2][0];
					d = this->locuslist[loc].freq[sample3][0];
				}
				x = (a - b) * (c - d);
			}
			this->cal_snpstatRedacc(gr, numsnp, x, w);			
		}
	}
}

////////////// Début Q1/Q2
double ParticleC::q1(int gr, int loc, int sample, bool bias) {
	double x = 0;
	auto np = static_cast<double>(samplesize(loc, sample));
	if (np > 1.0)
	{
		if (bias)
		{
			auto p1 = this->locuslist[loc].freq[sample][0];
			auto p2 = this->locuslist[loc].freq[sample][1];
			x = p1 * p1 + p2 * p2;
		}
		else if (this->grouplist[gr].type == 3)
		{ //Poolseq
			auto r1 = this->locuslist[loc].nreads1[sample];
			auto c1 = this->locuslist[loc].nreads[sample];
			auto r2 = c1 - r1;
			auto s1 = r1 * (r1 - 1);
			auto s2 = r2 * (r2 - 1);
			x = ((np / (c1 * (c1 - 1))) * (s1 + s2) - 1) / (np - 1);
		}
		else // snp
		{
			auto y1 = this->locuslist[loc].freq[sample][0] * np;
			auto y2 = this->locuslist[loc].freq[sample][1] * np;
			x = (y1 * (y1 - 1) + y2 * (y2 - 1)) / (np * (np - 1));
		}
	}
	return x;
}

double ParticleC::q2(int gr, int loc, const vector<int>& samp, bool bias) {
	double x = 0;
	auto sample = samp[0];
	auto sample1 = samp[1];
	auto n1 = static_cast<double>(samplesize(loc, sample));
	auto n2 = static_cast<double>(samplesize(loc, sample1));
	if ((n1 > 0) and (n2 > 0))
	{
		if (bias)
		{
			double p11 = locuslist[loc].freq[sample][0];
			double p21 = locuslist[loc].freq[sample1][0];
			double p12 = 1 - p11;
			double p22 = 1 - p21;
			x = (p11 * p12 + p21 * p22);
		}
		else if (grouplist[gr].type == 3)
		{ //Poolseq
			double r11 = locuslist[loc].nreads1[sample];
			double c1 = locuslist[loc].nreads[sample];
			double r12 = c1 - r11;
			double r21 = locuslist[loc].nreads1[sample1];
			double c2 = locuslist[loc].nreads[sample1];
			double r22 = c2 - r21;
			x = (r11 * r21 + r12 * r22) / (c1 * c2);
		}
		else
		{ // Genepop/SNP
			double y11 = locuslist[loc].freq[sample][0] * n1;
			double y21 = locuslist[loc].freq[sample1][0] * n2;
			double y12 = n1 - y11;
			double y22 = n2 - y21;
			x = (y11 * y21 + y12 * y22) / (n1 * n2);
		}
	}
	return x;
}

void ParticleC::cal_snhw(int gr, int numsnp) {
	int iloc, loc;
	long double x = 0, w;
	auto&& samp = grouplist[gr].sumstatsnp[numsnp].samp.get();
	grouplist[gr].sumstatsnp[numsnp].n = grouplist[gr].nloc;
	for (iloc = 0; iloc < grouplist[gr].nloc; iloc++) {
		loc = grouplist[gr].loc[iloc];
		w = locuslist[loc].weight;
		if (w > 0.0) {
			x = 1 - q1(gr, loc, samp[0], false);
			cal_snpstatRedacc(gr, numsnp, x, w);
		}
	}
}

void ParticleC::cal_snhb(int gr, int numsnp) {
	int iloc, loc;
	long double x = 0, w;
	auto&& samp = grouplist[gr].sumstatsnp[numsnp].samp.get();
	grouplist[gr].sumstatsnp[numsnp].n = grouplist[gr].nloc;
	for (iloc = 0; iloc < grouplist[gr].nloc; iloc++) {
		loc = grouplist[gr].loc[iloc];
		w = locuslist[loc].weight;
		if (w > 0.0) {
			x = 1 - q2(gr, loc, samp, false);
			cal_snpstatRedacc(gr, numsnp, x, w);
		}
	}
}

int ParticleC::findstat(string&& statname){
	int statnum = 0;
	bool statfound = false;
	while (statnum < stats.size() && !statfound) {
		statfound = stats[statnum++].name == statname;
	}
	return statfound ? statnum - 1 : -1;
}

void ParticleC::cal_snfsti(int gr, int numsnp)
{
	StatsnpC& stsnp = grouplist[gr].sumstatsnp[numsnp];
	auto sample = stsnp.samp.get()[0];
	long double Hbmoy = 0.0, Hw, Hwv;
	int Hbstat = findstat("HBm");
	int Hwstat = findstat("HWm");
	int Hwvstat = findstat("HWv");
	int nQb = 0;
	for(int i = 0; i < grouplist[gr].sumstat.size(); i++) {
		auto&& samp = grouplist[gr].sumstat[i].samp.get();
		auto curstat = grouplist[gr].sumstat[i].cat;
		auto val = grouplist[gr].sumstat[i].val;
		if (curstat == Hbstat) {
			Hbmoy += val;
			nQb++;				
		}
		if (samp.size() == 1 && samp[0] == sample) {
			if (curstat == Hwstat) Hw = val;
			if (curstat == Hwvstat) Hwv = val;
		}
	}
	Hbmoy /= static_cast<double>(nQb);
	stsnp.mx = 1.0 - Hw / Hbmoy;
	stsnp.mx2 = Hwv / (Hbmoy * Hbmoy);
	stsnp.sw = 1.0;
	stsnp.sw2 = 0.0;
}


////////////// Fin Stats SNP

long double ParticleC::cal_pid1p(int gr, int st) {
	int iloc, kloc, nt = 0, ni = 0, cat;
	auto&& samp = this->grouplist[gr].sumstat[st].samp.get();
	int sample = samp[0];
	for (iloc = 0; iloc < this->grouplist[gr].nloc; iloc++) {
		kloc = this->grouplist[gr].loc[iloc];
		cat = this->locuslist[kloc].type % 5;
		if (samplesize(kloc, sample) > 1) {
			for (int i = 0; i < dataobs.ssize[cat][sample] - 1; i++) {
				if (this->locuslist[kloc].haplomic[sample][i] != MICMISSING) {
					for (int j = i + 1; j < dataobs.ssize[cat][sample]; j++) {
						if (this->locuslist[kloc].haplomic[sample][j] != MICMISSING) {
							nt++;
							if (this->locuslist[kloc].haplomic[sample][i] == this->locuslist[kloc].haplomic[sample][j]) ni++;
						}
					}
				}
			}
		}
	}
	//cout<<ni<<" sur "<<nt<<"\n";
	if (nt > 0) return (long double)ni / (long double)nt;
	else return 0.0;
}

long double ParticleC::cal_nal1p(int gr, int st) {
	long double nal, nalm = 0.0;
	int iloc, kloc, nl = 0;
	auto&& samp = this->grouplist[gr].sumstat[st].samp.get();
	int sample = samp[0];
	//cout <<"groupe "<<gr<<"  cat "<<this->grouplist[gr].sumstat[st].cat<<"   sample "<<this->grouplist[gr].sumstat[st].samp;//<<"\n";
	//cout << "     nloc = " << this->grouplist[gr].nloc;// <<"\n";
	for (iloc = 0; iloc < this->grouplist[gr].nloc; iloc++) {
		kloc = this->grouplist[gr].loc[iloc];
		//cout << this->locuslist[loc].samplesize[stat.samp] <<"\n";
		nal = 0.0;
		if (samplesize(kloc, sample) > 0) {
			for (int k = 0; k < this->locuslist[kloc].nal; k++) {
				if (this->locuslist[kloc].freq[sample][k] > 0.00001) {
					nalm += 1.0;
					nal += 1.0;
				}
			}
			nl++;
		}
		//cout <<"locus "<<kloc<<"    nal="<<nal<<"   nl="<<nl<<"\n";//<<"    nmoy="<<  "\n";
	}
	if (nl > 0) nalm = nalm / (long double)nl;
	//cout<<"    nmoy="<<nalm<<  "\n";
	return nalm;
}

long double ParticleC::cal_nal2p(int gr, int st) {
	int k, iloc, kloc, nl = 0;
	long double nalm = 0.0;
	auto&& samp = this->grouplist[gr].sumstat[st].samp.get();
	int sample = samp[0];
	int sample1 = samp[1];
	for (iloc = 0; iloc < this->grouplist[gr].nloc; iloc++) {
		kloc = this->grouplist[gr].loc[iloc];
		if ((samplesize(kloc, sample) > 0)and (samplesize(kloc, sample1) > 0)) {
			for (int k = 0; k < this->locuslist[kloc].nal; k++) {
				if (this->locuslist[kloc].freq[sample][k] + this->locuslist[kloc].freq[sample1][k] > 0.000001) nalm += 1.0;
			}
			nl++;
		}
		else {
			if (samplesize(kloc, sample) > 0) {
				for (k = 0; k < this->locuslist[kloc].nal; k++) {
					if (this->locuslist[kloc].freq[sample][k] > 0.000001) nalm += 1.0;
				}
				nl++;
			}
			else {
				if (samplesize(kloc, sample1) > 0) {
					for (k = 0; k < this->locuslist[kloc].nal; k++) {
						if (this->locuslist[kloc].freq[sample1][k] > 0.000001) nalm += 1.0;
					}
					nl++;
				}
			}
		}
	}
	if (nl > 0) nalm = nalm / (long double)nl;
	return nalm;
}

long double ParticleC::cal_het1p(int gr, int st) {
	long double het, hetm = 0.0;
	int iloc, kloc, nl = 0, sasize;
	auto&& samp = this->grouplist[gr].sumstat[st].samp.get();
	int sample = samp[0];
	for (iloc = 0; iloc < this->grouplist[gr].nloc; iloc++) {
		kloc = this->grouplist[gr].loc[iloc];
		sasize = samplesize(kloc, sample);
		if (sasize > 1) {
			het = 1.0;
			for (int k = 0; k < this->locuslist[kloc].nal; k++) het -= sqr(this->locuslist[kloc].freq[sample][k]);
			het *= (long double)sasize / (long double)(sasize - 1);
			hetm += het;
			//printf("kloc=%2d   het = %13.10Lf\n",kloc,het);
			nl++;
		}
	}
	if (nl > 0) hetm = hetm / (long double)nl;
	//printf("\n heterozygotie = %13.10Lf\n\n",hetm);
	return hetm;
}

long double ParticleC::cal_het2p(int gr, int st) {
	long double het, hetm = 0.0;
	int iloc, loc, n1, n2, nt, nl = 0;
	auto&& samp = this->grouplist[gr].sumstat[st].samp.get();
	int sample = samp[0];
	int sample1 = samp[1];
	for (iloc = 0; iloc < this->grouplist[gr].nloc; iloc++) {
		loc = this->grouplist[gr].loc[iloc];
		n1 = samplesize(loc, sample);
		n2 = samplesize(loc, sample1);
		nt = n1 + n2;
		if (nt > 1) {
			het = 1.0;
			for (int k = 0; k < this->locuslist[loc].nal; k++) {
				het -= sqr(((long double)n1 * this->locuslist[loc].freq[sample][k] + (long double)n2 * this->locuslist[loc].freq[sample1][k]) / (long double)nt);
			}
			het *= (long double)nt / (long double)(nt - 1);
			hetm += het;
			nl++;
		}
	}
	if (nl > 0) return hetm / (long double)nl;
	else return 0.0;
}

long double ParticleC::cal_var1p(int gr, int st) {
	long double s, v, vm = 0.0, m, ms;
	int iloc, loc, n, nl = 0, cat;
	auto&& samp = this->grouplist[gr].sumstat[st].samp.get();
	int sample = samp[0];
	for (iloc = 0; iloc < this->grouplist[gr].nloc; iloc++) {
		loc = this->grouplist[gr].loc[iloc];
		cat = this->locuslist[loc].type % 5;
		v = 0.0;
		s = 0.0;
		n = 0;
		for (int i = 0; i < dataobs.ssize[cat][sample]; i++) {
			if (this->locuslist[loc].haplomic[sample][i] != MICMISSING) {
				n++;
				s += this->locuslist[loc].haplomic[sample][i];
				v += sqr(this->locuslist[loc].haplomic[sample][i]);
			}
		}
		m = (long double)n;
		ms = (long double)this->locuslist[loc].motif_size;
		if (n > 1) {
			vm += (v - sqr(s) / m) / (m - 1.0) / sqr(ms);
			nl++;
		}
		//cout<<"var["<<sample<<"]["<<loc<<"]="<<(v-sqr(s)/m)/(m-1.0)/sqr(ms)<<"   v="<<v<<"   s="<<s<<"   m="<<m<<"\n";
	}
	//printf("\nvm=%13.10Lf\n",vm/(long double)nl);
	if (nl > 0) return vm / (long double)nl;
	else return 0.0;
}

long double ParticleC::cal_var2p(int gr, int st) {
	long double s, v, vm = 0.0, m, ms;
	int iloc, loc, n, nl = 0, cat;
	auto&& samp = this->grouplist[gr].sumstat[st].samp.get();
	int sample = samp[0];
	int sample1 = samp[1];
	for (iloc = 0; iloc < this->grouplist[gr].nloc; iloc++) {
		loc = this->grouplist[gr].loc[iloc];
		cat = this->locuslist[loc].type % 5;
		v = 0.0;
		s = 0.0;
		n = 0;
		for (int i = 0; i < dataobs.ssize[cat][sample]; i++) {
			if (this->locuslist[loc].haplomic[sample][i] != MICMISSING) {
				n++;
				s += this->locuslist[loc].haplomic[sample][i];
				v += sqr(this->locuslist[loc].haplomic[sample][i]);
			}
		}
		for (int i = 0; i < dataobs.ssize[cat][sample1]; i++) {
			if (this->locuslist[loc].haplomic[sample1][i] != MICMISSING) {
				n++;
				s += this->locuslist[loc].haplomic[sample1][i];
				v += sqr(this->locuslist[loc].haplomic[sample1][i]);
			}
		}
		m = (long double)n;
		ms = (long double)this->locuslist[loc].motif_size;
		if (n > 1) {
			vm += (v - sqr(s) / m) / (m - 1.0) / sqr(ms);
			nl++;
		}
	}
	if (nl > 0) return vm / (long double)nl;
	else return 0.0;
}

long double ParticleC::cal_mgw1p(int gr, int st) {
	long double num = 0.0, den = 0.0, un = 1.0;
	int min, max;
	int iloc, loc;
	auto&& samp = this->grouplist[gr].sumstat[st].samp.get();
	int sample = samp[0];
	for (iloc = 0; iloc < this->grouplist[gr].nloc; iloc++) {
		loc = this->grouplist[gr].loc[iloc];
		if (samplesize(loc, sample) > 0) {
			for (int k = 0; k < this->locuslist[loc].nal; k++) {
				if (this->locuslist[loc].freq[sample][k] > 0.00001) num += 1.0;
			}
			min = 0;
			while (this->locuslist[loc].freq[sample][min] < 0.00001) min++;
			max = this->locuslist[loc].kmax - this->locuslist[loc].kmin;
			while (this->locuslist[loc].freq[sample][max] < 0.00001) max--;
			den += un + (long double)(max - min) / (long double)this->locuslist[loc].motif_size;
		}
	}
	if (den > 0) return num / den;
	else return 0.0;
}

long double ParticleC::cal_Fst2p(int gr, int st) {
	long double s1 = 0.0, s3 = 0.0, s1l, s2l, s3l, sni, sni2, sniA, sniAA;
	int al, pop, kpop, i, ig1, ig2, nA, AA, ni;
	long double s2A, MSG, MSI, MSP, s2G, s2I, s2P, nc;
	int iloc, loc, nn, ind, cat;
	auto&& samp = this->grouplist[gr].sumstat[st].samp.get();
	int sample = samp[0];
	int sample1 = samp[1];
	//cout << "cal_Fst2p\n";
	for (iloc = 0; iloc < this->grouplist[gr].nloc; iloc++) {
		loc = this->grouplist[gr].loc[iloc];
		cat = this->locuslist[loc].type % 5;
		//cout << "Fst loc = "<<loc<<"\n";
		s1l = 0.0;
		s2l = 0.0;
		s3l = 0.0;
		nc = 0;
		//		    cout << "\n\nlocus "<< loc << "\n";
		for (al = this->locuslist[loc].kmin; al <= this->locuslist[loc].kmax; al++) {
			//		    	cout << "\n allele "<< al << "\n";
			sni = 0;
			sni2 = 0;
			sniA = 0;
			sniAA = 0;
			s2A = 0.0;
			for (kpop = 0; kpop < 2; kpop++) {
				if (kpop == 0) pop = sample;
				else pop = sample1;
				nA = 0;
				AA = 0;
				ni = 0;
				ind = 0;
				//ind est le numéro de l'individu (haploïde ou diploïde ou même sans génotype si Y femelle)
				//i est le numéro de la copie du gène
				for (i = 0; i < dataobs.ssize[cat][pop];) {
					nn = calploidy(loc, pop, ind);
					ind++;
					switch (nn) {
						case 1: ig1 = this->locuslist[loc].haplomic[pop][i];
							i++;
							if (ig1 == al) {
								AA++;
								nA += 2;
							}
							if (ig1 != MICMISSING) ni++;

							break;
						case 2: ig1 = this->locuslist[loc].haplomic[pop][i];
							i++;
							ig2 = this->locuslist[loc].haplomic[pop][i];
							i++;
							if ((ig1 == al) and (ig2 == al)) AA++;
							if (ig1 == al) nA++;
							if (ig2 == al) nA++;
							if ((ig1 != MICMISSING)and (ig2 != MICMISSING)) ni++;
							break;
					}
				}
				sniA += (long double)nA;
				sniAA += (long double)AA;
				sni += (long double)ni;

				sni2 += sqr((long double)ni);
				s2A += sqr((long double)nA) / (long double)(2 * ni);
			}
			if (sni > 0) nc = sni - sni2 / sni;
			else nc = 0.0;
			//		    	cout << "Allele " << al <<"   sni = " << sni << "   sni2 = "<< sni2 << "   nc = "<< nc;
			//		    	cout << "   sniA = "<< sniA << "   sniAA = " << sniAA<< "   s2A = " << s2A << "\n";
			if (sni * nc > 0.0) {
				MSG = (0.5 * sniA - sniAA) / sni;
				MSI = (0.5 * sniA + sniAA - s2A) / (sni - 2.0);
				MSP = (s2A - 0.5 * sqr(sniA) / sni);
				s2G = MSG;
				s2I = 0.5 * (MSI - MSG);
				s2P = (MSP - MSI) / (2.0 * nc);
				s1l = s1l + s2P; //numerateur de theta
				s2l = s2l + s2P + s2I; //numerateur de F
				s3l = s3l + s2P + s2I + s2G; //denominateur de theta et de F
			}
		}
		//            cout << "s1l = "<< s1l<<"   s2l = "<<s2l<<"   s3l = "<< s3l << "\n";
		s1 += (long double)nc * s1l;
		s3 += (long double)nc * s3l;
	}
	if (s3 > 0.0) return s1 / s3;
	else return 0.0;
}

long double ParticleC::cal_lik2p(int gr, int st) {
	int pop, nal, k, i, ind, ig1, ig2;
	long double frt, a, b, num_lik, den_lik, lik;
	int iloc, loc, nn, nl = 0, cat, sasize, sasize1;
	auto&& samp = this->grouplist[gr].sumstat[st].samp.get();
	int sample = samp[0];
	int sample1 = samp[1];
	lik = 0.0;
	for (iloc = 0; iloc < this->grouplist[gr].nloc; iloc++) {
		loc = this->grouplist[gr].loc[iloc];
		cat = this->locuslist[loc].type % 5;
		sasize = samplesize(loc, sample);
		sasize1 = samplesize(loc, sample1);
		if (sasize * sasize1 > 0) {
			nl++;
			nal = 0;
			for (k = 0; k < this->locuslist[loc].nal; k++) {
				frt = 0.0;
				for (pop = 0; pop < this->nsample; pop++) frt += this->locuslist[loc].freq[pop][k];
				if (frt > 0.000001) nal++;
			}
			b = 1.0 / (long double)nal;
			a = 1.0 / (long double)dataobs.nind[sample];
			ind = 0;
			for (i = 0; i < dataobs.ssize[cat][sample]; i++) {
				nn = calploidy(loc, sample, ind);
				ind++;
				switch (nn) {
					case 1: ig1 = this->locuslist[loc].haplomic[sample][i];
						if (ig1 != MICMISSING) {
							num_lik = round(this->locuslist[loc].freq[sample1][ig1 - this->locuslist[loc].kmin] * sasize1) + b;
							den_lik = (long double) (sasize1 + 1);
							lik -= a * log10(num_lik / den_lik);
						}
						break;
					case 2: ig1 = this->locuslist[loc].haplomic[sample][i];
						ig2 = this->locuslist[loc].haplomic[sample][++i];
						if ((ig1 != MICMISSING)and (ig2 != MICMISSING)) {
							if (ig1 == ig2) {
								num_lik = 1.0 + b + arrondi(this->locuslist[loc].freq[sample1][ig1 - this->locuslist[loc].kmin] * sasize1);
								num_lik *= b + arrondi(this->locuslist[loc].freq[sample1][ig2 - this->locuslist[loc].kmin] * sasize1);
							}
							else {
								num_lik = 2.0 * (b + arrondi(this->locuslist[loc].freq[sample1][ig1 - this->locuslist[loc].kmin] * sasize1));
								num_lik *= b + arrondi(this->locuslist[loc].freq[sample1][ig2 - this->locuslist[loc].kmin] * sasize1);
							}
							den_lik = (long double) ((sasize1 + 2) * (sasize1 + 1));
							lik -= a * log10(num_lik / den_lik);
						}
						break;
				}
			}
		}
	}
	if (nl > 0) return lik / (long double)nl;
	else return 0.0;
}

long double ParticleC::cal_dmu2p(int gr, int st) {
	long double *moy, dmu2 = 0.0, s;
	int iloc, loc, pop, nl = 0, cat, sasize, sasize1;
	auto&& samp = this->grouplist[gr].sumstat[st].samp.get();
	int sample = samp[0];
	int sample1 = samp[1];
	moy = new long double[2];
	for (iloc = 0; iloc < this->grouplist[gr].nloc; iloc++) {
		loc = this->grouplist[gr].loc[iloc];
		cat = this->locuslist[loc].type % 5;
		sasize = samplesize(loc, sample);
		sasize1 = samplesize(loc, sample1);
		//cout<<"          samplesize["<<loc<<","<<sample<<"] = "<<sasize<<"\n";
		//cout<<"          samplesize["<<loc<<","<<sample1<<"] = "<<sasize1<<"\n";

		if (sasize * sasize1 > 0) {
			nl++;
			for (int kpop = 0; kpop < 2; kpop++) {
				if (kpop == 0) pop = sample;
				else pop = sample1;
				s = 0.0;
				moy[kpop] = 0.0;
				for (int i = 0; i < dataobs.ssize[cat][pop]; i++) {
					if (this->locuslist[loc].haplomic[pop][i] != MICMISSING) s += this->locuslist[loc].haplomic[pop][i];
				}
				moy[kpop] = s / (long double)samplesize(loc, pop);
			}
		}
		dmu2 += sqr((moy[1] - moy[0]) / (long double)this->locuslist[loc].motif_size);
		//cout<<"loc="<<loc<<"   moy["<<sample<<"]="<<moy[0]<<"   moy["<<sample1<<"]="<<moy[1]<< "       delta2="<<sqr((moy[1]-moy[0])/(long double)this->locuslist[loc].motif_size) ;
		//cout<<"   motifsize="<<this->locuslist[loc].motif_size<<"\n";
	}
	//if (loc==17) cout<<"   nl="<<nl<<"  dmu2="<<dmu2/(long double)nl<<"\n\n";
	delete []moy;
	if (nl > 0) return dmu2 / (long double)nl;
	else return 0.0;
}

long double ParticleC::cal_das2p(int gr, int st) {
	long double s = 0.0;
	int iloc, loc, nl = 0, cat;
	auto&& samp = this->grouplist[gr].sumstat[st].samp.get();
	int sample = samp[0];
	int sample1 = samp[1];
	for (iloc = 0; iloc < this->grouplist[gr].nloc; iloc++) {
		loc = this->grouplist[gr].loc[iloc];
		cat = this->locuslist[loc].type % 5;
		for (int i = 0; i < dataobs.ssize[cat][sample]; i++) {
			if (this->locuslist[loc].haplomic[sample][i] != MICMISSING) {
				for (int j = 0; j < dataobs.ssize[cat][sample1]; j++) {
					if (this->locuslist[loc].haplomic[sample1][j] != MICMISSING) {
						nl++;
						if (this->locuslist[loc].haplomic[sample][i] == this->locuslist[loc].haplomic[sample1][j]) s += 1.0;
					}
				}
			}
		}
	}
	if (nl > 0) return s / (long double)nl;
	else return 0.0;
}

pair<long double, long double> ParticleC::pente_lik(int gr, int st, int i0) {
	long double a, freq1, freq2, li0, delta, *li;
	int ig1, ig2, ind, loc, iloc, nn, cat;
	li = new long double[2];
	auto&& samp = this->grouplist[gr].sumstat[st].samp.get();
	int sample = samp[0];
	int sample1 = samp[1];
	int sample2 = samp[2];
	//cout<<"sample ="<<sample<<"   sample1="<<sample1<<"   sample2="<<sample2<<"\n";
	for (int rep = 0; rep < 2; rep++) {
		a = 0.001 * (long double)(i0 + rep);
		li[rep] = 0.0;
		//cout << "rep = " << rep << "   a = "<< a << "  lik = "<< li[rep] << "\n";
		for (iloc = 0; iloc < this->grouplist[gr].nloc; iloc++) {
			loc = this->grouplist[gr].loc[iloc];
			cat = this->locuslist[loc].type % 5;
			//cout << "\n\n locus "<< loc<< "\n";
			ind = 0;
			for (int i = 0; i < dataobs.ssize[cat][sample];) {
				nn = calploidy(loc, sample, ind);
				//cout <<"ploidie="<<nn<<"\n";
				ind++;
				switch (nn) {
					case 1: ig1 = this->locuslist[loc].haplomic[sample][i];
						i++;
						if (ig1 != MICMISSING) {
							freq1 = a * this->locuslist[loc].freq[sample1][ig1 - this->locuslist[loc].kmin] + (1.0 - a) * this->locuslist[loc].freq[sample2][ig1 - this->locuslist[loc].kmin];
							if (freq1 > 0.0) li[rep] += log(freq1);
						}
						break;
					case 2: ig1 = this->locuslist[loc].haplomic[sample][i];
						i++;
						ig2 = this->locuslist[loc].haplomic[sample][i];
						i++;
						if ((ig1 != MICMISSING)and (ig2 != MICMISSING)) {
							//cout<<"ig1="<<ig1<<"    ig2="<<ig2<<"   kmin="<<this->locuslist[loc].kmin<<"\n";
							//cout<<"freq[samp1][ig1]="<<this->locuslist[loc].freq[sample1][ig1-this->locuslist[loc].kmin]<<"\n";
							//cout<<"freq[samp2][ig1]="<<this->locuslist[loc].freq[sample2][ig1-this->locuslist[loc].kmin]<<"\n";
							//cout<<"freq[samp1][ig2]="<<this->locuslist[loc].freq[sample1][ig2-this->locuslist[loc].kmin]<<"\n";
							//cout<<"freq[samp2][ig2]="<<this->locuslist[loc].freq[sample2][ig2-this->locuslist[loc].kmin]<<"\n";
							if (ig1 == ig2) {
								freq1 = a * this->locuslist[loc].freq[sample1][ig1 - this->locuslist[loc].kmin] + (1.0 - a) * this->locuslist[loc].freq[sample2][ig1 - this->locuslist[loc].kmin];
								if (freq1 > 0.0) li[rep] += log(sqr(freq1));
							}
							else {
								freq1 = a * this->locuslist[loc].freq[sample1][ig1 - this->locuslist[loc].kmin] + (1.0 - a) * this->locuslist[loc].freq[sample2][ig1 - this->locuslist[loc].kmin];
								//cout<<"freq1 = "<<freq1<<"\n";
								freq2 = a * this->locuslist[loc].freq[sample1][ig2 - this->locuslist[loc].kmin] + (1.0 - a) * this->locuslist[loc].freq[sample2][ig2 - this->locuslist[loc].kmin];
								//cout<<"freq2 = "<<freq2<<"\n";
								if (freq1 * freq2 > 0.0) li[rep] += log(2.0 * freq1 * freq2);
								else {
									if (freq1 > 0.0) li[rep] += log(freq1);
									if (freq2 > 0.0) li[rep] += log(freq2);
								}
							}
						}
						// cout<<"coucou\n";
						break;
				}
				//cout << "lik = " << li[rep] << "\n";
			}
		}
		//cout << "rep = " << rep << "   a = "<< a << "   lik = "<< li[rep] << "\n";
	}
	//cout << "   li[i0] = "<< li[0] << "\n";
	//cout << "   li[i1] = "<< li[1] << "\n";
	li0 = li[0];
	delta = li[1] - li[0];
	//cout<<"delta ="<<delta<<"\n";
	delete []li;
	return make_pair(li0, delta);
}

long double ParticleC::cal_Aml3p(int gr, int st) {
	int i1 = 1, i2 = 998, i3;
	long double p1, p2, p3, lik1, lik2, lik3;
	pair<long double, long double> c;
	c = pente_lik(gr, st, i1);
	lik1 = c.first;
	p1 = c.second;
	c = pente_lik(gr, st, i2);
	lik2 = c.first;
	p2 = c.second;
	if (abs(lik1) + abs(lik2) < 1.0E-10) return -9999.0;
	if ((p1 < 0.0)and (p2 < 0.0)) return 0.0;
	if ((p1 > 0.0)and (p2 > 0.0)) return 1.0;
	do {
		i3 = (i1 + i2) / 2;
		c = pente_lik(gr, st, i3);
		lik3 = c.first;
		p3 = c.second;
		if (p1 * p3 < 0.0) {
			i2 = i3;
			p2 = p3;
			lik2 = lik3;
		}
		else {
			i1 = i3;
			p1 = p3;
			lik1 = lik3;
		}
	}
	while (abs(i2 - i1) > 1);
	if (lik1 > lik2) return 0.001 * (long double)i1;
	else return 0.001 * (long double)i2;
}

bool ParticleC::identseq(string seq1, string seq2, string appel, int kloc, int sample, int ind) {
	bool ident = true;
	int k = 0;
	if (seq1.length() != seq2.length()) {
		cout << "Erreur au locus " << kloc << " appelé par " << appel;
		cout << "     seq1=" << seq1 << "<(" << seq1.length() << ")   seq2=" << seq2 << "<(" << seq2.length() << ")\n";
		cout << "dnavar = " << this->locuslist[kloc].dnavar << "   sample=" << sample << "   ind=" << ind << "\n";
		exit(1);
	}
	if (seq1.length() == 0) return ident;
	//cout<<seq1<<"    "<<seq2<<"\n";
	while ((ident) and (k < (int)seq1.length())) {
		ident = ((seq1[k] == 'N')or (seq2[k] == 'N')or (seq1[k] == seq2[k]));
		if (ident) k++;
	}
	//if (ident) cout<<"ident\n"; else cout<<"different\n";
	return ident;
}

long double ParticleC::cal_nha1p(int gr, int st) {
	string* haplo;
	int iloc, kloc, j, nhl = 0, nhm = 0, nl = 0, cat;
	bool trouve;
	long double res = 0.0;
	auto&& samp = this->grouplist[gr].sumstat[st].samp.get();
	int sample = samp[0];
	for (iloc = 0; iloc < this->grouplist[gr].nloc; iloc++) {
		nhl = 0;
		kloc = this->grouplist[gr].loc[iloc];
		cat = this->locuslist[kloc].type % 5;
		haplo = new string[dataobs.ssize[cat][sample]];
		if (samplesize(kloc, sample) > 0) {
			nl++;
			if (this->locuslist[kloc].dnavar == 0) {
				haplo[nhl] = "";
				nhl++;
			}
			else {
				for (int i = 0; i < dataobs.ssize[cat][sample]; i++) {
					if (this->locuslist[kloc].haplodna[sample][i] != SEQMISSING) {
						//cout<<"not seqmissing nhl="<<nhl<<"\n";
						if (nhl == 0) {
							haplo[nhl] = this->locuslist[kloc].haplodnavar[sample][i];
							nhl++;
							//cout<<"nl="<<nl<<"   nhl="<<nhl<<"    "<<haplo[nhl-1]<<"\n";
						}
						else {
							trouve = false;
							j = 0;
							while ((not trouve)and (j < nhl)) {
								trouve = identseq(haplo[j], this->locuslist[kloc].haplodnavar[sample][i], "nha1p", kloc, sample, i);
								if (not trouve) j++;
							}
							if (trouve) {
								for (int k = 0; k < this->locuslist[kloc].dnavar; k++) if (haplo[j][k] == 'N') haplo[nhl][k] = this->locuslist[kloc].haplodnavar[sample][i][k];
							}
							else {
								haplo[nhl] = this->locuslist[kloc].haplodnavar[sample][i];
								nhl++;
							}
						}
					}
				}
			}
		}
		nhm += nhl;
		delete []haplo;
	}
	if (nl > 0) res = (long double)nhm / (long double)nl;
	return res;
}

int* ParticleC::cal_nsspl(int kloc, int sample, int* nssl, bool* OK) {
	char c0;
	bool ident;
	int j, nss = 0, cat = this->locuslist[kloc].type % 5;
	vector<int> nuvar;
	if (samplesize(kloc, sample) > 0) {
		*OK = true;
		if ((debuglevel == 9)and (kloc == 10)) cout << "OK=true  dnavar=" << this->locuslist[kloc].dnavar << "\n";
		if (this->locuslist[kloc].dnavar > 0) {
			for (j = 0; j < this->locuslist[kloc].dnavar; j++) {
				if ((debuglevel == 9)and (kloc == 10)) cout << "j=" << j << "  ss=" << dataobs.ssize[cat][sample] << "\n";
				c0 = 'z';
				ident = true;
				for (int i = 0; i < dataobs.ssize[cat][sample]; i++) {
					/*if((debuglevel==9)and(kloc==10)) {
					  cout<<this->locuslist[kloc].haplodnavar[sample][i]<<"   "<<this->locuslist[kloc].haplodnavar[sample][i][j]<<"   "<<c0;
					  cout<<"  "<<(this->locuslist[kloc].haplodna[sample][i]!=SEQMISSING);
					  cout<<"  "<<(this->locuslist[kloc].haplodnavar[sample][i][j]!='N');
					  cout<<"  "<<((this->locuslist[kloc].haplodna[sample][i]!=SEQMISSING)and(this->locuslist[kloc].haplodnavar[sample][i][j]!='N'))<<"\n";
					  }*/
					if ((this->locuslist[kloc].haplodna[sample][i] != SEQMISSING)and (this->locuslist[kloc].haplodnavar[sample][i][j] != 'N')) {
						if ((debuglevel == 9)and (kloc == 10)) cout << "c0=" << c0 << "\n";
						if (c0 == 'z') c0 = this->locuslist[kloc].haplodnavar[sample][i][j];
						else ident = (c0 == this->locuslist[kloc].haplodnavar[sample][i][j]);
						if ((debuglevel == 9)and (kloc == 10)) {
							if (not ident) cout << "nucleotide " << j << "   c0=" << c0 << "   c1=" << this->locuslist[kloc].haplodnavar[sample][i][j] << "   i=" << i << "\n";
						}
					}
					if (not ident) break;
				}
				if (not ident) {
					nss++;
					nuvar.push_back(j);
				}
			}
		}
	}
	else *OK = false;
	int* ss = new int[nss];
	for (j = 0; j < nss; j++) ss[j] = nuvar[j];
	if (not nuvar.empty()) nuvar.clear();
	*nssl = nss;
	//cout<<"locus "<<kloc<<"   nss="<<nss<<"\n";
	return ss;
}

long double ParticleC::cal_nss1p(int gr, int st) {
	int iloc, kloc = 0, nssl, nssm = 0, nl = 0, *ss, cat;
	long double res = 0.0;
	bool OK;
	auto&& samp = this->grouplist[gr].sumstat[st].samp.get();
	int sample = samp[0];
	for (iloc = 0; iloc < this->grouplist[gr].nloc; iloc++) {
		kloc = this->grouplist[gr].loc[iloc];
		cat = this->locuslist[kloc].type % 5;
		ss = this->cal_nsspl(kloc, sample, &nssl, &OK);
		//cout<<"nss1p   nssl="<<nssl<<"\n";
		//for (int k=0;k<nssl;k++) cout<<"  "<<ss[k];if(nssl>0) cout<<"\n";
		if (OK) {
			nl++;
			nssm += nssl;
		}
		delete []ss;
		if ((debuglevel == 9)and (kloc == 10)) {
			cout << "\ndnavar =" << this->locuslist[kloc].dnavar << "\n";
			for (int sa = 0; sa < this->nsample; sa++) {
				cout << " dans cal_nss1p   sample " << sa << "\n";
				for (int ind = 0; ind < dataobs.ssize[cat][sa]; ind++) cout << ind << "  " << this->locuslist[kloc].haplodna[sa][ind] << "\n";
			}
			cout << "nssl = " << nssl << "\n";
		}
	}
	if (nl > 0) res = (long double)nssm / (long double)nl;
	if ((debuglevel == 9)and (kloc == 10)) {
		cout << "nombre de sites segregeant = " << res << "\n";
		if (not this->dnatrue) exit(1);
	}
	return res;
}

long double ParticleC::cal_mpdpl(int kloc, int sample, int* nd) {
	int npdl = 0, di, k, ndd = 0, cat = this->locuslist[kloc].type % 5;
	long double res = 0.0;
	for (int i = 0; i < dataobs.ssize[cat][sample] - 1; i++) {
		for (int j = i + 1; j < dataobs.ssize[cat][sample]; j++) {
			if ((this->locuslist[kloc].haplodna[sample][i] != SEQMISSING)and (this->locuslist[kloc].haplodna[sample][j] != SEQMISSING)) {
				ndd++;
				di = 0;
				for (k = 0; k < this->locuslist[kloc].dnavar; k++) {
					if ((this->locuslist[kloc].haplodnavar[sample][i][k] != 'N')and (this->locuslist[kloc].haplodnavar[sample][j][k] != 'N')and (this->locuslist[kloc].haplodnavar[sample][i][k] != this->locuslist[kloc].haplodnavar[sample][j][k])) di++;
				}
				npdl += di;
			}
		}
	}
	if (ndd > 0) res = (long double)npdl / (long double)ndd;
	*nd = ndd;
	return res;
}

long double ParticleC::cal_mpd1p(int gr, int st) {
	int iloc, kloc, nd = 0, nl = 0;
	long double npdm = 0.0, res = 0.0, mpdp;
	auto&& samp = this->grouplist[gr].sumstat[st].samp.get();
	int sample = samp[0];
	for (iloc = 0; iloc < this->grouplist[gr].nloc; iloc++) {
		//npdl=0;
		kloc = this->grouplist[gr].loc[iloc];
		mpdp = this->cal_mpdpl(kloc, sample, &nd);
		if (nd > 0) {
			npdm += mpdp;
			nl++;
		}
	}
	if (nl > 0) res = (long double)npdm / (long double)nl;
	return res;

}

long double ParticleC::cal_vpd1p(int gr, int st) {
	int iloc, kloc, k, di, nl = 0, nd, cat;
	long double res = 0.0, mpd, spd, vpd = 0.0;
	auto&& samp = this->grouplist[gr].sumstat[st].samp.get();
	int sample = samp[0];
	for (iloc = 0; iloc < this->grouplist[gr].nloc; iloc++) {
		nd = 0;
		mpd = 0.0;
		spd = 0.0;
		kloc = this->grouplist[gr].loc[iloc];
		cat = this->locuslist[kloc].type % 5;
		for (int i = 0; i < dataobs.ssize[cat][sample] - 1; i++) {
			for (int j = i + 1; j < dataobs.ssize[cat][sample]; j++) {
				if ((this->locuslist[kloc].haplodna[sample][i] != SEQMISSING)and (this->locuslist[kloc].haplodna[sample][j] != SEQMISSING)) {
					nd++;
					di = 0;
					for (k = 0; k < this->locuslist[kloc].dnavar; k++) {
						if ((this->locuslist[kloc].haplodnavar[sample][i][k] != 'N')and (this->locuslist[kloc].haplodnavar[sample][j][k] != 'N')and (this->locuslist[kloc].haplodnavar[sample][i][k] != this->locuslist[kloc].haplodnavar[sample][j][k])) di++;
					}
					mpd += (long double)di;
					spd += (long double)di * (long double)di;
				}
			}
		}
		if (nd > 1) {
			spd = (spd - mpd * mpd / (long double)nd) / (long double)(nd - 1);
			nl++;
			vpd += spd;
		}
	}
	if (nl > 0) res = vpd / (long double)nl;
	return res;

}

long double ParticleC::cal_dta1pl(int kloc, int sample, bool* OKK) {
	long double a1, a2, b1, b2, c1, c2, e1, e2, S, pi;
	int nd, n = 0, *ss, kS, cat = this->locuslist[kloc].type % 5;
	bool OK;
	*OKK = true;
	if (this->locuslist[kloc].dnavar < 1) return 0.0;
	long double res = 0.0;
	for (int i = 0; i < dataobs.ssize[cat][sample]; i++) if (this->locuslist[kloc].haplodna[sample][i] != SEQMISSING) n++;
	if (n < 2) {
		*OKK = false;
		return 0.0;
	}
	a1 = 0;
	for (int i = 1; i < n; i++) a1 += 1.0 / (long double)i;
	a2 = 0;
	for (int i = 1; i < n; i++) a2 += 1.0 / (long double)(i * i);
	b1 = (long double)(n + 1) / (long double)(n - 1) / 3.0;
	b2 = 2.0 * ((long double)(n * n + n) + 3.0) / 9.0 / (long double)(n * n - n);
	c1 = b1 - 1.0 / a1;
	c2 = b2 - ((long double)(n + 2) / a1 / (long double)n) + (a2 / a1 / a1);
	e1 = c1 / a1;
	e2 = c2 / (a1 * a1 + a2);
	pi = cal_mpdpl(kloc, sample, &nd);
	ss = cal_nsspl(kloc, sample, &kS, &OK);
	S = (long double)kS;
	if ((nd > 0)and (OK)and (e1 * S + e2 * S * (S - 1.0) > 0.0)) res = (pi - S / a1) / sqrt(e1 * S + e2 * S * (S - 1.0));
	//cout<<"a1="<<a1<<"  a2="<<a2<<"b1="<<b1<<"  b2="<<b2<<"c1="<<c1<<"  c2="<<c2<<"e1="<<e1<<"  e2="<<e2<<"\n";
	//cout<<"nd="<<nd<<"   OK="<<OK<<"   pi="<<pi<<"   S="<<S<<"   res="<<res<<"\n";
	delete []ss;
	return res;
}

long double ParticleC::cal_dta1p(int gr, int st) {
	int iloc, kloc, nl = 0;
	long double res = 0.0, tal;
	bool OK;
	auto&& samp = this->grouplist[gr].sumstat[st].samp.get();
	int sample = samp[0];
	for (iloc = 0; iloc < this->grouplist[gr].nloc; iloc++) {
		kloc = this->grouplist[gr].loc[iloc];
		tal = this->cal_dta1pl(kloc, sample, &OK);
		if (OK) {
			nl++;
			res += tal;
		}
	}
	//cout<<"nl="<<nl<<"   res="<<res<<"\n";
	if (nl > 0) res /= nl;
	return res;
}

long double ParticleC::cal_pss1p(int gr, int st) {
	int iloc, kloc, nl = 0, nps = 0, **ssa, *nssa;
	bool trouve, OK;
	long double res = 0.0;
	auto&& samp = this->grouplist[gr].sumstat[st].samp.get();
	int sample = samp[0];
	for (iloc = 0; iloc < this->grouplist[gr].nloc; iloc++) {
		kloc = this->grouplist[gr].loc[iloc];
		ssa = new int*[this->nsample];
		nssa = new int[this->nsample];
		//on cherche les sites variables des différents échantillons
		for (int sa = 0; sa < this->nsample; sa++) ssa[sa] = this->cal_nsspl(kloc, sa, &nssa[sa], &OK);
		//cout<<"\nlocus "<<kloc<<"\n";
		for (int sa = 0; sa < this->nsample; sa++) {
			//cout<<"sample "<<sa+1<<"   nssa="<<nssa[sa]<<"\n";
			//for (int k=0;k<nssa[sa];k++) cout<<"   "<<ssa[sa][k];
			//cout<<"\n";
		}
		//on compte le nombre de sites variables de l'échantillon cible qui ne sont pas variables dans les autres échantillons
		nl++;
		for (int j = 0; j < nssa[sample]; j++) {
			trouve = false;
			for (int sa = 0; sa < this->nsample; sa++) {
				if ((sa != sample)and (nssa[sa] > 0)) {
					for (int k = 0; k < nssa[sa]; k++) {
						trouve = (ssa[sample][j] == ssa[sa][k]);
						if (trouve) break;
					}
				}
				if (trouve) break;
			}
			if (not trouve) nps++;
		}
		for (int sa = 0; sa < this->nsample; sa++) delete [] ssa[sa];
		delete []ssa;
		delete []nssa;
	}
	if (nl > 0) res = (long double)nps / (long double)nl;
	//cout<<"PSS_"<<this->grouplist[gr].sumstat[st].samp<<" = "<<res<<"   (nps="<<nps<<" nl="<<nl<<")\n";
	return res;
}

void ParticleC::afs(int sample, int iloc, int kloc) {
	int nf[4], ii, jj, cat;
	this->n_afs[sample][iloc] = 0;
	//cout<<"sample "<<sample<<"  locus "<<kloc<<"   dnavar = "<<this->locuslist[kloc].dnavar <<"\n";
	if (this->locuslist[kloc].dnavar == 0) this->n_afs[sample][iloc] = 0;
	else {
		for (int j = 0; j < this->locuslist[kloc].dnavar; j++) {
			//cout <<"site "<<j+1<<" sur "<<this->locuslist[kloc].dnavar<<"\n";
			for (int k = 0; k < 4; k++) nf[k] = 0;
			cat = this->locuslist[kloc].type % 5;
			for (int i = 0; i < dataobs.ssize[cat][sample]; i++) {
				//if (this->locuslist[kloc].haplodna[sample][i]!=SEQMISSING) {
				//cout<<this->locuslist[kloc].haplodnavar[sample][i][j]<<"\n";
				switch (this->locuslist[kloc].haplodnavar[sample][i][j]) {
					case 'A': nf[0]++;
						break;
					case 'C': nf[1]++;
						break;
					case 'G': nf[2]++;
						break;
					case 'T': nf[3]++;
						break;
				}
				//}
			}
			//tri des frequences dans l'ordre croissant
			for (int i = 0; i < 3; i++) {
				for (int k = i + 1; k < 4; k++)
					if (nf[i] > nf[k]) {
						ii = nf[i];
						nf[i] = nf[k];
						nf[k] = ii;
					}
			}
			//for (int i=0;i<4;i++) cout <<nf[i]<<"  ";cout<<"\n";
			//recherche de la plus petite fréquence non nulle
			jj = 0;
			while ((nf[jj] == 0)and (jj < 4)) jj++;
			if (jj < 3) {
				this->n_afs[sample][iloc]++;
				this->t_afs[sample][iloc].push_back(nf[jj]);
			}
		}
	}
	//cout<<"sample "<<sample<<"  locus "<<iloc<<"   n_afs = "<<this->n_afs[sample][iloc]<<"\n";
	//if (this->n_afs[sample][iloc]>0) {for (int i=0;i<n_afs[sample][iloc];i++) cout<<"  "<<t_afs[sample][iloc][i];cout<<"\n";}
	this->afsdone[sample][iloc] = true;
}

long double ParticleC::cal_mns1p(int gr, int st) {
	int iloc, kloc, nl = 0;
	long double res = 0.0;
	auto&& samp = this->grouplist[gr].sumstat[st].samp.get();
	int sample = samp[0];
	for (iloc = 0; iloc < this->grouplist[gr].nloc; iloc++) {
		kloc = this->grouplist[gr].loc[iloc];
		if (not this->afsdone[sample][iloc]) this->afs(sample, iloc, kloc);
		for (int i = 0; i < this->n_afs[sample][iloc]; i++) res += (long double)this->t_afs[sample][iloc][i] / (long double)this->n_afs[sample][iloc];
		nl++;
	}
	if (nl > 0) res /= (long double)nl;
	return res;
}

long double ParticleC::cal_vns1p(int gr, int st) {
	int iloc, kloc, nl = 0;
	long double res = 0.0, sx, sx2, a, v;
	auto&& samp = this->grouplist[gr].sumstat[st].samp.get();
	int sample = samp[0];
	for (iloc = 0; iloc < this->grouplist[gr].nloc; iloc++) {
		kloc = this->grouplist[gr].loc[iloc];
		//cout<<"Locus "<<kloc<<"\n";
		if (not afsdone[sample][iloc]) afs(sample, iloc, kloc);
		//cout<<"apres afs\n";
		if (n_afs[sample][iloc] > 1) {
			sx = 0.0;
			sx2 = 0.0;
			for (int i = 0; i < n_afs[sample][iloc]; i++) {
				a = (long double)t_afs[sample][iloc][i];
				sx += a;
				sx2 += a * a;
			}
			//cout<<"sample "<<sample<<"   locus "<<kloc<<"  sx="<<sx<<"  sx2="<<sx2<<"\n";
			a = (long double)n_afs[sample][iloc];
			v = (sx2 - sx * sx / a) / a;
		}
		else {
			v = 0.0;/*cout<<"vns="<<v<<"\n";*/
		}
		nl++;//cout<<"nl="<<nl<<"\n";
		//if (v>0.0) cout<<"sample "<<sample<<"   locus "<<kloc<<"  sx="<<sx<<"  sx2="<<sx2<<"  n="<<a<<"   vns="<<v<<"\n";
		if (v > 0.0) res += v;
	}
	if (nl > 0) res /= (long double)nl;
	//exit(1);
	return res;
}

long double ParticleC::cal_nha2p(int gr, int st) {
	string* haplo;
	int iloc, kloc, j, nhl = 0, nhm = 0, nl = 0, sample, dmax, cat;
	bool trouve;
	long double res = 0.0;
	auto&& samp = this->grouplist[gr].sumstat[st].samp.get();
	int samp0 = samp[0];
	int samp1 = samp[1];
	//cout <<"samples "<<samp0<<" & "<<samp1<<"\n";
	for (iloc = 0; iloc < this->grouplist[gr].nloc; iloc++) {
		nhl = 0;
		kloc = this->grouplist[gr].loc[iloc];
		cat = this->locuslist[kloc].type % 5;
		dmax = dataobs.ssize[cat][samp0] + dataobs.ssize[cat][samp1];
		haplo = new string[dmax];
		if (samplesize(kloc, samp0) + samplesize(kloc, samp1) > 0) {
			nl++;
			if (this->locuslist[kloc].dnavar == 0) {
				haplo[nhl] = "";
				nhl++;
			}
			else {
				for (int samp = 0; samp < 2; samp++) {
					if (samp == 0) sample = samp0;
					else sample = samp1;
					for (int i = 0; i < dataobs.ssize[cat][sample]; i++) {
						if (this->locuslist[kloc].haplodna[sample][i] != SEQMISSING) {
							if (nhl == 0) {
								haplo[nhl] = this->locuslist[kloc].haplodnavar[sample][i];
								nhl++;
							}
							else {
								trouve = false;
								j = 0;
								while ((not trouve)and (j < nhl)) {
									trouve = identseq(haplo[j], this->locuslist[kloc].haplodnavar[sample][i], "nah2p", kloc, sample, i);
									if (not trouve) j++;
								}
								if (trouve) {
									for (int k = 0; k < this->locuslist[kloc].dnavar; k++) if (haplo[j][k] == 'N') haplo[j][k] = this->locuslist[kloc].haplodnavar[sample][i][k];
								}
								else {
									haplo[nhl] = this->locuslist[kloc].haplodnavar[sample][i];
									nhl++;
								}
							}
						}
					}
				}
			}
		}
		nhm += nhl;
		delete [] haplo;
	}
	if (nl > 0) res = (long double)nhm / (long double)nl;
	return res;
}

int* ParticleC::cal_nss2pl(int kloc, int samp0, int samp1, int* nssl, bool* OK) {
	char c0;
	bool ident;
	int j, *ss, nss = 0, cat = this->locuslist[kloc].type % 5;
	vector<int> nuvar;
	if (samplesize(kloc, samp0) + samplesize(kloc, samp1) > 0) {
		*OK = true;
		if (this->locuslist[kloc].dnavar != 0) {
			for (j = 0; j < this->locuslist[kloc].dnavar; j++) {
				c0 = '\0';
				ident = true;
				for (int i = 0; i < dataobs.ssize[cat][samp0]; i++) {
					if ((this->locuslist[kloc].haplodna[samp0][i] != SEQMISSING)and (this->locuslist[kloc].haplodnavar[samp0][i][j] != 'N')) {
						if (c0 == '\0') c0 = this->locuslist[kloc].haplodnavar[samp0][i][j];
						else ident = (c0 == this->locuslist[kloc].haplodnavar[samp0][i][j]);
						//if (not ident) cout<<"nucleotide "<<j<<"   c0="<<c0<<"   c1="<<this->locuslist[kloc].haplodnavar[sample][i][j]<<"   i="<<i<<"\n";
					}
					if (not ident) break;
				}
				if (ident) {
					for (int i = 0; i < dataobs.ssize[cat][samp1]; i++) {
						if ((this->locuslist[kloc].haplodna[samp1][i] != SEQMISSING)and (this->locuslist[kloc].haplodnavar[samp1][i][j] != 'N')) {
							if (c0 == '\0') c0 = this->locuslist[kloc].haplodnavar[samp1][i][j];
							else ident = (c0 == this->locuslist[kloc].haplodnavar[samp1][i][j]);
							//if (not ident) cout<<"nucleotide "<<j<<"   c0="<<c0<<"   c1="<<this->locuslist[kloc].haplodnavar[sample][i][j]<<"   i="<<i<<"\n";
						}
						if (not ident) break;
					}
				}
				if (not ident) {
					nss++;
					nuvar.push_back(j);
				}
			}
		}
	}
	else *OK = false;
	ss = new int[nss];
	for (j = 0; j < nss; j++) ss[j] = nuvar[j];
	*nssl = nss;
	return ss;
}

long double ParticleC::cal_nss2p(int gr, int st) {
	int iloc, kloc, nssl, nssm = 0, nl = 0, *ss;
	ss = NULL;
	long double res = 0.0;
	bool OK;
	//cout<<"\n";
	auto&& samp = this->grouplist[gr].sumstat[st].samp.get();
	int samp0 = samp[0];
	int samp1 = samp[1];
	for (iloc = 0; iloc < this->grouplist[gr].nloc; iloc++) {
		kloc = this->grouplist[gr].loc[iloc];
		if (ss != NULL) {
			delete [] ss;
			ss = NULL;
		}
		ss = this->cal_nss2pl(kloc, samp0, samp1, &nssl, &OK);
		//cout<<"nssl = "<<nssl<<"\n";
		//for (int k=0;k<nssl;k++) cout<<"  "<<ss[k];if(nssl>0) cout<<"\n";
		if (OK) {
			nl++;
			nssm += nssl;
		}
		//cout<<"cumul nl="<<nl<<"   cumul nssm="<<nssm<<"\n";
		delete []ss;
		ss = NULL;
	}
	if (nl > 0) res = (long double)nssm / (long double)nl;
	return res;

}

long double ParticleC::cal_mpw2pl(int kloc, int samp0, int samp1, bool* OK) {
	int isamp, samp, di, mdl = 0, nd = 0, cat = this->locuslist[kloc].type % 5;
	if (samplesize(kloc, samp0) + samplesize(kloc, samp1) > 0) {
		*OK = true;
		if (this->locuslist[kloc].dnavar != 0) {
			for (isamp = 0; isamp < 2; isamp++) {
				if (isamp == 0) samp = samp0;
				else samp = samp1;
				//cout<<" sample "<<samp<<"   ss="<< this->locuslist[kloc].ss[samp]<<"\n";
				for (int i = 0; i < dataobs.ssize[cat][samp] - 1; i++) {
					if (this->locuslist[kloc].haplodna[samp][i] != SEQMISSING) {
						//cout <<"coucou i\n";
						for (int j = i + 1; j < dataobs.ssize[cat][samp]; j++) {
							if (this->locuslist[kloc].haplodna[samp][j] != SEQMISSING) {
								//cout<<"coucou j\n";
								nd++;
								di = 0;
								for (int k = 0; k < this->locuslist[kloc].dnavar; k++) {
									if ((this->locuslist[kloc].haplodnavar[samp][i][k] != 'N')and (this->locuslist[kloc].haplodnavar[samp][j][k] != 'N')and (this->locuslist[kloc].haplodnavar[samp][i][k] != this->locuslist[kloc].haplodnavar[samp][j][k])) di++;
								}
								mdl += di;
							}
							//cout<<"individus "<<i<<" et "<<j<<"   di="<<di<<"  nd="<<nd<<"   mdl="<<mdl<<"\n";
						}
					} //else cout <<"MISSING  ";
				}
			}
		}
	}
	else *OK = false;
	if (nd > 0) return (long double)mdl / (long double)nd;
	else return (long double)mdl;

}


long double ParticleC::cal_mpw2p(int gr, int st) {
	int iloc, kloc, nl = 0;
	long double res = 0.0, md = 0.0, mdl;
	bool OK;
	auto&& samp = this->grouplist[gr].sumstat[st].samp.get();
	int samp0 = samp[0];
	int samp1 = samp[1];
	for (iloc = 0; iloc < this->grouplist[gr].nloc; iloc++) {
		kloc = this->grouplist[gr].loc[iloc];
		mdl = cal_mpw2pl(kloc, samp0, samp1, &OK);
		if (OK) {
			nl++;
			md += mdl;
		}
	}
	if (nl > 0) res = md / (long double)nl;
	return res;
}

long double ParticleC::cal_mpb2pl(int kloc, int samp0, int samp1, bool* OK) {
	int di, mdl = 0, nd = 0, cat = this->locuslist[kloc].type % 5;
	if (samplesize(kloc, samp0) + samplesize(kloc, samp1)) {
		*OK = true;
		if (this->locuslist[kloc].dnavar != 0) {
			for (int i = 0; i < dataobs.ssize[cat][samp0]; i++) {
				if (this->locuslist[kloc].haplodna[samp0][i] != SEQMISSING) {
					//cout <<"coucou i\n";
					for (int j = 0; j < dataobs.ssize[cat][samp1]; j++) {
						if (this->locuslist[kloc].haplodna[samp1][j] != SEQMISSING) {
							//cout<<"coucou j\n";
							nd++;
							di = 0;
							for (int k = 0; k < this->locuslist[kloc].dnavar; k++) {
								if ((this->locuslist[kloc].haplodnavar[samp0][i][k] != 'N')and (this->locuslist[kloc].haplodnavar[samp1][j][k] != 'N')and (this->locuslist[kloc].haplodnavar[samp0][i][k] != this->locuslist[kloc].haplodnavar[samp1][j][k])) di++;
							}
							mdl += di;
						}
						//cout<<"individus "<<i<<" et "<<j<<"   di="<<di<<"  nd="<<nd<<"   mdl="<<mdl<<"\n";
					}
				} //else cout <<"MISSING  ";
			}
		}
	}
	else *OK = false;
	if (nd > 0) return (long double)mdl / (long double)nd;
	else return (long double)mdl;
}

long double ParticleC::cal_mpb2p(int gr, int st) {
	int iloc, kloc, nl = 0;
	long double res = 0.0, md = 0.0, mdl;
	bool OK;
	auto&& samp = this->grouplist[gr].sumstat[st].samp.get();
	int samp0 = samp[0];
	int samp1 = samp[1];
	for (iloc = 0; iloc < this->grouplist[gr].nloc; iloc++) {
		kloc = this->grouplist[gr].loc[iloc];
		mdl = cal_mpb2pl(kloc, samp0, samp1, &OK);
		if (OK) {
			nl++;
			md += mdl;
		}
	}
	if (nl > 0) res = md / (long double)nl;
	return res;

}

long double ParticleC::cal_fst2p(int gr, int st) {
	int iloc, kloc, nl = 0;
	long double res = 0.0, num = 0.0, den = 0.0, Hw, Hb;
	bool OKw, OKb;
	auto&& samp = this->grouplist[gr].sumstat[st].samp.get();
	int samp0 = samp[0];
	int samp1 = samp[1];
	for (iloc = 0; iloc < this->grouplist[gr].nloc; iloc++) {
		kloc = this->grouplist[gr].loc[iloc];
		Hw = cal_mpw2pl(kloc, samp0, samp1, &OKw);
		Hb = cal_mpb2pl(kloc, samp0, samp1, &OKb);
		if ((OKw)and (OKb)) {
			nl++;
			num += Hb - Hw;
			den += Hb;
		}
	}
	if (den > 0) res = num / den;
	return res;

}

void ParticleC::cal_freq(int gr, int st) {
	int iloc, kloc, nhaplo, isamp, sample, dmax, j, cat;
	long double d;
	auto&& samp = this->grouplist[gr].sumstat[st].samp.get();
	int samp0 = samp[0];
	int samp1 = samp[1];
	int samp2 = samp[2];
	string* haplo;
	bool trouve;
	for (iloc = 0; iloc < this->grouplist[gr].nloc; iloc++) {
		kloc = this->grouplist[gr].loc[iloc];
		cat = this->locuslist[kloc].type % 5;
		dmax = samplesize(kloc, samp0) + samplesize(kloc, samp1) + samplesize(kloc, samp2);
		//cout<<"Locus "<<kloc<<"   (cal_freq)\n";
		this->locuslist[kloc].freq = vector<vector<long double>>(3);
		this->locuslist[kloc].haplomic = vector<vector<int>>(3);
		this->locuslist[kloc].kmin = 100;
		nhaplo = 0;
		if (this->locuslist[kloc].dnavar == 0) { // construction pour ne pas deleter du vide (haplomic)
			for (isamp = 0; isamp < 3; isamp++) {
				if (isamp == 0) sample = samp0;
				else if (isamp == 1) sample = samp1;
				else sample = samp2;
				this->locuslist[kloc].freq[sample] = vector<long double>(1);
				this->locuslist[kloc].freq[sample][0] = 1.0;
				this->locuslist[kloc].haplomic[sample] = vector<int>(dataobs.ssize[cat][sample]);
				for (int j = 0; j < dataobs.ssize[cat][sample]; j++)
					this->locuslist[kloc].haplomic[sample][j] = this->locuslist[kloc].kmin;
			}
			//haplo[nhaplo] = "";
			nhaplo++;
		}
		else {
			haplo = new string[dmax];
			for (isamp = 0; isamp < 3; isamp++) {
				if (isamp == 0) sample = samp0;
				else if (isamp == 1) sample = samp1;
				else sample = samp2;
				for (int i = 0; i < dataobs.ssize[cat][sample]; i++) {
					if (this->locuslist[kloc].haplodna[sample][i] != SEQMISSING) {
						if (nhaplo == 0) {
							haplo[nhaplo] = this->locuslist[kloc].haplodnavar[sample][i];
							nhaplo++;
						}
						else {
							trouve = false;
							j = 0;
							while ((not trouve)and (j < nhaplo)) {
								trouve = identseq(haplo[j], this->locuslist[kloc].haplodnavar[sample][i], "cal_freq1", kloc, sample, i);
								if (not trouve) j++;
							}
							if (trouve) {
								for (int k = 0; k < this->locuslist[kloc].dnavar; k++) if (haplo[j][k] == 'N') haplo[nhaplo][k] = this->locuslist[kloc].haplodnavar[sample][i][k];
							}
							else {
								haplo[nhaplo] = this->locuslist[kloc].haplodnavar[sample][i];
								nhaplo++;
							}
						}
					}
				}
			}
			for (isamp = 0; isamp < 3; isamp++) {
				if (isamp == 0) sample = samp0;
				else if (isamp == 1) sample = samp1;
				else sample = samp2;
				this->locuslist[kloc].freq[sample] = vector<long double>(nhaplo);
				for (j = 0; j < nhaplo; j++) this->locuslist[kloc].freq[sample][j] = 0.0;
				d = 1.0 / (long double)samplesize(kloc, sample);
				this->locuslist[kloc].haplomic[sample] = vector<int>(dataobs.ssize[cat][sample]);
				for (int i = 0; i < dataobs.ssize[cat][sample]; i++) {
					if (this->locuslist[kloc].haplodna[sample][i] != SEQMISSING) {
						trouve = false;
						j = 0;
						while ((not trouve)and (j < nhaplo)) {
							trouve = identseq(haplo[j], this->locuslist[kloc].haplodnavar[sample][i], "cal_freq2", kloc, sample, i);
							if (not trouve) j++;
						}
						this->locuslist[kloc].freq[sample][j] += d;
						this->locuslist[kloc].haplomic[sample][i] = this->locuslist[kloc].kmin + j;
					}
					else this->locuslist[kloc].haplomic[sample][i] = MICMISSING;
				}
				//for (j=0;j<nhaplo;j++) cout<<"freq[loc]["<<j<<"]="<<this->locuslist[kloc].freq[sample][j]<<"\n";
				//cout<<"\n";
			}
			delete []haplo;
			//cout<<"cal_freq  nhaplo="<<nhaplo<<"\n";

		}
	}
}

long double ParticleC::cal_aml3p(int gr, int st) {
	int iloc, kloc, nlocutil = 0, *ss = NULL, nssl;
	long double p1, p2, p3, lik1, lik2, lik3;
	int i1 = 1, i2 = 998, i3;
	bool OK;
	long double res = 0.0;
	pair<long double, long double> c;
	//cout<<"\ndebut de cal_aml3p\n";
	auto&& samp = this->grouplist[gr].sumstat[st].samp.get();
	int samp1 = samp[1];
	int samp2 = samp[2];
	//cout<<"DEBUT D AML3P\n";


	for (iloc = 0; iloc < this->grouplist[gr].nloc; iloc++) {
		kloc = this->grouplist[gr].loc[iloc];
		ss = this->cal_nss2pl(kloc, samp1, samp2, &nssl, &OK);
		if ((OK)and (nssl > 0)) nlocutil++;
		delete [] ss;
	}
	if (nlocutil < 1) return 0.5;
	//cout<<"avant cal_freq\n";
	cal_freq(gr, st);
	c = pente_lik(gr, st, i1);
	lik1 = c.first;
	p1 = c.second;
	c = pente_lik(gr, st, i2);
	lik2 = c.first;
	p2 = c.second;
	//cout<<"apres les pente_lik  p1="<<p1<<"   p2="<<p2<<"   i1="<<i1<<"   i2="<<i2<<"\n";
	//cout<<"lik1="<<lik1<<"   lik2="<<lik2<<"\n";
	if (abs(lik1) + abs(lik2) < 1.0E-10) {
		libere_freq(gr);
		return -9999.0;
	}
	if ((p1 < 0.0)and (p2 < 0.0)) {
		libere_freq(gr);
		return 0.0;
	}
	if ((p1 > 0.0)and (p2 > 0.0)) {
		libere_freq(gr);
		return 1.0;
	}
	do {
		i3 = (i1 + i2) / 2;
		c = pente_lik(gr, st, i3);
		lik3 = c.first;
		p3 = c.second;
		//cout<<"p3="<<p3<<"   i3="<<i3<<"   lik3="<<lik3<<"\n";
		if (p1 * p3 < 0.0) {
			i2 = i3;
			p2 = p3;
			lik2 = lik3;
		}
		else {
			i1 = i3;
			p1 = p3;
			lik1 = lik3;
		}
	}
	while (abs(i2 - i1) > 1);
	if (lik1 > lik2) res = 0.001 * (long double)i1;
	else res = 0.001 * (long double)i2;
	//cout<<"fin de cal_aml3p res="<<res<<"\n";
	libere_freq(gr);
	//cout<<"apreslibere_freq\n";
	return res;
}

void ParticleC::cal_numvar(int gr) {
	int i, j, k, ns, pop, j0, cat = 0;
	int locus, nlocs = this->grouplist[gr].nloc;
	string ss;
	char* site;
	bool ident;
	//this->nvar = new int[nlocs];
	//this->numvar = new int*[nlocs];
	vector<int> nuvar;
	//cout <<"debut de cal_numvar dnatrue ="<<this->dnatrue<<"\n";
	if (not this->dnatrue) {
		for (int iloc = 0; iloc < nlocs; iloc++) {
			locus = this->grouplist[gr].loc[iloc];
			cat = this->locuslist[locus].type % 5;
			this->locuslist[locus].dnavar = 0;
			for (pop = 0; pop < this->nsample; pop++) {
				for (i = 0; i < dataobs.ssize[cat][pop]; i++) {
					if ((int)this->locuslist[locus].haplodna[pop][i].length() > this->locuslist[locus].dnavar) {
						this->locuslist[locus].dnavar = this->locuslist[locus].haplodna[pop][i].length();
						//cout<<"  >"<<this->locuslist[locus].haplodna[pop][i]<<"<  ("<<this->locuslist[locus].haplodna[pop][i].length()<<")\n";
					}
				}
				if (this->locuslist[locus].dnavar < 0) {
					//cout<<"Locus "<<locus<<"   sample "<<pop <<"   dnavar="<<this->locuslist[locus].dnavar<<"\n";
					//this->locuslist[locus].dnavar = 0;
					exit(1);
				}
			}
			//cout<<"coucou\n";
			this->locuslist[locus].haplodnavar = vector<vector<string>>(this->nsample);
			for (pop = 0; pop < this->nsample; pop++) {
				this->locuslist[locus].haplodnavar[pop] = vector<string>(dataobs.ssize[cat][pop]);
				for (i = 0; i < dataobs.ssize[cat][pop]; i++) {
					if ((int)this->locuslist[locus].haplodna[pop][i].length() == this->locuslist[locus].dnavar) {
						this->locuslist[locus].haplodnavar[pop][i] = this->locuslist[locus].haplodna[pop][i];
					}
					else {
						this->locuslist[locus].haplodnavar[pop][i] = "";
						for (int j = 0; j < this->locuslist[locus].dnavar; j++)this->locuslist[locus].haplodnavar[pop][i] += "N";
					}
				}
			}
		}
	}
	else {
		for (int iloc = 0; iloc < nlocs; iloc++) {
			locus = this->grouplist[gr].loc[iloc];
			cat = this->locuslist[locus].type % 5;
			this->locuslist[locus].haplodnavar = vector<vector<string>>(this->nsample);
			ns = 0;
			for (pop = 0; pop < this->nsample; pop++) ns += dataobs.ssize[cat][pop];
			site = new char[ns];
			this->locuslist[locus].dnavar = 0;
			if (this->locuslist[locus].dnalength > 0) {
				//cout <<"longseq = "<<this->locuslist[locus].dnalength<<"\n";;
				for (k = 0; k < this->locuslist[locus].dnalength; k++) {
					j = -1;
					for (pop = 0; pop < this->nsample; pop++) {
						for (i = 0; i < dataobs.ssize[cat][pop]; i++) {
							j++;
							if (this->locuslist[locus].haplodna[pop][i] == SEQMISSING) site[j] = 'N';
							else site[j] = this->locuslist[locus].haplodna[pop][i][k];
							//if ((locus==13)and(j<1000)) cout<<"locus"<<locus<<"   site["<<j<<"] = "<<site[j]<<"\n";
						}
					}
					//cout<<"fin de la boucle pop\n";
					j0 = 0;
					while ((site[j0] == 'N')and (j0 < ns - 1)) j0++; //recherche du premier nucléotide non N
					j = j0 + 1;
					ident = true;
					while ((ident)and (j < ns)) {
						ident = ((site[j] == 'N')or (site[j] == site[j0]));
						j++;
					}
					if (not ident) {
						this->locuslist[locus].dnavar++;
						nuvar.push_back(k);
					}
				}
			}
			//cout<<"locus "<<locus<<"   dnavar = "<<this->locuslist[locus].dnavar<<"\n";
			for (pop = 0; pop < this->nsample; pop++) {
				this->locuslist[locus].haplodnavar[pop] = vector<string>(dataobs.ssize[cat][pop]);
				for (i = 0; i < dataobs.ssize[cat][pop]; i++) {
					ss = "";
					for (int k = 0; k < this->locuslist[locus].dnavar; k++) {
						if (this->locuslist[locus].haplodna[pop][i] != SEQMISSING) {
							ss += this->locuslist[locus].haplodna[pop][i][nuvar[k]];
						}
						else {
							ss += "N";
						}
					}
					this->locuslist[locus].haplodnavar[pop][i] = ss;
					//for(int k=0;k<this->locuslist[locus].dnavar;k++) cout<<"   " <<nuvar[k]; cout<<"\n";
				}
			}
			if (not nuvar.empty()) nuvar.clear();
			delete []site;
			ss.clear();
		}
	}
}

/*  Numérotation des stat
 *  1 : nal			-1 : nha	    -13 : fst
 *  2 : het			-2 : nss            -14 : aml
 *  3 : var			-3 : mpd
 *  4 : MGW			-4 : vpd
 *  5 : Fst			-5 : dta
 *  6 : lik			-6 : pss
 *  7 : dm2			-7 : mns
 *  8 : n2P			-8 : vns
 *  9 : h2P			-9 : nh2
 * 10 : v2P		       -10 : ns2
 * 11 : das                    -11 : mp2
 * 12 : Aml                    -12 : mpb
 *
 *
 */

 ////////////// Sous-stats (mean/variance)

void ParticleC::cal_snpstatRedinit(int gr, int numsnp)
{
	StatsnpC& stsnp = this->grouplist[gr].sumstatsnp[numsnp];
	stsnp.sw = 0.0;
	stsnp.mx = 0.0;
	stsnp.sw0 = 0.0;
	stsnp.sw1 = 0.0;
	stsnp.sw12 = 0.0;
	stsnp.mx1 = 0.0;
	stsnp.mx2 = 0.0;
	stsnp.sw2 = 0.0;
	stsnp.mx12 = 0.0;
	stsnp.densw0 = 0.0;
	stsnp.densw1 = 0.0;
	stsnp.densw12 = 0.0;
	stsnp.denmx = 0.0;
	stsnp.denmx1 = 0.0;
	stsnp.denmx12 = 0.0;
}

void ParticleC::cal_snpstatRedacc(int gr, int numsnp, long double x , long double w) {
	StatsnpC& stsnp = this->grouplist[gr].sumstatsnp[numsnp];
	long double mo;
	stsnp.sw += w;
	if (x < 0.000000001) stsnp.sw0 += w;
	if (w > 0.0) {
		mo = stsnp.mx;
		stsnp.mx += (w / stsnp.sw) * (x - mo);
		stsnp.mx2 += w *(x - mo) * (x - stsnp.mx);
		stsnp.sw2 += w * w;
		if (x > 0.0) {
			stsnp.sw1 += w;
			stsnp.sw12 += w * w;
			mo = stsnp.mx1;
			stsnp.mx1 += (w / stsnp.sw1) * (x - mo);
			stsnp.mx12 +=  w * (x - mo) * (x - stsnp.mx1);
		}
	}
	
}
	
long double ParticleC::cal_p0L(int gr, int numsnp) {
	StatsnpC& stsnp = this->grouplist[gr].sumstatsnp[numsnp];	
	if (stsnp.sw < 0.000000001) return 0.0;
	return stsnp.sw0 / stsnp.sw;
}

//  long double ParticleC::cal_p0L(StatsnpC stsnp) {
// 	long double p0 = 0.0;
// 	if (stsnp.sw < 0.000000001) return 0.0;
// 	for (int i = 0; i < stsnp.n; i++) {
// 		//if ((i==0)or(i==stsnp.n-1)) cout<<"dans cal_p0L  stsnp.x["<<i<<"]="<<stsnp.x[i]<<"      stsnp.w["<<i<<"]="<<stsnp.w[i]<<"\n";
// 		//cout<<"dans cal_p0L  					stsnp.w["<<i<<"]="<<stsnp.w[i]<<"\n";
// 		if (stsnp.x[i] < 0.000000001) p0 += stsnp.w[i];
// 	}
// 	//exit(3);
// 	stsnp.w.clear();
// 	stsnp.x.clear();
// 	return p0 / stsnp.sw;
// }

long double ParticleC::cal_moyL0(int gr, int numsnp) {
	StatsnpC& stsnp = this->grouplist[gr].sumstatsnp[numsnp];
	if (stsnp.sw1 < 0.000000001) return 0.0;	
	return stsnp.mx1;
}


// long double ParticleC::cal_moyL0(StatsnpC stsnp) {
// 	long double sx = 0.0, sw = 0.0;
// 	for (int i = 0; i < stsnp.n; i++) {
// 		if (stsnp.w[i] > 0.0) {
// 			if (stsnp.x[i] > 0.0) {
// 				sx += stsnp.w[i] * stsnp.x[i];
// 				sw += stsnp.w[i];
// 			}
// 		}
// 	}
// 	stsnp.w.clear();
// 	stsnp.x.clear();
// 	if (sw > 0.0) return sx / sw;
// 	else return 0.0;
// }


long double ParticleC::cal_varL0(int gr, int numsnp) {
	StatsnpC& stsnp = grouplist[gr].sumstatsnp[numsnp];
	long double sw2diff = stsnp.sw1 * stsnp.sw1 - stsnp.sw12;
	if (sw2diff < 0.000000001) return 0.0;	
	return ( stsnp.mx12 * stsnp.sw1) / sw2diff;
}

long double ParticleC::cal_varL(int gr, int numsnp) {
	StatsnpC& stsnp = grouplist[gr].sumstatsnp[numsnp];
	long double sw2diff = stsnp.sw * stsnp.sw - stsnp.sw2;
	if (sw2diff < 0.000000001) return 0.0;	
	return ( stsnp.mx2 * stsnp.sw) / sw2diff;
}

// long double ParticleC::cal_varL0(StatsnpC stsnp) {
// 	long double sx = 0.0, sx2 = 0.0, sw = 0, sw2 = 0.0, wi;
// 	for (int i = 0; i < stsnp.n; i++) if (stsnp.x[i] > 0.0) sw += stsnp.w[i];
// 	if (sw == 0.0) return 0.0;
// 	for (int i = 0; i < stsnp.n; i++) {
// 		if (stsnp.w[i] > 0) {
// 			if (stsnp.x[i] > 0.0) {
// 				wi = stsnp.w[i] / sw;
// 				sx += wi * stsnp.x[i];
// 				sx2 += wi * (stsnp.x[i] * stsnp.x[i]);
// 				sw2 += wi * wi;
// 			}
// 		}
// 	}
// 	stsnp.w.clear();
// 	stsnp.x.clear();
// 	if (sw2 < 1.0) return (sx2 - sx * sx) / (1.0 - sw2);
// 	else return 0.0;
// }

long double ParticleC::cal_moyL(int gr, int numsnp) {
	StatsnpC& stsnp = this->grouplist[gr].sumstatsnp[numsnp];
	if (stsnp.sw < 0.000000001)	return 0.0;
	return stsnp.mx;
}


// long double ParticleC::cal_moyL(StatsnpC stsnp) {
// 	long double sx = 0.0, sw = 0.0;
// 	for (int i = 0; i < stsnp.n; i++) {
// 		if (stsnp.w[i] > 0) {
// 			sx += stsnp.w[i] * stsnp.x[i];
// 			sw += stsnp.w[i];
// 		}
// 	}
// 	if (sw > 0) return sx / sw;
// 	else return 0.0;
// }

////////////// Fin Meta-stats (mean/variance)

long double ParticleC::mQ(int n0, int n1, long double* x) {
	long double sx = 0;
	for (int i = n0; i < n1; i++) sx += x[i];
	//cout <<"x["<<n0<<"]="<<x[n0]<<"    x["<<n1-1<<"]="<<x[n1-1]<<"       ";
	if (n1 > n0) return sx / (long double)(n1 - n0);
	else return 0.0;
}

long double ParticleC::cal_medL(int n, long double* x) {
	long double sx = 0.0, nn = 0.0;
	for (int i = 0; i < n; i++)
		if ((x[i] != 0.0)and (x[i] != 1.0)) {
			sx += x[i];
			nn += 1.0;
		}
	if (nn > 0.0) return sx / nn;
	else return 0.0;
	/*int m=0,nn;
	  while(x[m]==0.0) m++;
	  nn=n-m;
	  //cout<<"ME    n="<<n<<"   n="<<m<<"   nn="<<nn<<"      ";
	  return mQ(m+nn/3,m+2*nn/3,x);*/
	/*if (nn==0) return 0.0;
	  if (nn==1) return x[n-1];
	  if (nn==2) return 0.5*(x[n-2]+x[n-1]);
	  if ((nn%2)==0) return 0.5*(x[m+nn/2-1]+x[m+nn/2]);
	  else return x[m+nn/2];*/
}

long double ParticleC::cal_qu1L(int n, long double* x) {
	long double p0 = 0.0, nn = (long double)n;
	for (int i = 0; i < n; i++) if (x[i] == 1.0) p0 += 1.0;
	return p0 / nn;
	/*int m=0,nn;
	  while(x[m]==0.0) m++;
	  nn=n-m;
	  //cout<<"Q1    n="<<n<<"   n="<<m<<"   nn="<<nn<<"      ";
	  return mQ(m,m+nn/3,x);*/
	/*if (nn<4) return 0.0;
	  if ((nn % 4)==0) return 0.5*(x[m+(nn/4)]+x[m+(nn/4)-1]);
	  else if ((nn % 2)==0) return x[m+(nn/4)];
	  else return (0.25*x[m+(nn/4)-1])+(0.75*x[m+(nn/4)]);*/
}

long double ParticleC::cal_qu3L(int n, long double* x) {
	int m = 0, nn;
	//long double a;
	while (x[m] == 0.0) m++;
	nn = n - m;
	//cout<<"Q3    n="<<n<<"   n="<<m<<"   nn="<<nn<<"      ";
	return mQ(m + 2 * nn / 3, n, x);
	/*if (nn<4) a=0.0;
	  else if ((nn % 4)==0) a=0.5*(x[m+(3*nn/4)]+x[m+(3*nn/4)-1]);
	  else if ((nn % 2)==0) a=x[m+(3*nn/4)];
	  else a= (0.25*x[m+(3*nn/4)+1])+(0.75*x[m+(3*nn/4)]);
	  if(a!=a) {cout<<"nn="<<nn<<"   m="<<m<<"\n";exit(3);}
	  return a;*/
}

void ParticleC::docalstat(int gr) {
	//cout << this->grouplist[gr].type<<"\n";
	if (this->grouplist[gr].type == 0) this->calfreq(gr);
	else if (this->grouplist[gr].type == 1) this->cal_numvar(gr);
	else if (this->grouplist[gr].type == 2) this->calfreqsnp(gr);
	else if (this->grouplist[gr].type == 3) this->calfreqpool(gr);
	//cout << "apres calfreq\n";
	for (int st = 0; st < this->grouplist[gr].nstat; st++) {
		if ((this->grouplist[gr].sumstat[st].cat == -7)or (this->grouplist[gr].sumstat[st].cat == -8)) {
			this->afsdone.resize(this->nsample);
			this->t_afs.resize(this->nsample);
			this->n_afs.resize(this->nsample);
			for (int sa = 0; sa < this->nsample; sa++) {
				this->afsdone[sa].resize(this->grouplist[gr].nloc);
				this->t_afs[sa].resize(this->grouplist[gr].nloc);
				this->n_afs[sa].resize(this->grouplist[gr].nloc);
				for (int loc = 0; loc < this->grouplist[gr].nloc; loc++) this->afsdone[sa][loc] = false;
			}
		}
	}
	int numsnp;
	//cout<<"Docalstat groupe "<<gr<<"    nstat="<<this->grouplist[gr].nstat<<"\n";
	for (int st = 0; st < this->grouplist[gr].nstat; st++) {
		/*	
		if (this->grouplist[gr].sumstat[st].cat<5)
		{cout <<" calcul de la stat "<<st<<"   cat="<<this->grouplist[gr].sumstat[st].cat<<"   group="<<gr<<"   samp = "<<this->grouplist[gr].sumstat[st].samp  <<"\n";fflush(stdin);}
		else if (this->grouplist[gr].sumstat[st].cat<12)
		{cout <<" calcul de la stat "<<st<<"   cat="<<this->grouplist[gr].sumstat[st].cat<<"   group="<<gr<<"   samp = "<<this->grouplist[gr].sumstat[st].samp <<"   samp1 = "<<this->grouplist[gr].sumstat[st].samp1 <<"\n";fflush(stdin);}
		else
		{cout <<" calcul de la stat "<<st<<"   cat="<<this->grouplist[gr].sumstat[st].cat<<"   group="<<gr<<"   samp = "<<this->grouplist[gr].sumstat[st].samp <<"   samp1 = "<<this->grouplist[gr].sumstat[st].samp1 <<"   samp2 = "<<this->grouplist[gr].sumstat[st].samp2<<"\n";fflush(stdin);}
		*/
		int categ;
		categ = this->grouplist[gr].sumstat[st].cat;
		//cout<<"categ = "<<categ<<"\n";
		auto&& onestat = stats[categ];
		if (checkStatType(this->grouplist[gr].type,onestat.t)) this->grouplist[gr].sumstat[st].val = onestat.st(*this,gr, st);


//		StatsnpC& stsnp = this->grouplist[gr].sumstatsnp[numsnp];
//		if (stsnp.sw != this->grouplist[gr].nloc) { cout << "categ " << categ << " numsnp " << numsnp << " stsnp.sw " << stsnp.sw  << " nloc : " << this->grouplist[gr].nloc << endl; }
	
		//cout << "      stat["<<st<<"]="<<this->grouplist[gr].sumstat[st].val<<"\n";fflush(stdin);
	}
	if (this->grouplist[gr].type == 1) liberednavar(gr);
	else liberefreq(gr);
	if (not this->afsdone.empty()) {
		for (int sa = 0; sa < this->nsample; sa++) {
			this->afsdone[sa].clear();
			for (int loc = 0; loc < this->grouplist[gr].nloc; loc++) if (not this->t_afs[sa][loc].empty()) this->t_afs[sa][loc].clear();
			this->t_afs[sa].clear();
			this->n_afs[sa].clear();
		}
		this->afsdone.clear();
		this->t_afs.clear();
		this->n_afs.clear();
	}
}
