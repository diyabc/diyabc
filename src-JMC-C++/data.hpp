/*
 * data.h
 *
 *  Created on: 24 janvier 2015
 *      Author: jmc
 *
*/
#pragma once

#include <string>
#include <vector>
#include <ciso646>

/**
*  Structure LocusC : définition de la structure LocusC
*/
class LocusC {
public:
	std::string name;
	int type; //0 à 14
	int groupe; //numero du groupe auquel appartient le locus
	double coeffcoal; // coefficient pour la coalescence (dépend du type de locus et du sexratio)
	std::vector<std::vector<long double>> freq;
	//Proprietes des locus sequences
	long double pi_A, pi_C, pi_G, pi_T;
	std::vector<long double> mutsit; //array of dnalength elements giving the relative probability of a mutation to a given site of the sequence
	std::vector<int> sitmut; //array of nsitmut dna sites that are changed through a mutation
	std::vector<int> sitmut2; //
	int dnalength, dnavar;
	std::vector<int> tabsit; //array of dnalength elements giving the number of a dna site;
	std::vector<std::vector<std::string>> haplodna; //array[sample][gene copy][nucleotide] tous les nucleotides de chaque individu sont mis à la suite les uns des autres
	std::vector<std::vector<std::string>> haplodnavar; //seulement les sites variab
	//Proprietes des locus microsatellites
	int mini, maxi, kmin, kmax, motif_size, motif_range, nal;
	double mut_rate, Pgeom, sni_rate, mus_rate, k1, k2;
	std::vector<std::vector<int>> haplomic; //array[sample][gene copy]
	//Propriétés des locus SNP
	bool firstime;
	std::vector<std::vector<short int>> haplosnp; //array[sample][gene copy] 0,1,9
	bool mono; //mono=true si un seul allèle dans l'échantillon global
	double weight; //poids du locus=1, sauf quand biais de recrutement
	int nsample, nmisssnp;
	std::vector<int> samplesize;
	std::vector<std::vector<short int>> ploidie;
	//Propriétés des "locus" PoolSeq
	std::vector<int> nreads1;
	std::vector<int> nreads;

	LocusC() {
		name = "";
		type = -1;
		groupe = -1;
		pi_A = pi_C = pi_G = pi_T = -1.0;
		dnalength = dnavar = -1;
		mini = maxi = kmin = kmax = -1;
		motif_range = motif_size = -1;
		mut_rate = Pgeom = sni_rate = mus_rate = k1 = k2 = 0.0;
		firstime = mono = true;
		weight = 1.0;
		nsample = 0;
		nmisssnp = 0;
	};

	~LocusC() {
		if (! mutsit.empty()) mutsit.clear();
		if (! sitmut.empty()) sitmut.clear();
		if (! sitmut2.empty()) sitmut2.clear();
		if (! tabsit.empty()) tabsit.clear();
		if (! freq.empty()) {
			int imax = (int)freq.size();
			for (int i = 0; i < imax; i++) {
				if (! freq[i].empty()) freq[i].clear();
			}
		}
		if (! haplomic.empty()) {
			int imax = (int)haplomic.size();
			for (int i = 0; i < imax; i++) {
				if (! haplomic[i].empty()) haplomic[i].clear();
			}
		}
		if (! haplodna.empty()) {
			int imax = (int)haplodna.size();
			for (int i = 0; i < imax; i++) {
				if (! haplodna[i].empty()) haplodna[i].clear();
			}
		}
		if (! haplodnavar.empty()) {
			int imax = (int)haplodnavar.size();
			for (int i = 0; i < imax; i++) {
				if (! haplodnavar[i].empty()) haplodnavar[i].clear();
			}
		}
		if (! haplosnp.empty()) {
			int imax = (int)haplosnp.size();
			for (int i = 0; i < imax; i++) {
				if (! haplosnp[i].empty()) haplosnp[i].clear();
			}
		}
		if (! ploidie.empty()) {
			int imax = (int)ploidie.size();
			for (int i = 0; i < imax; i++) {
				if (! ploidie[i].empty()) ploidie[i].clear();
			}
		}
		if (! samplesize.empty()) samplesize.clear();
		if (! nreads1.empty()) nreads1.clear();
		if (! nreads.empty()) nreads.clear();

	};

	LocusC& operator=(LocusC const& source);

	void libere(bool obs, int nsample);
};

class MissingHaplo {
public:
	int locus, sample, indiv;
	MissingHaplo& operator=(MissingHaplo const& source);
};

class MissingNuc {
public:
	int locus, sample, indiv, nuc;
	MissingNuc& operator=(MissingNuc const& source);
};


class DataC {
public:
	std::string message, title;
	std::vector<std::vector<std::string>> indivname;
	int nsample, nsample0, nloc, nmisshap, nmissnuc, nmisssnp, filetype, mrc;
	double sexratio, maf, lambdapoolseq;
	std::vector<MissingHaplo> misshap;
	std::vector<MissingNuc> missnuc;
	std::vector<bool> catexist;
	std::vector<std::vector<int>> ssize;//nombre de copies de gènes (manquantes incluses) par [locustype][sample], locustype variant de 0 à 4.
	std::vector<int> nind;
	std::vector<std::vector<int>> indivsexe;
	std::vector<LocusC> locus;

	/* Méthodes */

	DataC& operator=(DataC const& source);

	/**
	 * détermination du type de fichier de donnée
	 * return=-1
	 * return=0 si genepop
	 * return=1 si snp
     * return=2 si PoolSeq
	 */
	int testfile(std::string filename);
	/**
	 * lecture d'un fichier de donnée SNP et stockage des informations dans une structure DataC
	 */
	int readfilesnp(std::string filename);
	/**
	 * lecture d'un fichier de donnée PoolSeq et stockage des informations dans une structure DataC
	 */
	int readfilePoolSeq(std::string filename);
	/**
	 * supprime les locus monomorphes
	 */
	void purgelocmonomorphes();
	/**
	 * supprime les locus en dessous de la MAF ou du MRC
	 */
	void purgelocMRCPOOLSEQ();
	void missingdata();
	/**
	 * traitement des locus snp
	 */
	void do_snp(int loc);


	/**
	 * ecriture en binaire d'un fichier snp
	 */
	void ecribin(std::string filenamebin, std::string filenametxt);
	/**
	 * lecture en binaire d'un fichier snp
	 */
	void libin(std::string filenamebin);
	/**
	 * lecture d'un fichier de donnée GENEPOP et stockage des informations dans une structure DataC
	 */
	int readfile(std::string filename);
	/**
	 * traitement des génotypes microsat au locus loc
	 */
	void do_microsat(int loc);
	/**
	 * traitement des génotypes DNA sequence au locus loc
	 */
	void do_sequence(int loc);
	/**
	 * calcul du coefficient dans la formule de coalescence en fonction du type de locus
	 * 0:autosomal diploide, 1:autosomal haploïde, 2:X-linked, 3:Y-linked, 4:mitochondrial
	 */
	void cal_coeffcoal(int loc);

	void calcule_ss();

	void calploidie();
	/**
	 * chargement des données dans une structure DataC
	 */
	int loadfromfile(std::string filename);
};
