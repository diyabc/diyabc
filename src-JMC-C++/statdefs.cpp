#include <map>
#include <numeric>
#include <algorithm>

#include "statdefs.hpp"

/////////////

bool checkStatType(int grT, statType t) {
	return ((grT == 0) && (t == statType::MICROSAT)) ||
		((grT == 1) && (t == statType::DNA)) ||
		((grT == 2 || grT == 3) && (t == statType::SNP));
}


vector<vector<int>> getCombinations(int n, int r) {
    vector<bool> selector(n);
    fill(selector.begin(), selector.begin() + r, true);
    vector<int> v(n);
    iota(v.begin(),v.end(),0);
    vector<int> s(r);
    vector<vector<int>> res;
    do
    {
        copy_if(v.begin(), v.end(), s.begin(), [&](const int &i) { return selector[i]; });
        res.push_back(s);
    } while (prev_permutation(selector.begin(), selector.end()));
    return res;
} 

vector<vector<int>> subGetArrangements(vector<int>& s) {
    vector<vector<int>> res;
    int r = s.size();
    do {
        res.push_back(s);
    } while (next_permutation(s.begin(), s.end()));
    return res;
}

bool halfsortedbypairs(const vector<int>& v) {
    bool res = true;
    int n = v.size();
    for(int i = n-1; i > 0; i -= 2) {
        res = res && (v[i-1] <= v[i]);
        if ((i-2) > 0) res = res && (v[i-3] <= v[i-1]);
    }
    return res;
}

vector<vector<int>> subGetArrangementsHalSortedByPairs(vector<int>& s) {
    vector<vector<int>> res;
    int r = s.size();
    do {
        if (halfsortedbypairs(s)) res.push_back(s);
    } while (next_permutation(s.begin(), s.end()));
    return res;
}

vector<vector<int>> getArrangementsHalfSortedByPairs(int n, int r) {
    auto combis = getCombinations(n, r);
    vector<vector<vector<int>>> arrs;
    transform(combis.begin(),combis.end(),back_inserter(arrs),
            subGetArrangementsHalSortedByPairs);
    auto arrsF = flatten(arrs);
    return arrsF;
}

vector<vector<int>> getArrangements(int n, int r) {
    auto combis = getCombinations(n, r);
    vector<vector<vector<int>>> arrs;
    transform(combis.begin(),combis.end(),back_inserter(arrs),
            subGetArrangements);
    auto arrsF = flatten(arrs);
    return arrsF;
}

// Types of sort for getting Pops combinations

map<sortArr, combAlgFun> sortAlgos 
{
    { COMB, combAlgFun(getCombinations) },
    { HALF, combAlgFun(getArrangementsHalfSortedByPairs) },
    { ALL,  combAlgFun(getArrangementsHalfSortedByPairs )}
};


map<string, calstatn> accnames =
    {
        { "P0", calstatn(&ParticleC::cal_p0L) },
        { "M1", calstatn(&ParticleC::cal_moyL0) },
        { "V1", calstatn(&ParticleC::cal_varL0) },
        { "MO", calstatn(&ParticleC::cal_moyL) },
        { "MF", calstatn(&ParticleC::cal_fstmoyL)},
        { "MQ", calstatn(&ParticleC::cal_fsti)}
    };


vector<statn> microsat_statns = {
    statn{ "PID" , 1, sortArr::HALF, statType::MICROSAT, calstatn(&ParticleC::cal_pid1p) }, // 0
    statn{ "NAL" , 2, sortArr::HALF, statType::MICROSAT, calstatn(&ParticleC::cal_nal1p) }, // 1 
    statn{ "HET" , 2, sortArr::HALF, statType::MICROSAT, calstatn(&ParticleC::cal_het1p) }, // 2
    statn{ "VAR" , 2, sortArr::HALF, statType::MICROSAT, calstatn(&ParticleC::cal_var1p) }, // 3
    statn{ "MGW" , 2, sortArr::HALF, statType::MICROSAT, calstatn(&ParticleC::cal_mgw1p) }, // 4
    statn{ "N2P" , 2, sortArr::HALF, statType::MICROSAT, calstatn(&ParticleC::cal_nal2p) }, // 5
    statn{ "H2P" , 2, sortArr::HALF, statType::MICROSAT, calstatn(&ParticleC::cal_het2p) }, // 6
    statn{ "V2P" , 2, sortArr::HALF, statType::MICROSAT, calstatn(&ParticleC::cal_var2p) }, // 7
    statn{ "FST" , 2, sortArr::HALF, statType::MICROSAT, calstatn(&ParticleC::cal_Fst2p) }, // 8
    statn{ "LIK" , 2, sortArr::ALL, statType::MICROSAT, calstatn(&ParticleC::cal_lik2p) }, // 9
    statn{ "DAS" , 2, sortArr::HALF, statType::MICROSAT, calstatn(&ParticleC::cal_das2p) }, // 10
    statn{ "DM2" , 2, sortArr::HALF, statType::MICROSAT, calstatn(&ParticleC::cal_dmu2p) }, // 11
    statn{ "AML" , 3, sortArr::HALF, statType::MICROSAT, calstatn(&ParticleC::cal_Aml3p) }  // 12
};

vector<statn> dna_statns {
    statn{ "NHA" , 1, sortArr::HALF, statType::DNA, calstatn(&ParticleC::cal_nha1p) }, // -1
    statn{ "NSS" , 1, sortArr::HALF, statType::DNA, calstatn(&ParticleC::cal_nss1p) }, // -2
    statn{ "MPD" , 1, sortArr::HALF, statType::DNA, calstatn(&ParticleC::cal_mpd1p) }, // -3
    statn{ "VPD" , 1, sortArr::HALF, statType::DNA, calstatn(&ParticleC::cal_vpd1p) }, // -4
    statn{ "DTA" , 1, sortArr::HALF, statType::DNA, calstatn(&ParticleC::cal_dta1p) }, // -5
    statn{ "PSS" , 1, sortArr::HALF, statType::DNA, calstatn(&ParticleC::cal_pss1p) }, // -6
    statn{ "MNS" , 1, sortArr::HALF, statType::DNA, calstatn(&ParticleC::cal_mns1p) }, // -7
    statn{ "VNS" , 1, sortArr::HALF, statType::DNA, calstatn(&ParticleC::cal_vns1p) }, // -8
    statn{ "NH2" , 1, sortArr::HALF, statType::DNA, calstatn(&ParticleC::cal_nha2p) }, // -9
    statn{ "NS2" , 2, sortArr::HALF, statType::DNA, calstatn(&ParticleC::cal_nss2p) }, // -10
    statn{ "MP2" , 2, sortArr::HALF, statType::DNA, calstatn(&ParticleC::cal_mpw2p) }, // -11
    statn{ "MPB" , 2, sortArr::HALF, statType::DNA, calstatn(&ParticleC::cal_mpb2p) }, // -12
    statn{ "HST" , 2, sortArr::HALF, statType::DNA, calstatn(&ParticleC::cal_fst2p) }, // -13
    statn{ "SML" , 3, sortArr::HALF, statType::DNA, calstatn(&ParticleC::cal_aml3p) }  // -14
};


vector<snpstatn> snp_statns {
    snpstatn{ "QW", 1, sortArr::HALF, // Q1 : 21, 22
                wrapstatn(&ParticleC::cal_snq1), 
                { "V1", "MO" }
            }, 
    // snpstatn{ "H", 1, sortArr::HALF, // HET: 23, 24, 25, 26
    //             wrapstatn(&ParticleC::cal_snhet), 
    //             { "P0", "M1", "V1", "MO" }
    //         }, 
    snpstatn{ "H", 1, sortArr::HALF, // HET: 23, 24, 25, 26
                wrapstatn(&ParticleC::cal_snhet), 
                { "MO" }
            }, 
    // snpstatn{ "F", 2, sortArr::HALF, // FST (biaisée) : 27, 28, 29, 30
    //             wrapstatn(&ParticleC::cal_snfst), 
    //             { "P0", "M1", "V1", "MO" }
    //         }, 
    // snpstatn{ "NEI", 2, sortArr::HALF, // NEI : 31, 32, 33, 34
    //             wrapstatn(&ParticleC::cal_snnei), 
    //             { "P0", "M1", "V1", "MO" }
    //         }, 
    snpstatn{ "QB", 2, sortArr::HALF, // Q2 : 35, 36, 37, 38
                wrapstatn(&ParticleC::cal_snq2), 
                { "P0", "M1", "V1", "MO" }
            }, 
    // snpstatn{ "A", 3, sortArr::HALF, // AML : 39, 40, 41, 42
    //             wrapstatn(&ParticleC::cal_snaml), 
    //             { "P0", "M1", "V1", "MO" }
    //         }, 
    // snpstatn{ "F3", 3, sortArr::HALF, // F3 : 43, 44, 45, 46
    //             wrapstatn(&ParticleC::cal_snf3r), 
    //             { "P0", "M1", "V1", "MO" }
    //         }, 
    // snpstatn{ "F4", 4, sortArr::HALF, // F4 : 47, 48, 49, 50
    //             wrapstatn(&ParticleC::cal_snf4r), 
    //             { "P0", "M1", "V1", "MO" }
    //         },
    snpstatn{ "FSTI", 1, sortArr::COMB, // Fst mono pop },
                wrapstatn(&ParticleC::ParticleC::cal_snfsti),
                { "MQ" }
            },
    snpstatn{ "FST2", 2, sortArr::COMB, // Fst bi pop}
                wrapstatn(bind(&ParticleC::cal_snfstd,_1,_2,_3,2)),
                { "MF" }
            },
    snpstatn{ "FST3", 3, sortArr::COMB, // Fst tri pop}
                wrapstatn(bind(&ParticleC::cal_snfstd,_1,_2,_3,3)),
                { "MF" }
            },
    snpstatn{ "FST4", 4, sortArr::COMB, // Fst quadri pop}
                wrapstatn(bind(&ParticleC::cal_snfstd,_1,_2,_3,4)),
                { "MF" }
            },
    snpstatn{ "FSTG", 0, sortArr::COMB, // Fst Globale}
                wrapstatn(bind(&ParticleC::cal_snfstd,_1,_2,_3,0)),
                { "MF" }
            }

};

vector<statn> mapStat(snpstatn& st) {
    vector<statn> res;
    for(auto&& s: st.accs) {
        calstatn ms = [&s,&st] (ParticleC& p, int gr, int numst) {
            int numsnp = p.grouplist[gr].sumstat[numst].numsnp; 
            if (not p.grouplist[gr].sumstatsnp[numsnp].defined) {
                p.cal_snpstatRedinit(gr, numsnp);
                st.st(p, gr, numsnp); 
                p.grouplist[gr].sumstatsnp[numsnp].defined = true;
            }
            return accnames[s](p, gr, numsnp); 
        }; 
        res.push_back(statn {st.name + s, st.npop, st.comb, statType::SNP, ms, ref(st) });
    }
    return res;
}
