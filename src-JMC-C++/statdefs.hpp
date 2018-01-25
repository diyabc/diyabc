#pragma once

#include <string>
#include <functional>
#include <vector>
#include <map>

#include "particuleC.hpp"

using namespace std;

typedef function<double(ParticleC&,int,int)> calstat;

struct stat {
    string name;
    int npop;
    calstat st;
};

vector<stat> microsat_stats = {
    stat{ "PID" , 1, calstat(&ParticleC::cal_pid1p) }, // 0
    stat{ "NAL" , 2, calstat(&ParticleC::cal_nal1p) }, // 1 
    stat{ "HET" , 2, calstat(&ParticleC::cal_het1p) }, // 2
    stat{ "VAR" , 2, calstat(&ParticleC::cal_var1p) }, // 3
    stat{ "MGW" , 2, calstat(&ParticleC::cal_mgw1p) }, // 4
    stat{ "N2P" , 2, calstat(&ParticleC::cal_nal2p) }, // 5
    stat{ "H2P" , 2, calstat(&ParticleC::cal_het2p) }, // 6
    stat{ "V2P" , 2, calstat(&ParticleC::cal_var2p) }, // 7
    stat{ "FST" , 2, calstat(&ParticleC::cal_Fst2p) }, // 8
    stat{ "LIK" , 2, calstat(&ParticleC::cal_lik2p) }, // 9
    stat{ "DAS" , 2, calstat(&ParticleC::cal_das2p) }, // 10
    stat{ "DM2" , 2, calstat(&ParticleC::cal_dmu2p) }, // 11
    stat{ "AML" , 3, calstat(&ParticleC::cal_Aml3p) }  // 12
};

vector<stat> dna_stats {
    stat{ "NHA" , 1, calstat(&ParticleC::cal_nha1p) }, // -1
    stat{ "NSS" , 1, calstat(&ParticleC::cal_nss1p) }, // -2
    stat{ "MPD" , 1, calstat(&ParticleC::cal_mpd1p) }, // -3
    stat{ "VPD" , 1, calstat(&ParticleC::cal_vpd1p) }, // -4
    stat{ "DTA" , 1, calstat(&ParticleC::cal_dta1p) }, // -5
    stat{ "PSS" , 1, calstat(&ParticleC::cal_pss1p) }, // -6
    stat{ "MNS" , 1, calstat(&ParticleC::cal_mns1p) }, // -7
    stat{ "VNS" , 1, calstat(&ParticleC::cal_vns1p) }, // -8
    stat{ "NH2" , 1, calstat(&ParticleC::cal_nha2p) }, // -9
    stat{ "NS2" , 2, calstat(&ParticleC::cal_nss2p) }, // -10
    stat{ "MP2" , 2, calstat(&ParticleC::cal_mpw2p) }, // -11
    stat{ "MPB" , 2, calstat(&ParticleC::cal_mpb2p) }, // -12
    stat{ "HST" , 2, calstat(&ParticleC::cal_fst2p) }, // -13
    stat{ "SML" , 3, calstat(&ParticleC::cal_aml3p) }  // -14
};

typedef function<void(ParticleC&,int,int)> wrapstat;

map<string,calstat> accnames =
    {
        { "P0", calstat(&ParticleC::cal_p0L) },
        { "M1", calstat(&ParticleC::cal_moyL0) },
        { "V1", calstat(&ParticleC::cal_varL0) },
        { "MO", calstat(&ParticleC::cal_moyL) }
    };

struct snpstat {
    string name;
    int npop;
    wrapstat st;
    vector<string> accs;
};

vector<snpstat> snp_stats {
    snpstat{ "Q", 1, // Q1 : 21, 22
                wrapstat(&ParticleC::cal_snhet), 
                { "V1", "MO" }
            }, 
    snpstat{ "H", 1, // HET: 23, 24, 25, 26
                wrapstat(&ParticleC::cal_snhet), 
                { "P0", "M1", "V1", "MO" }
            }, 
    snpstat{ "F", 2, // FST (biaisée) : 27, 28, 29, 30
                wrapstat(&ParticleC::cal_snfst), 
                { "P0", "M1", "V1", "MO" }
            }, 
    snpstat{ "N", 2, // NEI : 31, 32, 33, 34
                wrapstat(&ParticleC::cal_snnei), 
                { "P0", "M1", "V1", "MO" }
            }, 
    snpstat{ "L", 2, // Q2 : 35, 36, 37, 38
                wrapstat(&ParticleC::cal_snq2), 
                { "P0", "M1", "V1", "MO" }
            }, 
    snpstat{ "A", 3, // AML : 39, 40, 41, 42
                wrapstat(&ParticleC::cal_snaml), 
                { "P0", "M1", "V1", "MO" }
            }, 
    snpstat{ "R", 3, // F3 : 43, 44, 45, 46
                wrapstat(&ParticleC::cal_snf3r), 
                { "P0", "M1", "V1", "MO" }
            }, 
    snpstat{ "Z", 4, // F4 : 47, 48, 49, 50
                wrapstat(&ParticleC::cal_snf4r), 
                { "P0", "M1", "V1", "MO" }
            },
    snpstat{ "snFST2", 2, // Fst bi pop}
                wrapstat(&Particle::cal_snfst2),
                { "P0", "M1", "V1", "MO" }
            },
    snpstat{ "snFST3", 3, // Fst tri pop}
                wrapstat(&Particle::cal_snfst3),
                { "P0", "M1", "V1", "MO" }
            },
    snpstat{ "snFST2", 4, // Fst quadri pop}
                wrapstat(&Particle::cal_snfst4),
                { "P0", "M1", "V1", "MO" }
            },
    snpstat{ "snFSTGlob", 0, // Fst Globale}
                wrapstat(&Particle::cal_snfstglob),
                { "P0", "M1", "V1", "MO" }
            }

};