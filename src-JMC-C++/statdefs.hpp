#pragma once

#include <string>
#include <functional>
#include <vector>
#include <optional>

#include "particuleC.hpp"

using namespace std;
using namespace std::placeholders;

template <typename T>
std::vector<T> flatten(const std::vector<std::vector<T>>& v) {
    std::size_t total_size = 0;
    for (const auto& sub : v)
        total_size += sub.size(); // I wish there was a transform_accumulate
    std::vector<T> result;
    result.reserve(total_size);
    for (const auto& sub : v)
        result.insert(result.end(), sub.begin(), sub.end());
    return result;
}

vector<vector<int>> getCombinations(int n, int r);
vector<vector<int>> getArrangementsHalfSortedByPairs(int n, int r);
vector<vector<int>> getArrangements(int n, int r);

enum sortArr {
    COMB, // only combinations
    HALF, // half sorted by pairs
    ALL // every arrangements
};

enum statType {
    MICROSAT,
    DNA,
    SNP
};

bool checkStatType(int grT, statType t);

typedef function<vector<vector<int>>(int,int)> combAlgFun;
typedef function<double(ParticleC&,int,int)> calstatn;
typedef function<void(ParticleC&,int,int)> wrapstatn;


struct snpstatn {
    string name;
    int npop;
    sortArr comb; // Get half-ordered arrangements (true) or just combinations (false)
    wrapstatn st;
    vector<string> accs;

    friend bool operator==(const snpstatn& l, const snpstatn& r)
    {
        return l.name == r.name;
    }
};


struct statn {
    string name;
    int npop;
    sortArr comb; // Get half-ordered arrangements (false) or just combinations (true)
    statType t;
    calstatn st;
    optional<reference_wrapper<snpstatn>> snps;

    friend bool operator==(const statn& l, const statn& r)
    {
        return l.name == r.name;
    }
};

vector<statn> mapStat(snpstatn& st);

