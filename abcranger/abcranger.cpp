#include <iostream>
#include <fstream>
#include <iterator>
#include <regex>
#include <algorithm>
#include <iomanip>
#include <map>
#include <unordered_set>

#include "DataFloat.h"

using namespace std;

template <typename InputIt>
void print(InputIt it, InputIt end_it)
{
    while (it != end_it) {
        const string scenn(*it++);
        if (it == end_it) { break; }
        const string scendesc(*it++);
        cout << left << setw(28) << scenn
             << " : " << scendesc << endl;
    }
}
int N = 0;

uint32_t ntohl(uint32_t netlong)
{
    union {
        uint16_t num;
        uint8_t bytes[2];
    } static endian_test = { .bytes = { 0x01, 0x00 }};

    if (endian_test.num == 0x0001) {
        netlong = (netlong << 24) | ((netlong & 0xFF00ul) << 8) |
            ((netlong & 0xFF0000ul) >> 8) | (netlong >> 24);
    }

    return netlong;
}

int main() {
    ifstream headerStream("headerRF.txt",ios::in);
    headerStream >> noskipws;
    const std::string hS(istream_iterator<char>{headerStream}, {});
    
    const regex scen_re("\\bscenario\\s+(\\d+)\\s+.*\\n((?:(?!(?:scenario|\\n)).*\\n)+)");
    sregex_token_iterator itscen(begin(hS), end(hS), scen_re, {1,2});
    //    print(itscen, {});

    auto nscenh = distance(itscen, {})/2;
    vector<string> scendesc(nscenh);
    sregex_token_iterator endregexp;
    
    auto it = itscen;
    while (it != endregexp) {
        const string num(*it++);
        const string desc(*it++);
        scendesc[stoi(num)-1] = desc;
    }
    //cout << scendesc[0];
    const regex nparamtot("historical parameters priors \\((\\d+)\\D");
    smatch base_match;
    regex_search(begin(hS), end(hS),base_match,nparamtot);
    auto nparamtoth = stoi(base_match[1]);
//    cout << nparamtoth << endl;
    const regex reparamlist("\\bhistorical parameters priors.*\\n((?:\\w+\\W[^\\n]*\\n){" + to_string(nparamtoth) + "})");
    regex_search(begin(hS), end(hS),base_match,reparamlist);
    //cout << base_match[1] << endl;
    const string paramlistmatch = base_match[1];
    const regex reparam("(\\w+)\\W+\\w\\W+\\w\\w\\[([^,\\]]+),([^,\\]]+)[,\\]][^\\n]*\\n");
    sregex_token_iterator reparamit(begin(paramlistmatch),end(paramlistmatch), reparam, {1,2,3});
    it = reparamit;
    int reali = 1;
    map<string,double> paramdesc;
    while (it != endregexp) {
        const string param = *it++;
        const double mini = stod(*it++);
        const double maxi = stod(*it++);
        if (maxi != 0.0) {
            if ((maxi-mini)/maxi > 0.000001) {
                paramdesc[param] = reali;
                reali++;
            }
        }
//        cout << "blah : " << (*it++) << " blih : " << (*it++) << " blouh: " << (*it++) << endl;
    }
//    cout << paramdesc["BD4"] << endl;
    int realparamtot = reali - 1;
    vector<vector<int>> parambyscenh(nscenh);
    const regex splitre("\\W");
    for(int i = 0; i < nscenh; i++) {
        sregex_token_iterator it(begin(scendesc[i]),end(scendesc[i]),splitre,-1);
        for(; it != endregexp; ++it) {
            const string term = *it;
            if (paramdesc.count(term) > 0) {
                const int nterm = paramdesc[term];
                if (find(parambyscenh[i].begin(),parambyscenh[i].end(),nterm)  == parambyscenh[i].end())
                    parambyscenh[i].push_back(nterm);
            } 
        }
    }
    // for (auto scen : parambyscenh) {
    //     for(auto p : scen) cout << p << " ";
    //     cout << endl;       
    // }
    const regex restatsname("\\n\\s*\\nscenario\\s+");
    regex_search(begin(hS), end(hS),base_match,restatsname);
    const string allstatsname = base_match.suffix();
    const regex splitre2("\\s+");
    vector<string> allcolspre;
    for(sregex_token_iterator it(allstatsname.begin(),allstatsname.end(),splitre2,-1); it != endregexp; it++)
        allcolspre.push_back(*it);

    ifstream reftableStream("reftableRF.bin",ios::in|ios::binary);
    int realnrec;
    reftableStream.read(reinterpret_cast<char *>(&realnrec),sizeof(realnrec));
    int nrec = N > 0 ? min(realnrec,N) : realnrec;
    int nscen;
    reftableStream.read(reinterpret_cast<char *>(&nscen),sizeof(nscen));    
    vector<int> nrecscen(nscen);
    for(auto& r : nrecscen) reftableStream.read(reinterpret_cast<char *>(&r),sizeof(r));
    vector<int> nparam(nscen);
    for(auto& r : nparam) reftableStream.read(reinterpret_cast<char *>(&r),sizeof(r));
    int nstat;
    reftableStream.read(reinterpret_cast<char *>(&nstat),sizeof(nstat));    
    vector<string> paramsname { allcolspre.begin(), allcolspre.begin() + (allcolspre.size() - nstat) };
    vector<string> statsname  { allcolspre.begin()+ (allcolspre.size() - nstat),  allcolspre.end() };
    int nmutparams = paramsname.size() - realparamtot;
    vector<float> params(nrec * paramsname.size(),NAN);
    // for(auto r : statsname) cout << r << endl;
    vector<float> stats(nrec * nstat,NAN);
    vector<int> scenarios(nrec);
    for(auto i = 0; i < nrec; i++) {
        int scen;
        reftableStream.read(reinterpret_cast<char *>(&scen),4);
        scenarios[i] = scen;
        scen--;
        vector<float> lparam(nparam[scen]);
        for(auto& r: lparam) {
            reftableStream.read(reinterpret_cast<char *>(&r),sizeof(r));
        }
        for(auto j = 0; j < parambyscenh[scen].size(); j++)
            params[i * paramsname.size() + parambyscenh[scen][j] - 1] = lparam[j];
        for(auto j = 0; j < nmutparams; j++)
            params[i * paramsname.size() + realparamtot + j - 1] = lparam[nparam[scen] - nmutparams + j - 1];
        for_each(stats.begin() + (i * nstat), stats.begin() + (i + 1) * nstat, 
                 [&](float& r) { reftableStream.read(reinterpret_cast<char *>(&r),4);});        
    }
    cout << endl;
    cout.flush();
}