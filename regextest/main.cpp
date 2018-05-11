#include <iostream>
#include <regex>
#include <string>

using namespace std;

string filestr = R"#(
foo
NS1 N UN[1000,30000,15000,7500]
NS2 N UN[1000,30000,15000,7500]
NS3 N LU[10,1000,15000,7500]
NS4 N UN[1000,30000,15000,7500]
NS5 N UN[1000,30000,0.0,0.0]
t1 T UN[40,45,1500,750]
BD1 T UN[0,5,0.0,0.0]
NF1 N LU[2,1000,0.0,0.0]
tm1 T UN[0,3000,0.0,0.0]
t3 T UN[85,320,1500,750]
tm3 T UN[0,3000,0.0,0.0]
t4 T UN[68,73,1500,750]
BD4 T UN[0,5,0.0,0.0]
NF4 N LU[2,1000,0.0,0.0]
tm4 T UN[0,3000,0.0,0.0]
t5 T UN[60,65,0.0,0.0]
BD5 T UN[0,5,0.0,0.0]
NF5 N LU[2,1000,0.0,0.0]
tm5 T UN[0,3000,0.0,0.0]
texp T UN[1000,30000,0.0,0.0]
Nanc N LU[1000,30000,0.0,0.0]
ra4 A UN[0.1,0.9,0.0,0.0]
)#";

int main()
{

    string reparamlistrestr = R"#(foo\n(?:.*\n){22})#";
    const regex reparamlist(reparamlistrestr);
    sregex_token_iterator endregexp;

    // sregex_token_iterator itparamlist(begin(filestr),end(filestr),reparamlist,{1});
    // auto i = 0;
    // for(auto it = itparamlist; it != endregexp; it++) {
    //     cout << i++ << " : " << (*it);
    // }
    
    smatch base_match;
    if (regex_search(filestr,base_match,reparamlist,regex_constants::match_not_null)) {
        auto i = 0;
        for(auto it= next(base_match.begin()); it != base_match.end(); it++) {
            cout << i++ << " : " << (*it) << endl;
        }
        cout << endl;
        // cout << "result : " << base_match[1] << endl;
    } else { 
        cout << "not matched" << endl;
    }
    cout.flush();
    return 0;
}
