#include <vector>
#include <iostream>

#include "abcranger.hpp"


int main() {
    vector<float> preds(100);

    for(auto &p : preds) {
        p = abcranger("headerRF.txt","reftableRF.bin","statobsRF.txt",500,0,0);
    }
    cout << endl;
    for(auto p: preds) cout << p << " ";
    cout << endl;
}
