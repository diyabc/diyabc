/*
 * randforest.h
 *
 *  Created on: 28 février 2014
 *      Author: jmcornuet
 */

#pragma once

#include <string>
#include <vector>

#include "randomgenerator.hpp"

struct VMC {
    int ind;
    double x;
};

struct VMD {
    std::string name;
    double x;
};

class NodeRC {
public:
    int pere, filsG, filsD, nvar, nsets, npassages, nsetG, nsetD, model, imax;
    double cutval, disval, disval2, delta, modmoy;
    bool terminal, gauche;
    vector<int> indvar;
    vector<int> numset;
    vector<int> numsetG;
    vector<int> numsetD;
    int regle3(vector<VMC> vm, vector<int> b, MwcGen& mw);
    double calGinimin(vector<VMC>& vm, double& cutvalc, vector<int> modfreq) const;
    double calGini(vector<VMC>& vm, double cutval) const;
    double calvarmin(vector<VMC>& vm, double& cutvalc) const;
    double calvarmoy(vector<VMC>& vm, double cutval) const;
    int getdisval(MwcGen& mw);
    double getdisval2(MwcGen& mw);

    ~NodeRC() {
        if (not indvar.empty()) indvar.clear();
        if (not numset.empty()) numset.clear();
        if (not numsetG.empty()) numsetG.clear();
        if (not numsetD.empty()) numsetD.clear();
    }

    NodeRC& operator=(NodeRC const& source);
};

class TreeC {
public:
    int nnodes, nsets, nvar;
    MwcGen mw;
    bool fin;

    vector<int> numset;
    vector<int> score;
    vector<int> index;
    vector<NodeRC> node;
    vector<bool> varused;
    vector<bool> sim_participe;

    void initree();
    int infermodel(const vector<double>& stat);
    void deletree();
    void buildtree1(int seed, int i, int rep);
    void buildtree2(int seed, int i, int rep);
    void estim();
    double inferobs(vector<double>& stat);

    ~TreeC() {
        if (not numset.empty()) numset.clear();
        //if (not indsel.empty()) indsel.clear();
        if (not score.empty()) score.clear();
        if (not index.empty()) index.clear();
        if (not node.empty()) node.clear();
        if (not varused.empty()) varused.clear();
        if (not sim_participe.empty()) sim_participe.clear();
    }

    TreeC& operator=(TreeC const& source);
};

class RFC {
public:
    int ntrees, ntot, nsets, nstat, nmodel, nvar, nsel, nstatclass, nbootsamp;

    vector<int> model;
    vector<vector<double>> vote;
    vector<double> varimp;
    vector<vector<double>> stat;
    vector<vector<double>> importance;
    vector<TreeC> tree;
    vector<double> statobs;
    vector<int> bootsamp;
    vector<string> statname;
    vector<bool> bienestime;

    void readstat(bool LD);

    ~RFC() {
        if (not model.empty()) model.clear();
        if (not vote.empty()) vote.clear();
        if (not varimp.empty()) varimp.clear();
        if (not stat.empty()) stat.clear();
        if (not tree.empty()) tree.clear();
        if (not statobs.empty()) statobs.clear();
        if (not bootsamp.empty()) bootsamp.clear();
        if (not statname.empty()) statname.clear();
    }

    RFC& operator=(RFC const& source);
};

void dorandfor(std::string opt, int seed);
void var_importance3(int rep);

