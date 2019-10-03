/*
 * utils.cpp
 *
 *  Created on: Sept 1, 2018
 *      Author: Lukasz
 */

#include "utils.h"

using namespace std;

std::vector<std::string> Utils::getSample(string coverage_file, string tree_file, int size, int burnin)
{
    Parser* p = new Parser();
    MCNode* root = p->from_file(tree_file);    
    delete p;

    Supermatrix* sm = new Supermatrix();
    sm->readFile(coverage_file);

    std::vector<std::string> sample = getSample(sm, root, size, burnin);

    return sample;
}

std::vector<std::string> Utils::getSample(Supermatrix* sm, string newick, int size, int burnin)
{
    Parser* p = new Parser();
    MCNode* root = p->from_newick(newick);
    delete p;

    return getSample(sm, root, size, burnin);    
}

std::vector<std::string> Utils::getSample(Supermatrix* sm,  MCNode* root, int size, int burnin)
{
    SRW* srw = new MHRW();
    
    srw->sm = sm;
    srw->sampling_freq = 1;
    srw->sample_size = 1;

    srw->burnin = burnin;

    srw->command = "print_distinct";

    vector<string> sample;

    for (int i = 0; i < size; i++) {    
        srw->state = root->clone();
        srw->state->resetNNICache();        

        string tree = srw->getSample()->at(0);
        sample.push_back(tree);
    }

    delete srw;
    delete sm;
    delete root;

    return sample;
}
