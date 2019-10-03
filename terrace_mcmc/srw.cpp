/*
 * srw.cpp
 *
 *  Created on: Sept 1, 2018
 *      Author: Lukasz
 */

#include <iostream>
#include <string>
#include "srw.h"

using namespace std;

void SRW::processSampleAddition()
{
    static string root_label = state->getFirstLeaf();

    state->getOnTerraceNNIs(sm);

    int nni_degree = state->num_nnis;

    MCNode* clone = state->clone();
    clone->root(clone, root_label);
    string nwk = clone->printSorted(true);//!print_unrooted);
    sample.push_back(nwk);

    if (print_degrees) {
        cout << nwk << " " << nni_degree << endl;
        //cout << nwk << " " << nni_degree << " " << clone->getTerraceId(sm) << endl;
    }

    delete clone;
}

vector<string>* SRW::getSample()
{
    sample.clear();

    for (int i = 0; i < burnin; i++) {
        this->transition();
    }

    for (int i = 0; i < sample_size * sampling_freq; i++) {
        if (i % sampling_freq == 0) {
            processSampleAddition();            
        }

        this->transition();
    }

    return &sample;
}

void SRW::transition()
{
    state->getOnTerraceNNIs(sm);

    if (state->num_nnis == 0) {
        return;
    }

    int ind = (rand() % state->num_nnis);

    MCNode* nni[4];

    for (int i =0; i<4; i++) {
        nni[i] = state->nni_cache[ind*4 + i];    
    }

    bool whichNNI = (rand() % 2 == 0);

    state->performNNI(nni, whichNNI);
}

SRW::~SRW()
{
    if (state != NULL)
        delete this->state;
}
