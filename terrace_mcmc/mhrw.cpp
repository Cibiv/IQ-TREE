/*
 * mhrw.cpp
 *
 *  Created on: Sept 1, 2018
 *      Author: Lukasz
 */
#include "mhrw.h"
#include "iostream"
#include <random>
#include <map>

using namespace std;

void MHRW::transition() 
{
    static random_device seed;
    static mt19937 engine( seed( ) );

    static uniform_int_distribution<int> upto2(0, 1);
    static uniform_real_distribution<> dis(0.0, 1.0);

    state->getOnTerraceNNIs(sm);

    if (state->num_nnis == 0) {
        return;
    }

    int deg_current = state->num_nnis * 2;
    uniform_int_distribution<int> choose(0, state->num_nnis - 1);

    int ind = choose(engine) * 4;
    MCNode* nni[4];    

    for (int i =0; i<4; i++) {
        nni[i] = state->nni_cache[ind + i];    
    }

    bool nni_type = (upto2(engine) % 2 == 0);
    state->performNNI(nni, nni_type);

    state->getOnTerraceNNIs(sm);
    int deg = state->num_nnis * 2;

    double p = dis(engine);
    double frac = (double)deg_current/deg;    

    if (p <= frac) {
        // all good, transitioning to new state                
    } else {
        // revert
        state->performNNI(nni, nni_type);
    }    
}
