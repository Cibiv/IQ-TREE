/*
 * srw.h
 *
 *  Created on: Sept 1, 2018
 *      Author: Lukasz
 */
#pragma once

#include <string>
#include <set>
#include <vector>
#include "node.h"
#include "supermatrix.h"

class SRW {
public:    
    std::vector<std::string>* getSample();

    std::vector<std::string> sample;

    std::set<std::string> distinct_nwk;

    //std::vector<int> degree_counts;

    MCNode* state;

    Supermatrix* sm;

    int burnin = 0;

    int sample_size = 1;

    int sampling_freq = 1;

    std::string command;

    bool print_degrees = false;

    bool print_trees = false;

    bool print_unrooted = false;

    ~SRW();

protected:
    void processSampleAddition();

    void virtual transition();
};

