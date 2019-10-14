/*
 * supermatrix.h
 *
 *  Created on: Sept 1, 2018
 *      Author: Lukasz
 */
#pragma once

#include <string>
#include <set>
#include <bitset>
#include <map>
#include <vector>
#include "node.h"

// get rid of these fixed taxa and loci number
const int M_TAXA = 1320;
const int M_LOCI = 1200;


/*
 *  Presence/abscence matrix per taxon/loci
 */


class Supermatrix {
public:
    Supermatrix();

    void readFile(std::string filename);    // the format is one line per loci with names of taxa "with" available sequences

    std::map<int, std::string> name;
    std::map<std::string, int> index;

    std::map<int, std::bitset<M_TAXA>*> taxon_by_locus;
    std::map<int, std::bitset<M_LOCI>*> locus_by_taxon;

    std::set<std::string> getPartitionTaxa(int i);
    std::string print();

    int total_taxa;
    int total_loci;

    ~Supermatrix();
};
