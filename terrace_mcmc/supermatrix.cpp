/*
 * supermatrix.cpp
 *
 *  Created on: Sept 1, 2018
 *      Author: Lukasz
 */

#include <iostream>
#include <fstream>
#include <sstream>
#include "supermatrix.h"

using namespace std;

Supermatrix::Supermatrix()
{
    
};

std::string Supermatrix::print()
{
    string s = "";

    int counts[MAX_LOCI]; 

    for (int j=0; j<taxon_by_locus.size(); j++) {
        counts[j] = 0;
    }

    for (int i=0; i<locus_by_taxon.size(); i++) { // for all taxa
        for (int j=0; j<taxon_by_locus.size(); j++) { // for all loci
            if (locus_by_taxon[i]->test(j)) { 
                counts[j]++;
                s += "1";
            } else {
                s += "0";
            }

            s += " ";
        }

        s += name[i];
        s += "\n";
    }

    for (int j=0; j<taxon_by_locus.size(); j++) {
        s += to_string(counts[j]);
        s += " ";
    }

    s += "TOTAL";
    s += "\n";

    return s;
};

set<string> Supermatrix::getPartitionTaxa(int k)
{
    set<string> taxa;

    //cout << "Partition " << k << " taxa: " << endl;

    for (int i=0; i<total_taxa; i++) {
        if (taxon_by_locus[k]->test(i)) {            
            taxa.insert(name[i]);            
            //cout << " - : " << name[i] << endl;
        }
    }

    return taxa;
}

void Supermatrix::readFile(string filename)
{
    ifstream infile(filename.c_str());

    string line;

    getline(infile, line);

    istringstream is(line);

    is >> total_taxa;
    is >> total_loci;

    int ind = 0;

    //cout << "TT " << total_taxa << endl;
    //cout << "TL " << total_loci << endl;
    //cout << "FN " << filename << endl;

    while (getline(infile, line))
    {
        stringstream linestream;
        linestream << line;

        vector<string> elements;
        string el;

        while (linestream >> el) 
        {            
            elements.push_back(el);
        }

        string taxon = elements[elements.size()-1];

        index[taxon] = ind;
        name[ind] = taxon;

        if (locus_by_taxon[ind] == NULL) {
            locus_by_taxon[ind] = new bitset<MAX_LOCI>();
        }    

        for (int i=0; i<elements.size()-1; i++) {
            if (elements[i] == "1") {
                if (taxon_by_locus[i] == NULL) {
                    taxon_by_locus[i] = new bitset<MAX_TAXA>();
                }        

                taxon_by_locus[i]->set(ind);  
                locus_by_taxon[ind]->set(i);
            }
        }

        ind++;
    }
};

Supermatrix::~Supermatrix()
{
    for (int i=0; i<taxon_by_locus.size(); i++) {
        delete taxon_by_locus[i];
    }

    for (int i=0; i<locus_by_taxon.size(); i++) {
        delete locus_by_taxon[i];
    }    
}
