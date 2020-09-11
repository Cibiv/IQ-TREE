//
//  presenceabsencematrix.hpp
//  iqtree
//
//  Created by Olga on 10.09.20.
//

#ifndef presenceabsencematrix_hpp
#define presenceabsencematrix_hpp

#include <stdio.h>
#include "utils/tools.h"
#include "tree/mtree.h"

class PresenceAbsenceMatrix {
    
public:
    
    vector<IntVector> pr_ab_matrix;
    vector<string> taxa_names;
    
    int taxa_num;
    int part_num;

    /*
     *  Reading a presence-absence matrix for supermatrix
     */
    void read_pr_ab_matrix(const char *infile);
    void read_pr_ab_matrix(istream &in);
    
    void print_pr_ab_matrix();
};

vector<IntVector> getSubMatrix(vector<IntVector> pr_ab_complete, vector<string> taxa_names, MTree* tree);

#endif /* presenceabsencematrix_hpp */
