//
//  generateAllTreesFromTerrace.cpp
//  iqtree
//
//  Created by Olga on 10.02.20.
//

#include "generateAllTreesFromTerrace.hpp"
#include "terrace_iq.hpp"
#include "terracecheck.hpp"

void generateTerraceTrees(Terrace_IQ* terrace, const char *userOutFile){
    
    vector<string> terraceTrees;
    
    int maxIterations = 1000; // Since the algorithm is exponential, it is neccessary to set a hard limit on the number of iterations to perform and output whatever trees are available from a terrace from these number of iterations. You have to check in simulations what's the appropriate threshold
    
    int count_iter = 0;
    
    // identify an induced tree to start from. With the largest number of taxa?
    int i = 0;
    MTree tree;
    tree.copyTree(terrace->induced_part_trees[i]);
    
    vector<IntVector> pr_ab_init;
    vector<string> taxa_names_init;
    
    //YOU STOPPED HERE!!!!!!!!!!!!!!!
    
    // getting pr_ab_matrix and taxa_names for the initial tree
    
    Terrace_IQ terrace_init(pr_ab_init,taxa_names_init,&tree);
    
    // Perform first branch and taxon mappings
    
    // Identify taxa to be inserted and order them by some smart way
    
    // Start generation of trees by inserting the first missing taxon
    
    
    // Print out all generated trees into a file
    ofstream out_file;
    string file_name = userOutFile;
    file_name += ".terrace_trees";
    out_file.open(file_name);
    
    for(int i = 0; i < terraceTrees.size(); i++){
        out_file<<terraceTrees[i]<<endl;
    }
    out_file.close();
}
