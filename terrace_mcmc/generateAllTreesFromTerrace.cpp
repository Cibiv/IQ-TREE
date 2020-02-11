//
//  generateAllTreesFromTerrace.cpp
//  iqtree
//
//  Created by Olga on 10.02.20.
//

#include "generateAllTreesFromTerrace.hpp"
#include "terracecheck.hpp"

void generateTerraceTrees(Terrace_IQ* terrace, const char *userOutFile){
    
    vector<string> terraceTrees;
    
    int maxIterations = 1000; // Since the algorithm is exponential, it is neccessary to set a hard limit on the number of iterations to perform and output whatever trees are available from a terrace from these number of iterations. You have to check in simulations what's the appropriate threshold
    
    int count_iter = 0;
    
    MTree tree;
    vector<MTree*> induced_trees;
    
    // identify an induced tree to start from. With the largest number of taxa?
    int i = 0;
    tree.copyTree(terrace->induced_part_trees[i]);
    // here terrace->aln should be substituted by a sub-aln for taxa present in tree
    getALLInducedPartitionTrees(&tree,terrace->aln,induced_trees);
    
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
