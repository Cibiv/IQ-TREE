//
//  terraceanalysis.cpp
//  iqtree
//
//  Created by Olga on 10.09.20.
//

#include "terraceanalysis.hpp"
#include "tree/terracetree.hpp"

void runterraceanalysis(Params p){
    
    if(p.user_file && p.pr_ab_matrix){
        MTree *tree_rep = new MTree(p.user_file,p.is_rooted);
        
        
        IntVector pr_ab_matrix = read_pr_ab_matrix(p.pr_ab_matrix);
        TerraceTree *tree = new TerraceTree(pr_ab_matrix);
    }
    
};

presence_absence_matrix read_pr_ab_matrix(const char *file){
    
};
