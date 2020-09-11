//
//  terraceanalysis.cpp
//  iqtree
//
//  Created by Olga on 10.09.20.
//

#include "terraceanalysis.hpp"
#include "terrace/terracetree.hpp"
#include "terrace/presenceabsencematrix.hpp"

void runterraceanalysis(Params &params){
    
    if(params.user_file && params.pr_ab_matrix){
        
        PresenceAbsenceMatrix matrix;
        matrix.read_pr_ab_matrix(params.pr_ab_matrix);
        matrix.print_pr_ab_matrix();
        
        //MTree *tree_rep = new MTree(p.user_file,p.is_rooted);
        
        
//        IntVector pr_ab_matrix = read_pr_ab_matrix(p.pr_ab_matrix);
//        TerraceTree *tree = new TerraceTree(pr_ab_matrix);
    }
    
};

