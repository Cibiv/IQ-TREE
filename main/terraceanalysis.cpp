//
//  terraceanalysis.cpp
//  iqtree
//
//  Created by Olga on 10.09.20.
//

#include "terraceanalysis.hpp"
#include "terrace/terracenode.hpp"
#include "terrace/terracetree.hpp"
#include "terrace/terrace.hpp"
#include "terrace/presenceabsencematrix.hpp"

void runterraceanalysis(Params &params){
    
    if(params.user_file && params.pr_ab_matrix){
        
        Terrace *terrace = new Terrace(params.user_file,params.is_rooted,params.pr_ab_matrix);
        
        terrace->linkTrees();
        
        terrace->printMapInfo();
        
        /**
        PresenceAbsenceMatrix matrix;
        matrix.read_pr_ab_matrix(params.pr_ab_matrix);
        matrix.print_pr_ab_matrix();
        
        TerraceTree tree;
        tree.readTree(params.user_file,params.is_rooted);
        cout<<"Terrace representative tree:"<<endl;
        tree.printTree(cout);
        cout<<endl<<endl;
         */
        
        
        
       /**
        Neighbor *nei1 = tree.root->neighbors[0]->node->neighbors[0];
        Neighbor *nei2 = tree.root->neighbors[0]->node->neighbors[1];
        Neighbor *nei3 = tree.root->neighbors[0]->node->neighbors[2];
        
        ((TerraceNeighbor*)nei1)->link_neighbors.push_back(nei2);
        ((TerraceNeighbor*)nei1)->link_neighbors.push_back(nei3);
        
        ((TerraceNeighbor*)nei1)->taxa_to_insert.push_back(nei2->node);
        ((TerraceNeighbor*)nei1)->taxa_to_insert.push_back(nei3->node);
        
        Node *dad = tree.root->neighbors[0]->node;
        ((TerraceNeighbor*)nei1)->printInfo(dad);
        */
        
    } else {
        
        PresenceAbsenceMatrix matrix;
        matrix.read_pr_ab_matrix(params.pr_ab_matrix);
        matrix.print_pr_ab_matrix();
    }
    
};

