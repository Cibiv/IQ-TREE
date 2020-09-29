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
        
        /*  INFO:
         *  This is an actual terrace:
         *  - main tree - is terrace representative tree
         *  - induced partition trees define the terrace
         *
         *  Create an auxiliary terrace:
         *  - main tree - is the initial tree to be expanded by adding taxa to obtain a tree from the considered terrace (or to get into dead end)
         *  - induced trees - the common subtrees of the main tree and higher level induced partition tree, respectivelly per partition
         *
         *  There will be also a vector of terraces with just one partition. Per partition each terrace is a pair of high and lower induced partition trees
         *  - main tree - high level partition tree
         *  - induced tree - a common subtree between "initial to be expanded main tree" and high level induced tree
         */
        Terrace *terrace = new Terrace(params.user_file,params.is_rooted,params.pr_ab_matrix);
        
        //terrace->linkTrees();
        //terrace->printMapInfo();
        //terrace-> printBackMapInfo();
        
        /*  TODO:
         *  1. get submatrix (for testing purposes, any submatrix of original presence-absence matrix)
         *  2. get initial tree (for testing purposes, any subtree of the representative tree?)
         *  3. get low-level induced partition trees (i.e. create a sub-terrace)
         *  4. get a vector of terraces for high- and low-level induced partition trees
         */
        
        
        
        
        
        
        
        
        
        
        
        // BELOW stuff was only used for testing. I think, you can delete it.
        
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

