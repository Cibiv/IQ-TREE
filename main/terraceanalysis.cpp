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
    
    int i=0, j=0;
    
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
        //terrace->linkTrees(true,true);
        //terrace->printMapInfo();
        //terrace-> printBackMapInfo();
        
        vector<string> taxa_names_sub;
        
        //exit(0);
        
        // TODO: At this point taxa_names_sub is just a test example. Later on get initial tree by a smart choice.
        // chosen for testing, will not produce all trees from a terrace
        /*taxa_names_sub.push_back(terrace->matrix->taxa_names[2]);
        taxa_names_sub.push_back(terrace->matrix->taxa_names[5]);
        taxa_names_sub.push_back(terrace->matrix->taxa_names[6]);
        taxa_names_sub.push_back(terrace->matrix->taxa_names[8]);*/
        
        // for 9 taxa example
        /*taxa_names_sub.push_back(terrace->matrix->taxa_names[0]);
        taxa_names_sub.push_back(terrace->matrix->taxa_names[1]);
        taxa_names_sub.push_back(terrace->matrix->taxa_names[2]);
        taxa_names_sub.push_back(terrace->matrix->taxa_names[3]);
        taxa_names_sub.push_back(terrace->matrix->taxa_names[5]);*/
        
        
        // for 10 taxa example
        /*taxa_names_sub.push_back(terrace->matrix->taxa_names[0]);
        taxa_names_sub.push_back(terrace->matrix->taxa_names[1]);
        taxa_names_sub.push_back(terrace->matrix->taxa_names[2]);
        taxa_names_sub.push_back(terrace->matrix->taxa_names[3]);
        taxa_names_sub.push_back(terrace->matrix->taxa_names[6]);
        taxa_names_sub.push_back(terrace->matrix->taxa_names[7]);*/
        
        // for 15 taxa example cov_15_1 and tree_15_1
        /*taxa_names_sub.push_back(terrace->matrix->taxa_names[0]);
        taxa_names_sub.push_back(terrace->matrix->taxa_names[1]);
        taxa_names_sub.push_back(terrace->matrix->taxa_names[2]);
        taxa_names_sub.push_back(terrace->matrix->taxa_names[3]);
        taxa_names_sub.push_back(terrace->matrix->taxa_names[4]);
        taxa_names_sub.push_back(terrace->matrix->taxa_names[5]);
        taxa_names_sub.push_back(terrace->matrix->taxa_names[7]);
        taxa_names_sub.push_back(terrace->matrix->taxa_names[8]);
        taxa_names_sub.push_back(terrace->matrix->taxa_names[9]);*/
        
        // for 6 taxa example
        /*taxa_names_sub.push_back(terrace->matrix->taxa_names[2]);
        taxa_names_sub.push_back(terrace->matrix->taxa_names[3]);
        taxa_names_sub.push_back(terrace->matrix->taxa_names[4]);
        taxa_names_sub.push_back(terrace->matrix->taxa_names[5]);*/
        
        // for cov_7_1
        /*taxa_names_sub.push_back(terrace->matrix->taxa_names[3]);
        taxa_names_sub.push_back(terrace->matrix->taxa_names[4]);
        taxa_names_sub.push_back(terrace->matrix->taxa_names[5]);
        taxa_names_sub.push_back(terrace->matrix->taxa_names[6]);*/
        
        int init_part = 2;
        for(i=0; i<terrace->taxa_num; i++){
            if(terrace->matrix->pr_ab_matrix[i][init_part]==1){
                taxa_names_sub.push_back(terrace->matrix->taxa_names[i]);
                cout<<"TAXON "<<terrace->matrix->taxa_names[i]<<endl;
            }
        }
        
        cout<<endl;
        
        // Creating a subterrace: submatrix and an initial tree
        PresenceAbsenceMatrix *submatrix = new PresenceAbsenceMatrix();
        terrace->matrix->getSubPrAbMatrix(taxa_names_sub, submatrix);
        
        TerraceTree tree_init;
        tree_init.copyTree_byTaxonNames(terrace,taxa_names_sub);
        //tree_init.drawTree(cout, WT_BR_SCALE | WT_TAXON_ID | WT_NEWLINE);
        Terrace *init_terrace = new Terrace(tree_init, submatrix);
        
        //cout<<"taxa_names.size() = "<<init_terrace->matrix->taxa_names.size()<<endl;
        //cout<<"pr_ab_matrix.size() = "<<init_terrace->matrix->pr_ab_matrix.size()<<endl;
        //cout<<"pr_ab_matrix[i].size() = "<<init_terrace->matrix->pr_ab_matrix[0].size()<<endl;
        //init_terrace->matrix->print_pr_ab_matrix();
        
        init_terrace->linkTrees(true, false); // branch_back_map, taxon_back_map; in this case you only want to map branches
        //init_terrace->printMapInfo();
        //init_terrace-> printBackMapInfo();
        
        vector<Terrace*> part_tree_pairs;
        
        /*
        NodeVector aux_taxon_nodes;
        vector<TerraceTree*> aux_induced_part_trees;
        IntVector parts;
        bool back_branch_map = false, back_taxon_map = true;
        
        for(i=0; i<terrace->part_num; i++){
            parts.clear();
            parts.push_back(i);
            aux_taxon_nodes.clear();
            
            // 1. get submatrix for taxa, which are common for the initial tree and for the top level induced partition tree. common taxa, are the leaves of the low level induced tree (assert: they should be with 1's!!!)
            init_terrace->induced_trees[i]->getTaxa(aux_taxon_nodes);
            PresenceAbsenceMatrix *aux_submatrix = new PresenceAbsenceMatrix();
            init_terrace->matrix->getSubPrAbMatrix(aux_taxon_nodes, aux_submatrix,&parts);
            
            // 2. now update matrix with taxa, which occur only on the top level induced partition tree and not on the initial tree. Simply add 0, because it does not occur on the common subtree, i.e. low level induced partition tree
            aux_taxon_nodes.clear();
            terrace->induced_trees[i]->getTaxa(aux_taxon_nodes);
            
            parts.clear(); // auxiliary use of the variable, just need an integer vec with 0 value.
            parts.push_back(0);
            
            for(j=0; j<aux_taxon_nodes.size(); j++){
                if(aux_submatrix->findTaxonID(aux_taxon_nodes[j]->name) == -1){
                    aux_submatrix->pr_ab_matrix.push_back(parts);
                    aux_submatrix->taxa_names.push_back(aux_taxon_nodes[j]->name);
                    aux_submatrix->taxa_num +=1;
                }
            }
            
            //aux_submatrix->print_pr_ab_matrix();
            //terrace->induced_trees[i]->printTree(cout,WT_BR_LEN_ROUNDING + WT_NEWLINE);
            //cout<<endl;
            
            aux_induced_part_trees.clear();
            aux_induced_part_trees.push_back(init_terrace->induced_trees[i]);
            Terrace *aux_terrace = new Terrace(*(terrace->induced_trees[i]), aux_submatrix, aux_induced_part_trees);
            aux_terrace->linkTrees(back_branch_map, back_taxon_map);
            //aux_terrace->printMapInfo();
            //aux_terrace->printBackMapInfo();
            
            part_tree_pairs.push_back(aux_terrace);
            
        }
         */
        
        init_terrace->create_Top_Low_Part_Tree_Pairs(part_tree_pairs, terrace);
        /*for(i=0; i<part_tree_pairs.size(); i++){
            part_tree_pairs[i]->printMapInfo();
            part_tree_pairs[i]->printBackMapInfo();
        }*/
        
        
        
        /*  TODO:
         
         *  DONE: 1. get submatrix (for testing purposes, any submatrix of original presence-absence matrix)
         *  DONE: 2. get initial tree (for testing purposes, any subtree of the representative tree?)
         *  DONE: 3. get low-level induced partition trees (i.e. create a sub-terrace)
         *  DONE: -> branch and taxon maps look fine, but you need to introduce flags for separating them
         *  DONE: 4. get a vector of terraces for high- and low-level induced partition trees
         *      ABER: -> not thouroughly tested
         *  DONE: 5. perform branch mapping from subparent to low-level induced partition trees
         *  DONE: 6. perform taxon mappings from high- to low-level induced partition trees
         
         *  7. order taxa - who should be inserted first?
         *  8. get feasible branches for a taxon from all low-level induced partition trees
         
         *  WORK IN PROGRESS.... 9. insert taxon
                    TODO: at the moment there is a bit of cheating: insert taxa, clear all information about link neighbours, re-link with the new taxa, move on
         
         *  10. update all low-level induced partition trees
         *  11. update branch and taxon maps (is there a fast way, without performing a tree traversal again?)
         *  12. insert next taxon .... if no possibilities to insert, return without a tree
         *  13. if no taxa to insert -> you got a tree from a terrace. Add a tree to a list (write it to the file). Go one level up, check next possibility to insert final leaf.
         */

        
        // Get a list of taxa to be inserted.
        vector<string> list_taxa_to_insert;
        for(i=0; i<terrace->taxa_num; i++){
            if(init_terrace->matrix->findTaxonID(terrace->matrix->taxa_names[i])==-1){
               list_taxa_to_insert.push_back(terrace->matrix->taxa_names[i]);
            }
        }
        
        cout<<"Taxa to insert on initial tree: "<<endl;
        for(i=0; i<list_taxa_to_insert.size(); i++){
            cout<<i<<":"<<list_taxa_to_insert[i]<<endl;
        }
        
        cout<<endl<<"+++++++++++++++++++++++++++++++++++++++++++++++++++++++++"<<endl;
        cout<<endl<<"READY TO GENERATE TREES FROM A TERRACE"<<endl;
        cout<<endl<<"+++++++++++++++++++++++++++++++++++++++++++++++++++++++++"<<endl;
        init_terrace->print_ALL_DATA(part_tree_pairs);
        
        bool progress_status = true;
        init_terrace->generateTerraceTrees(terrace, part_tree_pairs, &list_taxa_to_insert, 0, &progress_status);
        
        
        // ==================================================================================================================================================
        
        /* insert taxon (here should be a loop afterwards, at the moment just test everything for this test example)
            - find allowed branches:
                - find node
                + FORLOOP over part
                    - get link_neighbor for part for nei_1
                    - get link_neigbour for part for nei_2
                    - get back_maps
                + get allowed branches by the overlap of branches from all partitions
        */
        /*vector<TerraceNeighbor*> nei1_vec, nei2_vec;
        
        string taxon_name;
        int id;
        
        for(i=0; i<list_taxa_to_insert.size(); i++){
            
            taxon_name = list_taxa_to_insert[i];
            // TODO: do not use FORLOOP for taxa, but a recursive function, bitte
            nei1_vec.clear();
            nei2_vec.clear();
            
            init_terrace->getAllowedBranches(taxon_name, part_tree_pairs, &nei1_vec, &nei2_vec);
            
            if(!nei1_vec.empty()){
                for(j=0; j<nei1_vec.size(); j++){
                    init_terrace->extendNewTaxon(taxon_name,(TerraceNode*)nei1_vec[j]->node,(TerraceNode*)nei2_vec[j]->node,part_tree_pairs); // insert a taxon on induced partition trees, where it occurs -> DONE in extendNewTaxon
                    id = terrace->matrix->findTaxonID(taxon_name);
                    assert(id!=-1);
                    init_terrace->matrix->extend_by_new_taxa(list_taxa_to_insert[i], terrace->matrix->pr_ab_matrix[id]);
                    
                    for(i=0; i<terrace->part_num; i++){
                        // update matrices of top-low part tree pairs
                        if(part_tree_pairs[i]->findLeafName(taxon_name)){
                            //IntVector aux_ident_vec;
                            //aux_ident_vec.push_back(1);
                            //BUG: you do not have to add a new taxon to the matrix. Just update the entry from 0 to 1.
                            //part_tree_pairs[i]->matrix->extend_by_new_taxa(taxon_name, aux_ident_vec);
                            int taxon_matrix_id = part_tree_pairs[i]->matrix->findTaxonID(taxon_name);
                            part_tree_pairs[i]->matrix->pr_ab_matrix[taxon_matrix_id][0]=1;
                        }
                        
                    }
                    
                    // re-link
                    cout<<endl<<endl<<"-----------------------------------"<<endl<<" AFTER TAXON INSERTION"<<endl<<"-----------------------------------"<<endl<<endl;
                    cout<<endl<<"PARENT TREE:"<<endl<<endl;
                    init_terrace->drawTree(cout, WT_BR_SCALE | WT_INT_NODE | WT_TAXON_ID | WT_NEWLINE);
                    init_terrace->linkTrees(true, false);
                    //init_terrace->printMapInfo();
                    //init_terrace-> printBackMapInfo();
                    for(i=0; i<terrace->part_num; i++){
                        cout<<"====================================================="<<endl;
                        cout<<endl<<"PARTITION TREE "<<i<<":"<<endl<<endl;
                        cout<<"====================================================="<<endl;
                        //if(init_terrace->induced_trees[i]->leafNum>2){
                        //    init_terrace->induced_trees[i]->drawTree(cout, WT_BR_SCALE | WT_INT_NODE | WT_TAXON_ID | WT_NEWLINE);
                        //}
                        //part_tree_pairs[i]->drawTree(cout, WT_BR_SCALE | WT_INT_NODE | WT_TAXON_ID | WT_NEWLINE);
                        part_tree_pairs[i]->linkTrees(false, true); //INFO/CHECK: I think you do not need back_taxon_map
                        part_tree_pairs[i]->printMapInfo();
                        part_tree_pairs[i]->printBackMapInfo();
                        cout<<endl;
                    }
                    
                    // repeat
                    
                    // TODO: write function to DELETE a taxon
                    
                    break; // for testing
                
                }
            } else {
                cout<<"For a given taxon "<<list_taxa_to_insert[i]<<" there are no allowed branches.. Dead end.."<<endl;
                break; //for testing, should be return in a separate function, I guess
            }
            break; // for testing
        }
        
         */
        
        // ==================================================================================================================================================
        
        
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

