//
//  terrace_iq.cpp
//  iqtree
//
//  Created by Olga on 07.11.19.
//

#include "terrace_iq.hpp"
#include "terracecheck.hpp"

Terrace_IQ::Terrace_IQ(){
    
}

Terrace_IQ::Terrace_IQ(Alignment* alignment, MTree* tree){
    
    aln = alignment;
    n_part = ((SuperAlignment*)aln)->partitions.size();
    getALLInducedPartitionTrees(tree, aln, induced_part_trees);
    
}

Terrace_IQ::Terrace_IQ(Alignment* alignment, vector<MTree*> trees){
    
    aln = alignment;
    n_part = ((SuperAlignment*)aln)->partitions.size();
    vector<MTree*>::iterator it;
    
    for(it = trees.begin(); it!=trees.end(); it++)
        induced_part_trees.push_back((*it));
    
}

void Terrace_IQ::setInducedPartitionTrees(vector<MTree*> trees){

    induced_part_trees.clear();
    
    vector<MTree*>::iterator it;
    for(it = trees.begin(); it!=trees.end(); it++)
        induced_part_trees.push_back((*it));
}

bool Terrace_IQ::checkTree(MTree* tree){
    
    vector<MTree*> induced_trees_query;
    getALLInducedPartitionTrees(tree, aln, induced_trees_query);
    
    vector<string> trees, taxa_names;
    
    for(int part=0; part<n_part; part++){
        
        trees.clear();
        trees.push_back(getTreeTopologyString(induced_part_trees[part]));
        //trees.push_back(getTreeTopologyString(induced_trees_query[part]));
        
        taxa_names.clear();
        induced_part_trees[part]->getTaxaName(taxa_names);
        
        MTreeSet tree_set;
        tree_set.init(trees,taxa_names,induced_part_trees[part]->rooted);
        //tree_set[0]->printTree(cout,WT_SORT_TAXA|WT_BR_LEN_ROUNDING|WT_NEWLINE);
        
        MTreeSet tree_set_1;
        trees.clear();
        trees.push_back(getTreeTopologyString(induced_trees_query[part]));
        tree_set_1.init(trees,taxa_names,induced_part_trees[part]->rooted);
        //tree_set_1[0]->printTree(cout,WT_SORT_TAXA|WT_BR_LEN_ROUNDING|WT_NEWLINE);
        
        int *rfdist;
        int n=1;
        rfdist = new int [n*n];
        memset(rfdist, 0, n*n* sizeof(int));
        tree_set.computeRFDist(rfdist,&tree_set_1);

        //cout<<"Partition "<<part+1<<": RF-distance between induced trees is "<<rfdist[0]<<endl;
        
        if(rfdist[0]!=0){
            cout<<"---------------------------------------------"<<endl;
            cout<<"A query tree does not belong to the terrace:"<<endl;
            //tree->printTree(cout,WT_SORT_TAXA|WT_BR_LEN_ROUNDING|WT_NEWLINE);
            
            cout<<"The first encountered pair of induced partition trees that do not match. Partition "<<part+1<<endl;
            tree_set[0]->printTree(cout,WT_SORT_TAXA|WT_BR_LEN_ROUNDING|WT_NEWLINE);
            tree_set_1[0]->printTree(cout,WT_SORT_TAXA|WT_BR_LEN_ROUNDING|WT_NEWLINE);

            tree_set.clear();
            
            delete [] rfdist;
            return false;
        }
        delete [] rfdist;
        tree_set.clear();
        
    }
    cout<<"---------------------------------------------"<<endl;
    cout<<"A query tree is ON the terrace"<<endl;

    return true;
}

Terrace_IQ::~Terrace_IQ(){
    for (vector<MTree*>::reverse_iterator it = induced_part_trees.rbegin(); it != induced_part_trees.rend(); it++) {
        MTree *tree = *it;
        delete tree;
    }
}
