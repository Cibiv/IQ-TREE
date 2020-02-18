//
//  terracecheck.cpp
//  iqtree
//
//  Created by Olga on 05.11.19.
//

#include "terracecheck.hpp"
#include "terrace_iq.hpp"
#include "tree/mtreeset.h"
#include "alignment/superalignment.h"
#include "tree/iqtree.h"
#include "tree/phylosupertree.h"


void runTerraceCheck(Params &params){
    
    string out_terrace = params.out_prefix;
    out_terrace += ".terrace_ON";
    
    string out_not_terrace = params.out_prefix;
    out_not_terrace += ".terrace_OUT";
    
    //string out_induced = params.terrace_rep_file;
    //out_induced += ".induced_trees";
    
    ofstream out_1,out_2,out_3;
    out_1.open(out_terrace.c_str());
    out_2.open(out_not_terrace.c_str());
    //out_3.open(out_induced.c_str());

    Alignment *aln;
    aln = new SuperAlignment(params);
    
    /*
     * Suggestion:
     * - function to check if one tree belongs to the terrace with tree_1 (representative tree)
     *      - create collections of MTreeSets for induced partition trees
     *      - use computeRFdist on these sets to check if the induced partition trees are pairwise identical
     *      - for each input tree from the set, output bool if belongs to the terrace
     *      - output a separate file with all trees that belong to the same terrace as tree_1
     *      - optionally, compute pairwise rf distances for these trees to form a terrace subgraph
     */
    
    MTree tree(params.terrace_rep_file, params.is_rooted);
    if(tree.rooted){
        cout<<"This feature is only available for unrooted trees.. Converning to unrooted tree.."<<endl;
        exit(0);
        //((PhyloTree*) &tree)->convertToUnrooted();
    }
    //ASSERT(checkTaxonSetsAlnTree(aln,&tree) && "Not all of the taxa appear in the alignment!");
    cout<<"=============================================="<<endl;
    cout<<"Terrace is represented by the following tree"<<endl;
    tree.printTree(cout,WT_SORT_TAXA | WT_BR_LEN_ROUNDING | WT_NEWLINE);
    cout<<"=============================================="<<endl;
    out_1<<getTreeTopologyString(&tree)<<endl;
    Terrace_IQ terrace(aln,&tree);
    //terrace.printInducedTrees(out_3);
    
    cout<<"Corresponding induced partition trees:"<<endl;
    for(vector<MTree*>::iterator it2 = terrace.induced_part_trees.begin(); it2!=terrace.induced_part_trees.end();it2++){
        (*it2)->printTree(cout,WT_BR_LEN_ROUNDING | WT_NEWLINE);
    }
    //cout<<"=============================================="<<endl;
    
    /* Get a set of trees to be tested, whether on the same terrace as tree */
    MTreeSet tree_set_query(params.terrace_query_set, params.is_rooted, params.tree_burnin, params.tree_max_count);
    MTreeSet::iterator it;
    tree_set_query.checkConsistency();

    int nseq = aln->getNSeq();
    for(it=tree_set_query.begin(); it!=tree_set_query.end(); it++){
    if((*it)->rooted){
        cout<<"This feature is only available for unrooted trees.. Converning to unrooted tree.."<<endl;
        exit(0);
         convertToUnrooted_MTree(*it, nseq);
        //((PhyloTree*) &tree)->convertToUnrooted();
    }
    }
    
    /* Check if the tree is on the terrace */
    int i = 0;
    for(it = tree_set_query.begin(); it != tree_set_query.end(); it++){
        i+=1;
        cout<<"=============================================="<<endl;
        cout<<"Tree "<<i<<" from the query set: "<<endl;
        //ASSERT(checkTaxonSetsAlnTree(aln,(*it)) && "Not all of the taxa appear in the alignment!");
        if(terrace.checkTree((*it))){
            out_1<<getTreeTopologyString((*it))<<endl;
        } else {
            out_2<<getTreeTopologyString((*it))<<endl;
        }
    }
    
    out_1.close();
    out_2.close();
    //out_3.close();
    
}

string getTreeTopologyString(MTree* tree){
        stringstream tree_stream;
        tree->printTree(tree_stream, WT_BR_LEN_ROUNDING + WT_SORT_TAXA);
        return tree_stream.str();
}

MTree* getInducedPartitionTree(MTree* tree_complete, Alignment* aln, int part){
    
    MTree* induced_tree = new MTree();
    
    /*
     *  Get taxon_set (0-1) string, 1's for taxa to remain on the tree. The order must be the same as on the tree.
     */
    
    int n_taxa = tree_complete->getNumTaxa();
    vector<uint32_t> check_int;
    check_int.resize(n_taxa);
    
    string taxa_set = "";
    string taxon_name;
    
    NodeVector taxa_nodes;
    NodeVector::iterator it2;
    
    tree_complete->getTaxa(taxa_nodes);
    
    for(it2=taxa_nodes.begin();it2!=taxa_nodes.end();it2++){
        taxon_name=(*it2)->name;
        assert(aln->getSeqID(taxon_name) != -1 && "Not all of the taxa appear in the alignment!");
        check_int[(*it2)->id] = (((SuperAlignment*)aln)->taxa_index[aln->getSeqID(taxon_name)][part] >= 0)? 1 : 0;
        //cout<<"PART: "<<part<<" | TAXON ID: "<<(*it2)->id<<" | TAXON NAME: " <<taxon_name<<" | TAXON INDEX: "<<((SuperAlignment*)aln)->taxa_index[aln->getSeqID(taxon_name)][part]<<endl;
    }
    
    taxa_set.insert(taxa_set.begin(), check_int.begin(), check_int.end());
    induced_tree->copyTree(tree_complete,taxa_set);
    
    //tree_complete->printTree(cout,WT_BR_LEN_ROUNDING | WT_NEWLINE);
    //induced_tree->printTree(cout,WT_BR_LEN_ROUNDING | WT_NEWLINE);
    
    return induced_tree;
}



void getALLInducedPartitionTrees(MTree* tree, Alignment* aln,vector<MTree*> &induced_part_trees){
    
    int n_part = ((SuperAlignment*)aln)->partitions.size();
    
    //cout<<"GETTING ALL INDUCED PARTITION TREES for tree: "<<endl;
    //cout<<getTreeTopologyString(tree)<<endl;
    //cout<<"induced trees:"<<endl;
    for(int part=0; part<n_part; part++){
        induced_part_trees.push_back(getInducedPartitionTree(tree, aln, part));
        //induced_part_trees[part]->printTree(cout,WT_BR_LEN_ROUNDING | WT_NEWLINE);
    }
}

int terrace_check_two_trees(vector<MTree*> &induced_part_trees_rep, MTree* tree_query){

    int result = -1;     // 0 - not on the terrace; 1 - is a part of a terrace;
    return result;
}

bool checkTaxonSetsAlnTree(Alignment* aln, MTree* tree){
    string taxon_name;
    
    NodeVector taxa_nodes;
    NodeVector::iterator it2;
    
    tree->getTaxa(taxa_nodes);
    for(it2=taxa_nodes.begin();it2!=taxa_nodes.end();it2++){
        taxon_name=(*it2)->name;
        if(aln->getSeqID(taxon_name) == -1){
            tree->printTree(cout,WT_BR_LEN_ROUNDING | WT_NEWLINE | WT_SORT_TAXA);
            cout<<"Taxon "<<aln->getSeqID(taxon_name)<<" does not appear in the alignment! Exiting.."<<endl;
            return false;
        }
    }
    
    return true;
}

void convertToUnrooted_MTree(MTree* tree,int nseq){
    //ASSERT(rooted);
    //if (aln)
     //   ASSERT(leafNum == aln->getNSeq()+1);
    ASSERT(tree->root);
    ASSERT(tree->root->isLeaf());
    //&& tree->root->id == tree->leafNum-1
    Node *node = tree->root->neighbors[0]->node;
    Node *taxon = tree->findFirstTaxon();
    
    tree->rooted = false;
    tree->leafNum--;
    
    // delete root node
    if (node->degree() == 3) {
        // delete and join adjacent branches
        Node *node1 = NULL, *node2 = NULL;
        double len = 0.0;
        FOR_NEIGHBOR_IT(node, tree->root, it) {
            if (!node1) node1 = (*it)->node; else node2 = (*it)->node;
            len += (*it)->length;
        }
        node1->updateNeighbor(node, node2, len);
        node2->updateNeighbor(node, node1, len);
        delete node;
    } else {
        // only delete root node
        auto it = node->findNeighborIt(tree->root);
        delete *it;
        node->neighbors.erase(it);
        
    }
    
    delete tree->root;
    // set a temporary taxon so that tree traversal works
    tree->root = taxon;
    
    //if (params)
     //   setRootNode(params->root);
    
    tree->initializeTree();
    //    computeBranchDirection();
}

void runPrintInduced(Params &params){
    string out_terrace = params.out_prefix;
    out_terrace += ".induced_trees";

    ofstream out;
    out.open(out_terrace.c_str());
    
    Alignment *aln;
    aln = new SuperAlignment(params);
    
    MTree tree(params.terrace_print_induced_trees, params.is_rooted);
    
    vector<MTree*> induced_trees;
    
    getALLInducedPartitionTrees(&tree, aln, induced_trees);
    
    vector<MTree*>::iterator it;
    for(it = induced_trees.begin(); it != induced_trees.end(); it++){
        (*it)->printTree(out,WT_SORT_TAXA | WT_BR_LEN_ROUNDING | WT_NEWLINE);
    }
    
    out.close();
}

vector<IntVector> getSubMatrix(vector<IntVector> pr_ab_complete, vector<string> taxa_names, MTree* tree){
    int id;
    NodeVector taxa_nodes;
    NodeVector::iterator it;
    string taxon_name;
    tree->getTaxa(taxa_nodes);
    for(it=taxa_nodes.begin(); it<taxa_nodes.end(); it++){
        taxon_name=(*it)->name;
        //id=getTaxonID_in_pr_ab_m(taxon_name);
    }
}



