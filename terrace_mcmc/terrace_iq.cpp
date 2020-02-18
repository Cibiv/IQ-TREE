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
    representative_tree = tree;
    n_part = ((SuperAlignment*)aln)->partitions.size();
    getALLInducedPartitionTrees(tree, aln, induced_part_trees);
    
    // set taxa_num, taxa_names, and pr_ab_matrix
    initTaxaNamesPresenceAbsenceM();
    
}

Terrace_IQ::Terrace_IQ(Alignment* alignment, vector<MTree*> trees){
    
    aln = alignment;
    representative_tree = nullptr;
    n_part = ((SuperAlignment*)aln)->partitions.size();
    vector<MTree*>::iterator it;
    
    for(it = trees.begin(); it!=trees.end(); it++)
        induced_part_trees.push_back((*it));
    
    // set taxa_num, taxa_names, and pr_ab_matrix
    initTaxaNamesPresenceAbsenceM();
    
}

Terrace_IQ::Terrace_IQ(const char* file_presence_absence, MTree* tree){
    
    aln = nullptr;
    representative_tree = tree;
    // The function below intinilises n_part, taxa_num, taxa_names, pr_ab_matrix
    readPresenceAbsenceMatrix(file_presence_absence);
    getALLInducedPartitionTreesM();
}

Terrace_IQ::Terrace_IQ(vector<IntVector> matrix, vector<string> names, MTree* tree){
    // you should be careful with odering of rows and entries in matrix and names!!!
    // the first row in matrix corresponds to the first taxon in names
    aln = nullptr;
    representative_tree = tree;
    n_part=matrix[0].size();
    taxa_num=matrix.size();
    for(int i=0; i<taxa_num; i++){
        pr_ab_matrix.push_back(matrix[i]);
        taxa_names.push_back(names[i]);
    }
    getALLInducedPartitionTreesM();
}

void Terrace_IQ::readPresenceAbsenceMatrix(const char *infile){
    ifstream in;
    //cout<<endl<<"-----------------------------------------------------"<<endl;
    try {
        in.exceptions(ios::failbit | ios::badbit);
        in.open(infile);
        in.exceptions(ios::badbit);
        readPresenceAbsenceMatrix(in);
        in.close();
    } catch (const char* str) {
        outError(str);
    } catch (ios::failure) {
        outError(ERR_READ_INPUT, infile);
    }
}

void Terrace_IQ::readPresenceAbsenceMatrix(istream &in) {
    
    string str_rest,name;
    
    if(!(in>>taxa_num)) throw "The first line should start with a number of taxa followed by the number of partitions!";
    if(!(in>>n_part)) throw "The first line should start with a number of taxa followed by the number of partitions!";
    //getline(in,str_rest);
    
    int i=0,j=0;
    for(i=0; i<taxa_num; i++){
        if(!(in>>name)) throw "Each line should start with a taxon name!";
        taxa_names.push_back(name);

        IntVector vec(n_part, -1);
        for(j=0; j<n_part; j++){
            if(!(in >> vec[j])) throw "Could not read a matrix entry! For each species make sure there are as many entries as the number of partitions specified in the first line of the file. Moreover, presence-absence matrix should only contain 0, 1!";
            if(vec[j] < 0) throw "Error: A negative entry! Presence-absence matrix should only contain 0, 1!";
            if(vec[j] > 1) throw "Error: The entry is greater than 1! Presence-absence matrix should only contain 0, 1!";
        }
        pr_ab_matrix.push_back(vec);
    }
}


void Terrace_IQ::setInducedPartitionTrees(vector<MTree*> trees){

    induced_part_trees.clear();
    
    vector<MTree*>::iterator it;
    for(it = trees.begin(); it!=trees.end(); it++)
        induced_part_trees.push_back((*it));
}

void Terrace_IQ::initTaxaNamesPresenceAbsenceM(){
    
    taxa_num = ((SuperAlignment*)aln)->getNSeq();
    for(int id=0; id<taxa_num; id++){
        IntVector vec(n_part, 0);
        for(int part=0; part<n_part; part++){
            if(((SuperAlignment*)aln)->taxa_index[id][part]>=0){
                vec[part]=1;
            }
        }
        pr_ab_matrix.push_back(vec);
        taxa_names.push_back(aln->getSeqName(id));
    }
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

void Terrace_IQ::printInducedTrees(ostream &out){
    vector <MTree*>::iterator it;
    for(it = induced_part_trees.begin(); it != induced_part_trees.end(); it++){
        (*it)->printTree(out,WT_SORT_TAXA | WT_BR_LEN_ROUNDING | WT_NEWLINE);
    }
}

Terrace_IQ::~Terrace_IQ(){
    for (vector<MTree*>::reverse_iterator it = induced_part_trees.rbegin(); it != induced_part_trees.rend(); it++) {
        MTree *tree = *it;
        delete tree;
    }
    // why don't you delete an aln pointer?
}

int Terrace_IQ::getTaxonID_in_pr_ab_m(string taxon_name){
    
    for(int i=0; i<taxa_num; i++){
        if(taxa_names[i]==taxon_name){
            return i;
        }
    }
    
    return -1;
}

void Terrace_IQ::getALLInducedPartitionTreesM(){
    
    //cout<<"GETTING ALL INDUCED PARTITION TREES for tree: "<<endl;
    //cout<<getTreeTopologyString(tree)<<endl;
    //cout<<"induced trees:"<<endl;
    
    int id;
    string taxa_set = "";
    string taxon_name;
    NodeVector taxa_nodes;
    NodeVector::iterator it2;
    vector<uint32_t> check_int;
    check_int.resize(taxa_num);
    
    
    assert(representative_tree != nullptr && "Sorry, a terrace does not have a representative tree, I cannot get induced partition trees using pr_ab_matrix!");
    representative_tree->getTaxa(taxa_nodes);
    
    for(int part=0; part<n_part; part++){
        MTree* induced_tree = new MTree();
        for(it2=taxa_nodes.begin();it2!=taxa_nodes.end();it2++){
            taxon_name=(*it2)->name;
            //id=this->getTaxonID_in_pr_ab_m(taxon_name);
            // you need to check the consistency between taxa nodes on a terrace and in pr_ab_matrix, do this when you read pr_ab_matrix from a file!
            assert(id != -1 && "Not all of the taxa appear in the pr_ab_matrix!");
            check_int[(*it2)->id] = pr_ab_matrix[id][part];
        }
        
        taxa_set.insert(taxa_set.begin(), check_int.begin(), check_int.end());
        induced_tree->copyTree(representative_tree,taxa_set);
        induced_part_trees.push_back(induced_tree);
        //induced_part_trees[part]->printTree(cout,WT_BR_LEN_ROUNDING | WT_NEWLINE);
    }
}
