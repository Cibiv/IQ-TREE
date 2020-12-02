//
//  presenceabsencematrix.cpp
//  iqtree
//
//  Created by Olga on 10.09.20.
//

#include "presenceabsencematrix.hpp"
#include "alignment/superalignment.h"
#include "tree/node.h"
#include "tree/mtree.h"

void PresenceAbsenceMatrix::get_from_alignment(Params &params){
    int i, j;
    
    // Input is an alignment with partition info
    SuperAlignment* alignment = new SuperAlignment(params);
    IntVector pattern_0_1;
    
    taxa_num = alignment->getNSeq();
    part_num = alignment->getNSite();
    
    for(i=0; i<taxa_num; i++){
        pattern_0_1.clear();
        for(j=0; j<part_num; j++){
            if(alignment->taxa_index[i][j]==-1){
                pattern_0_1.push_back(0);
            }else{
                pattern_0_1.push_back(1);
            }
        }
        pr_ab_matrix.push_back(pattern_0_1);
        taxa_names.push_back(alignment->getSeqName(i));
    }

    init();
}

void PresenceAbsenceMatrix::read_pr_ab_matrix(const char *infile){
    ifstream in;
    //cout<<endl<<"-----------------------------------------------------"<<endl;
    try {
        in.exceptions(ios::failbit | ios::badbit);
        in.open(infile);
        in.exceptions(ios::badbit);
        read_pr_ab_matrix(in);
        in.close();
    } catch (const char* str) {
        outError(str);
    } catch (ios::failure) {
        outError(ERR_READ_INPUT, infile);
    }
}

void PresenceAbsenceMatrix::read_pr_ab_matrix(istream &in){
    string str_rest,name;
    
    if(!(in>>taxa_num)) throw "The first line should start with a number of taxa followed by the number of partitions!";
    if(!(in>>part_num)) throw "The first line should start with a number of taxa followed by the number of partitions!";
    
    int i=0,j=0;
    for(i=0; i<taxa_num; i++){
        if(!(in>>name)) throw "Each line should start with a taxon name!";
        if(name == "0" or name == "1") throw "Each line should start with a taxon name! 0 and 1 are not allowed as taxon names.";
        //cout<<"Taxon "<<i<<": "<<name<<endl;
        taxa_names.push_back(name);
        vector<int> vec(part_num, -1);
        for(j=0; j<part_num; j++){
            if(!(in >> vec[j])) throw "Could not read a matrix entry! For each species make sure there are as many entries as the number of partitions specified in the first line of the file. Moreover, presence-absence matrix should only contain 0, 1!";
            if(vec[j] < 0) throw "Error: A negative entry! Presence-absence matrix should only contain 0, 1!";
            if(vec[j] > 1) throw "Error: The entry is greater than 1! Presence-absence matrix should only contain 0, 1!";
        }
        pr_ab_matrix.push_back(vec);
    }
    
    init();
};

void PresenceAbsenceMatrix::print_pr_ab_matrix(){
    
    cout<<"Presence-absence matrix:"<<endl;
    for(int i=0; i<taxa_num; i++){
        cout<<taxa_names[i]<<" ";
        for(int j=0; j<part_num; j++){
            cout<<pr_ab_matrix[i][j]<<" ";
        }
        cout<<endl;
    }
    cout<<endl;
};


int PresenceAbsenceMatrix::findTaxonID(string taxon_name){
    int id;
    for(id=0; id<taxa_names.size(); id++){
        if(taxa_names[id]==taxon_name){
            //cout<<"GET_TAXON_ID: taxa_name = "<<taxon_name<<" | taxa_id = "<<id<<endl;
            return id;
        }
    }
    return -1;
};

void PresenceAbsenceMatrix::init(){
    flag_reorderAccordingToTree = false;
    missing_percent = 0.0;
}

void PresenceAbsenceMatrix::getPartTaxa(int part, MTree *tree, MTree *part_tree, NodeVector &part_taxa){
    
    NodeVector taxa_nodes;
    tree->getTaxa(taxa_nodes);
    
    part_taxa.resize(taxa_num);
    
    //if(!flag_reorderAccordingToTree){
    //    reorderAccordingToTree(taxa_nodes);
    //}
    
    Node *node;
    int taxon_matrix_id;
    for(NodeVector::iterator it=taxa_nodes.begin(); it<taxa_nodes.end(); it++){
        //cout<<(*it)->name<<" id = "<<(*it)->id<<endl;
        taxon_matrix_id = findTaxonID((*it)->name);
        //cout<<"TAXON_MATRIX_ID:"<<taxon_matrix_id<<endl;
        if(pr_ab_matrix[taxon_matrix_id][part] == 1){
            node = part_tree->findLeafName((*it)->name);
            //cout<<"PREPARING PART_TAXA: part = "<<part<<"|leaf_id = "<<(*it)->id<<"|leaf_name = "<<(*it)->name<<endl;
            assert(node && "ERROR: The leaf is not found on partition tree!");
            part_taxa[taxon_matrix_id]=node;
        }
    }
    
    bool check = false;
    if(check){
        //cout<<"GET_taxa_nodes_info...."<<endl;
        cout<<"Presence-absence info:"<<endl;
        for(int i=0; i<taxa_num; i++){
            cout<<pr_ab_matrix[i][part]<<" ";
        }
        cout<<endl;
        cout<<"Taxon names info:"<<endl;
        for(int i=0; i<taxa_num; i++){
            cout<<taxa_names[i]<<" ";
        }
        cout<<endl;
        
        cout<<"Partition taxa info:"<<endl;
        for(int i=0; i<taxa_num; i++){
            if(part_taxa[i]){
                cout<<part_taxa[i]->name<<"("<<part_taxa[i]->id<<") ";
            }else{
                cout<<"NA ";
            }
        }
        cout<<endl;
    }
}

void PresenceAbsenceMatrix::reorderAccordingToTree(NodeVector taxa_nodes){

    // WARNING: when adding new taxa, this function is not helpful, because the ids of new taxa (at the current setting, as of 06.10.20) are larger than the number of taxa (id of new taxon is set to the number of nodes, which is then increased by 1, when a taxon is added)
    
    //cout<<"BEFORE reordering according to the tree:"<<endl;
    //print_pr_ab_matrix();
    
    int i=0;
    int id=-1;
    
    vector<IntVector> aux_matrix;
    vector<string> aux_names;
    
    aux_matrix.resize(taxa_num);
    aux_names.resize(taxa_num);
    
    for(i=0; i<taxa_nodes.size();i++){
        id=findTaxonID(taxa_nodes[i]->name);
        //cout<<taxa_nodes[i]->name<<" id = "<<id<<endl;
        aux_matrix[taxa_nodes[i]->id]=pr_ab_matrix[id];
        aux_names[taxa_nodes[i]->id]=taxa_names[id];
    }
    
    pr_ab_matrix.clear();
    taxa_names.clear();
    
    for(i=0; i<taxa_num; i++){
        pr_ab_matrix.push_back(aux_matrix[i]);
        taxa_names.push_back(aux_names[i]);
    }
    
    //cout<<"AFTER reordering according to the tree:"<<endl;
    //print_pr_ab_matrix();
}

vector<IntVector> getSubMatrix(vector<IntVector> pr_ab_complete, vector<string> taxa_names, MTree* tree){
    
    //int id;
    vector<IntVector> sub_matrix;
    NodeVector taxa_nodes;
    NodeVector::iterator it;
    string taxon_name;
    tree->getTaxa(taxa_nodes);
    for(it=taxa_nodes.begin(); it<taxa_nodes.end(); it++){
        taxon_name=(*it)->name;
        //id=getTaxonID_in_pr_ab_m(taxon_name);
    }
    
    return sub_matrix;
};

void PresenceAbsenceMatrix::getSubPrAbMatrix(vector<string> taxa_names_subset, PresenceAbsenceMatrix *submatrix, IntVector *parts){
 
    vector<string> not_found_taxon_names;
    bool found = false;
    
    int i,j,h;
    for(i=0; i<taxa_names_subset.size(); i++){
        //cout<<i<<": subset_taxon = "<<taxa_names_subset[i]<<endl;
        found = false;
        for(j=0; j<taxa_names.size(); j++){
            //cout<<j<<": taxon_name = "<<taxa_names[j]<<endl;
            if(taxa_names_subset[i]==taxa_names[j]){
                //cout<<"MATCH"<<endl;
                if(parts){
                    IntVector taxon_coverage;
                    for(IntVector::iterator k=parts->begin(); k<parts->end(); k++){
                        h = (*k);
                        taxon_coverage.push_back(pr_ab_matrix[j][h]);
                    }
                    submatrix->pr_ab_matrix.push_back(taxon_coverage);
                }else{
                    submatrix->pr_ab_matrix.push_back(pr_ab_matrix[j]);
                }
                submatrix->taxa_names.push_back(taxa_names[j]);
                found = true;
                break;
            }
        }
        if(!found){
            cout<<"Taxon "<<taxa_names_subset[i]<<" is not found in the presence-absence matrix..."<<endl;
            not_found_taxon_names.push_back(taxa_names_subset[i]);
        }
    }
    
    bool print_info = false;
    if(not_found_taxon_names.size()<taxa_names_subset.size()){
        submatrix->taxa_num = submatrix->taxa_names.size();
        submatrix->part_num = submatrix->pr_ab_matrix[0].size();
        if(print_info){
            cout<<"INFO: original matrix."<<endl;
            print_pr_ab_matrix();
            cout<<endl<<"INFO: a submatrix for "<<submatrix->taxa_num<<" taxa was extracted."<<endl;
            submatrix->print_pr_ab_matrix();
        }
    }
}

void PresenceAbsenceMatrix::getSubPrAbMatrix(NodeVector taxon_nodes, PresenceAbsenceMatrix *submatrix, IntVector *parts){

    vector<string> taxon_names;
    for(NodeVector::iterator it = taxon_nodes.begin(); it<taxon_nodes.end(); it++){
        taxon_names.push_back((*it)->name);
    }
    
    getSubPrAbMatrix(taxon_names,submatrix,parts);
}

void PresenceAbsenceMatrix::extend_by_new_taxa(string taxon_name, IntVector pr_ab_pattern){
    
    taxa_names.push_back(taxon_name);
    pr_ab_matrix.push_back(pr_ab_pattern);
    
    taxa_num+=1;
    
    flag_reorderAccordingToTree = false;
}

void PresenceAbsenceMatrix::remove_taxon(string taxon_name){
    
    //cout<<"REMOVING taxon "<<taxon_name<<" from matrix."<<endl;
    int id = findTaxonID(taxon_name);
    
    /*cout<<"BEFORE:"<<endl;
    print_pr_ab_matrix();
    cout<<"MATRIX   dim: "<<pr_ab_matrix.size()<<"x"<<pr_ab_matrix[0].size()<<endl;
    cout<<"TAXAname dim:"<<taxa_names.size()<<endl;*/
    
    pr_ab_matrix.erase(pr_ab_matrix.begin()+id);
    taxa_names.erase(taxa_names.begin()+id);
    
    taxa_num-=1;
    
    /*cout<<"AFTER:"<<endl;
    print_pr_ab_matrix();
    cout<<"MATRIX   dim: "<<pr_ab_matrix.size()<<"x"<<pr_ab_matrix[0].size()<<endl;
    cout<<"TAXAname dim:"<<taxa_names.size()<<endl;*/
    
    flag_reorderAccordingToTree = false;
    
}

void PresenceAbsenceMatrix::getINFO_init_tree_taxon_order(vector<string> &taxa_names_sub, vector<string> &list_taxa_to_insert){
    
    
    // TODO!!!!
    // starting tree is the most important stuff
    // one again think about the best choice
    // - the one with the largest overlap or the largest?
    // - include or not to include unique taxa on the initial tree?
    
    
    //cout<<endl<<"=================================================="<<endl<<"INITIAL tree and Taxon Order:"<<endl<<endl;
    int i,j,k;
    
    IntVector part_cov, taxon_cov, uniq_taxa;
    part_cov.resize(part_num,0);
    taxon_cov.resize(taxa_num,0);
    uniq_taxa.resize(part_num,0);
    
    for(i=0; i<part_num;i++){
        for(j=0; j<taxa_num; j++){
            taxon_cov[j]+=pr_ab_matrix[j][i];
        }
    }
    
    for(i=0; i<part_num; i++){
        for(j=0; j<taxa_num; j++){
            if(taxon_cov[j]!=1){ // taxa covered only by one partition
                part_cov[i]+=pr_ab_matrix[j][i]; // per partition you only sum up taxa, which are present in more than one partition
            }
            if(taxon_cov[j]==1 && pr_ab_matrix[j][i]==1){
                uniq_taxa[i]+=1;
            }
        }
    }
    
    // order partitions by taxon coverage
    vector<int> ordered_partitions;
    int id;
    ordered_partitions.push_back(0);
    
    bool inserted;
    for(i=1; i<part_num; i++){
        inserted = FALSE;
        for(j=0; j<ordered_partitions.size(); j++){ // TODO: maybe implement something more efficient
            id=ordered_partitions[j];
            //if(part_cov[i]+uniq_taxa[i]>=part_cov[id]+uniq_taxa[id]){
              if(part_cov[i]>=part_cov[id]){
                ordered_partitions.insert(ordered_partitions.begin()+j, i);
                inserted = TRUE;
                break;
            }
        }
        if(!inserted){
            ordered_partitions.push_back(i);
        }
    }
    
    //cout<<"Partition with maximum coverage based on non-spesific taxa is "<<ordered_partitions[0]<<" with "<<part_cov[ordered_partitions[0]]<<" taxa."<<endl;
    
    // Vector for part id stores its order
    vector<int> part_order;
    part_order.resize(part_num);
    for(i=0; i<ordered_partitions.size(); i++){
        part_order[ordered_partitions[i]]=i;
    }
    
    /*for(i=0; i<part_num; i++){
        cout<<"Partition "<<i<<": ordered as "<<part_order[i]<<" with coverage of "<<part_cov[i]<<endl;
    }*/
    
    int part_max = ordered_partitions[0];
    
    vector<IntVector> ordered_taxa_ids;
    ordered_taxa_ids.resize(part_num+1);
    
    // collecting taxon names from partition with the largest subtree excluding unique taxa
    cout<<"Partition "<<part_max+1<<" is chosen for the initial tree."<<endl;
    for(j=0; j<taxa_num; j++){
        //if(pr_ab_matrix[j][part_max]==1 && taxon_cov[j]>1){
        if(pr_ab_matrix[j][part_max]==1){
            taxa_names_sub.push_back(taxa_names[j]);
        } else {
            // ordering other taxa by the number of partitions they are covered by
            for(i=1;i<part_num+1;i++){
                if(taxon_cov[j]==i){
                    ordered_taxa_ids[i].push_back(j);
                    break;
                }
            }
        }
    }
    
    IntVector ordered_ids;
    
    for(i=1; i<part_num; i++){
        if(ordered_taxa_ids[i].size()>1){
            //cout<<endl<<"ANALYSING COVERAGE SIZE "<<i<<endl;
            vector<IntVector> cov_within_cat;
            cov_within_cat.resize(ordered_taxa_ids[i].size());
            ordered_ids.clear();
            
            for(j=0; j<ordered_taxa_ids[i].size(); j++){
                for(k=0;k<part_num;k++){
                    if(pr_ab_matrix[ordered_taxa_ids[i][j]][k]==1){
                        cov_within_cat[j].push_back(part_order[k]);
                    }
                }
                sort(cov_within_cat[j].begin(), cov_within_cat[j].end());
                
                //cout<<"TAXON_COVERAGE_INFO: taxon_"<<ordered_taxa_ids[i][j]<<"|"<<taxa_names[ordered_taxa_ids[i][j]];
                //for(k=0; k<cov_within_cat[j].size(); k++){
                //    cout<<" "<<cov_within_cat[j][k];
                //}
                //cout<<endl;
                
                
            }
            
            // REORDER taxa within each coverage category
            orderTaxaByCoverage(ordered_taxa_ids[i],cov_within_cat,ordered_ids);
            
            ordered_taxa_ids[i].clear();
            ordered_taxa_ids[i] = ordered_ids;
            
            /*cout<<"REORDERED taxa:"<<endl;
            for(k=0; k<ordered_taxa_ids[i].size(); k++){
                cout<<k+1<<"|taxon "<<ordered_taxa_ids[i][k]<<"|"<<taxa_names[ordered_taxa_ids[i][k]]<<endl;
            }*/
            
        }
    }
    
    for(i=ordered_taxa_ids.size()-1; i>-1; i--){
        if(!ordered_taxa_ids[i].empty()){
            //cout<<"Taxa with "<<i<<" gene coverage: "<<endl;
            for(j=0; j<ordered_taxa_ids[i].size(); j++){
                //cout<<" "<<j<<": "<<ordered_taxa_by_coverage[i][j]<<endl;
                list_taxa_to_insert.push_back(taxa_names[ordered_taxa_ids[i][j]]);
            }
            //cout<<endl;
        }
    }
    
    cout<<endl<<"FINAL taxon order:"<<endl;
    for(j=0; j<list_taxa_to_insert.size(); j++){
        cout<<j<<": "<<list_taxa_to_insert[j]<<endl;
    }
    
    //cout<<endl<<"=================================================="<<endl;
}


void PresenceAbsenceMatrix::orderTaxaByCoverage(vector<int> &taxon_ids, vector<IntVector> &coverage_info, IntVector &ordered_taxa){
    
    int i,j,k;
    
    IntVector ordered_taxa_ids;
    ordered_taxa_ids.push_back(0);
    bool inserted_check;
    
    /*for(i=0; i<coverage_info.size(); i++){
        for(j=0; j<coverage_info[i].size(); j++){
            cout<<" "<<coverage_info[i][j];
        }
        cout<<endl;
    }*/
    
    for(i=1; i<taxon_ids.size(); i++){
        //cout<<"CHECKING taxon "<<i<<" with id "<<taxon_ids[i]<<endl;
        inserted_check = FALSE;
        
        for(j=0; j<ordered_taxa_ids.size(); j++){
            //cout<<"comparing with taxon "<<j<<" with id "<<taxon_ids[ordered_taxa_ids[j]]<<endl;
            int ident=0;
            for(k=0; k<coverage_info[0].size(); k++){
                //cout<<"partition_order_info:"<<coverage_info[i][k]<<" vs. "<<coverage_info[ordered_taxa_ids[j]][k]<<endl;
                if(coverage_info[i][k]<coverage_info[ordered_taxa_ids[j]][k]){
                    // insert
                    ordered_taxa_ids.insert(ordered_taxa_ids.begin()+j,i);
                    inserted_check = TRUE;
                    //cout<<"inserted"<<endl;
                    break;
                } else if(coverage_info[i][k]>coverage_info[ordered_taxa_ids[j]][k]){
                    // do not insert, check next taxon
                    //cout<<"breaking, check next taxon.."<<endl;
                    break;
                } else {
                    ident+=1;
                    // TODO: maybe do something more elaborate with taxa having identical pattern of partition coverage, such that you do not have to check next taxon with n taxa with the same partition order..
                }
            }
            
            if(inserted_check){
                break;
            }
        }
        
        // if the taxon is last one, insert to the end of the vector
        if(!inserted_check){
            ordered_taxa_ids.push_back(i);
        }
    }
    
    for(i=0; i<taxon_ids.size(); i++){
        ordered_taxa.push_back(taxon_ids[ordered_taxa_ids[i]]);
    }
    
}

void PresenceAbsenceMatrix::percent_missing(){
    
    missing_percent = 0.0;
    int i,j,a = 0;
    for(i=0; i<taxa_num; i++){
        for(j=0; j<part_num; j++){
            a+= pr_ab_matrix[i][j];
        }
    }

    double total = taxa_num*part_num;
    missing_percent=(total-a)/total*100.0;
}











