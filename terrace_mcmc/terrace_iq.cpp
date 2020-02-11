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

Terrace_IQ::Terrace_IQ(const char* file_presence_absence, MTree* tree){
    
    aln = nullptr;
    
    
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
    
    int taxon_num;
    
    if(!(in>>taxon_num)) throw "The first line should start with a number of taxa followed by the number of partitions!";
    
    
    
    string name;
    if(!(in>>name)) throw "Each line should start with a taxon name!";
    this->taxa_names.push_back(name);
    
    
    
    
    
    
    int i=0,j=0;
    if(!(in>>SpeciesNUM)) throw "The first line must contain the number of species in this Food Web!!";
    string str_rest, speciesName;
    getline(in, str_rest);

    
    vector<double*> MM;
    for(i=0;i<SpeciesNUM;i++){
        MM.push_back(new double [SpeciesNUM]);
    }
    
    nvar = (TaxaNUM > SpeciesNUM) ? TaxaNUM : SpeciesNUM;
    for(i=0;i<nvar;i++){
        DAG.push_back(new double [nvar]);
        for(j=0; j<nvar; j++){
            DAG[i][j] = 0.0;
        }
    }
    i = 0;
    j = 0;
    if(rooted){
        while(i != SpeciesNUM-1){
            if(!(in >> speciesName)) throw "Each line should start with a species name!";
            dagNames.push_back(speciesName);
            j = 0;
            while(j != SpeciesNUM-1){
                if(!(in >> MM[i][j])) throw "Could not read matrix entry! For each species make sure there are as many entries as the number of species specified in the file. Only square matrices are accepted.";
                if(MM[i][j] < 0) throw "The Food Web matrix should not contain negative values.Use either 0, 1 or a positive number to indicate the portion of diet.";
                j++;
            }
            MM[i][SpeciesNUM-1] = 0;
            i++;
        }
        for(j=0; j<SpeciesNUM; j++)
            MM[SpeciesNUM-1][j] = 0;
        dagNames.push_back("_root");
    } else {
        while(i != SpeciesNUM){
            if(!(in >> speciesName)) throw "Each line should start with a species name!";
            dagNames.push_back(speciesName);
            j = 0;
            while(j != SpeciesNUM){
                if(!(in >> MM[i][j])) throw "Could not read matrix entry! For each species make sure there are as many entries as the number of species specified in the file. Only square matrices are accepted.";
                if(MM[i][j] < 0) throw "The Food Web matrix should not contain negative values.Use either 0, 1 or a positive number to indicate the portion of diet.";
                j++;
            }
            i++;
        }
    }
    
    /* ---------------------------------------------------------------------------------------------------------
     * Input data
     * ---------------------------------------------------------------------------------------------------------*/
    if(verbose_mode == VB_MAX){
        cout<<endl<<"Food web is defined by the following matrix"<<endl;
        for(i=0;i<SpeciesNUM;i++) {
            cout<<dagNames[i]<<"\t";
            for(j=0;j<SpeciesNUM;j++)
                cout<<MM[i][j]<<"\t";
            cout<<endl;
        }
        // Species in the food web and their ids
        for(i=0; i<SpeciesNUM;i++)
            cout<<"["<<i<<"] "<<dagNames[i]<<endl;
    }
    /* ---------------------------------------------------------------------------------------------------------
     * Processing the input data
     * ---------------------------------------------------------------------------------------------------------*/
    //Ignoring cannibalism -------------------------------------------------------------------------
    int cannibals=0;
    for(i=0;i<SpeciesNUM;i++)
        if(MM[i][i]!=0){
            cannibals++;
            if(weighted){
                if(cannibals == 1){
                    cout<<"------------------------------------------"<<endl;
                    cout<<"    Cannibal species         link weight  "<<endl;
                    cout<<"------------------------------------------"<<endl;
                }
                cout.width(30);
                cout<<left<<dagNames[i];
                cout<<" | "<<MM[i][i]<<endl;
            }else{
                if(cannibals == 1){
                    cout<<"-----------------------------"<<endl;
                    cout<<"       Cannibal species      "<<endl;
                    cout<<"-----------------------------"<<endl;
                }
                cout<<dagNames[i]<<endl;
            }
            MM[i][i]=0;
        }
    if(cannibals!=0){
        cout<<endl<<"Deleted "<<cannibals;
        if(cannibals == 1)
            cout<<" cannibalistic link."<<endl;
        else
            cout<<" cannibalistic links."<<endl;
    }
    
    //Check whether the graph is acyclic or not
    Graph g(SpeciesNUM);
    for(i=0; i<SpeciesNUM; i++)
        for(j=0; j<SpeciesNUM; j++)
            if(MM[i][j]>0)
                g.addEdge(i,j);
    if(g.isCyclic()){
        if(cannibals != 0)
            cout<<endl<<"ERROR: Even after deleting cannibalistic links, there are still some cycles present."<<endl;
        else
            cout<<endl<<"ERROR: ";
        cout<<"Cyclic food webs are not supported. Delete the links which cause cycles and run again."<<endl;
        cout<<"SOLUTION:"<<endl;
        cout<<"Detect species in the cycle and choose one link to be deleted in order to break the cycle."<<endl;
        cout<<"One possibility is to delete the link with least weight. This can be done by setting the corresponding value in the matrix to 0."<<endl;
        exit(0);
    }
    
    // The number of links -------------------------------------------------------------------------
    linksNUM = 0;
    for(i = 0; i<SpeciesNUM; i++)
        for(j = 0; j<SpeciesNUM; j++)
            if(MM[i][j]>0)
                linksNUM++;
    
    //Rescaling the diet if necessary --------------------------------------------------------------
    if(weighted){
        int dietReScaled = 0;
        vector<double> colsum;
        //cout<<"Food web is weighted."<<endl;
        for(j=0;j<SpeciesNUM;j++){
            
            colsum.push_back(0);
            for(i=0;i<SpeciesNUM;i++)
                colsum[j]=colsum[j]+MM[i][j];
            if(colsum[j]!=1 && colsum[j]!=0){
                dietReScaled++;
                //cout<<"    WARNING: rescaled diet composition of species "<<j<<". Column sum = "<<colsum[j]<<endl;
                for(i=0;i<SpeciesNUM;i++)
                    MM[i][j]=MM[i][j]/colsum[j];
            }
            colsum[j]=0;
            //for(i=0;i<SpeciesNUM;i++)
            //    colsum[j]=colsum[j]+MM[i][j];
            //cout<<j<<"  Column sum = "<<colsum[j]<<endl;
        }
        cout<<"Rescaled diet composition of "<<dietReScaled<<" species."<<endl;
    }else{
        for(i=0; i<SpeciesNUM; i++)
            for(j=0; j<SpeciesNUM; j++)
                if( MM[i][j] > 0)
                    MM[i][j] = 1;
        //cout<<"Since the -eco option was chosen, the entries of Food Web matrix will be converted to 0/1 [not prey / prey]. You can use -ecoW option to account for the Diet Composition."<<endl;
    }
    
    // Technical: in case of rooted trees, we check which species are basal ones, i.e. for which check = 0, and set them to "feed on" root M[i,j] = 1
    if(rooted){
        vector<double> check;
        for(j=0;j<SpeciesNUM-1;j++){
            check.push_back(0);
            for(i=0;i<SpeciesNUM-1;i++)
                check[j]=check[j]+MM[i][j];
            if(check[j]==0)
                MM[SpeciesNUM-1][j]=1;
        }
    }
    
    //Detecting which species are not present in either FoodWeb or Tree/SplitNetwork-----------------
    detectMissingSpecies();
    
    //Check whether all the species from initialTaxa set are actually present on Tree/SplitSys or in Food Web
    checkInitialTaxa();
    
    // Synchronization of species in Tree/SplitSys and species in FoodWeb ---------------------------
    synchronizeSpecies();
    
    for(i=0; i<SpeciesNUM; i++){
        for(j=0; j<SpeciesNUM; j++){
            DAG[phylo_order[i]][phylo_order[j]]=MM[i][j];
        }
    }
    
    for(i=SpeciesNUM-1;i>=0;i--)
        delete[] MM[i];
    
    if(verbose_mode == VB_MAX){
        // Print info about synchronization
        cout<<endl<<"Synchronization:"<<endl;
        cout<<"PhyloInfo id | FoodWeb id, name"<<endl;
        for(i=0; i<SpeciesNUM; i++){
            cout<<"["<<phylo_order[i]<<"] | ["<<i<<"] "<<dagNames[i]<<endl;
        }
        cout<<"PhyloInfo id | name"<<endl;
        for(i=0; i<TaxaNUM;i++){
            cout<<"["<<i<<"] "<<findNodeID(i)->name<<endl;
        }
        
        // Input data after processing: cannibalism, rescaling, reordering
        cout<<endl<<"Food web is defined by the following matrix"<<endl;
        for(i=0;i<nvar;i++) {
            if(findFoodWebID(i) != -1)
                cout<<dagNames[findFoodWebID(i)]<<"\t";
            else
                cout<<"\t\t";
            for(j=0;j<nvar;j++)
                cout<<DAG[i][j]<<"\t";
            cout<<endl;
        }
    }
    /* ---------------------------------------------------------------------------------------------------------
     * Filling out taxaDAG vector: node corresponds to taxa, neighbors to preys, length (node-neighbor) to weight
     * ---------------------------------------------------------------------------------------------------------*/
    vector<int> vec2;//the value of vec[j] is the height of the species in the DAG
    taxaDAG.resize(nvar,NULL);
    for(j=0;j<nvar;j++){
        taxaDAG[j] = newNode(j,j);
        //cout<<"taxonDAG[j="<<j+1<<"]->id="<<taxaDAG[j]->id<<endl;
    }
    
    for(j=0;j<nvar;j++){
        for(i=0;i<nvar;i++)
            if(DAG[i][j]>0){
                //cout<<"cheking matrix"<<i<<j<<endl;
                taxaDAG[j]->addNeighbor(taxaDAG[i], DAG[i][j], taxaDAG[i]->id);
                //cout<<"neighbors[i="<<taxaDAG[j]->degree()-1<<"]->id="<<taxaDAG[j]->neighbors[taxaDAG[j]->degree()-1]->node->id<<endl;
            }
        //cout<<endl;
    }
    
    /* ---------------------------------------------------------------------------------------------------------
     * Defining levels in the Food Web based on the longest food chain of predators
     * ---------------------------------------------------------------------------------------------------------*/
    for(j=0;j<nvar;j++){
        levelDAG.push_back(0);
        if(taxaDAG[j]->degree()>0)
            vec2.push_back(1);
        else
            vec2.push_back(0);
        //          if(taxaDAG[j]->degree()>0){
        //              cout<<"Children of taxonDAG[j="<<j<<"]->id="<<taxaDAG[j]->id<<":"<<endl;
        //             for(i=0;i<taxaDAG[j]->degree();i++)
        //                 cout<<"taxaDAG["<<j<<"]->neighbors["<<i<<"]->node->id "<<taxaDAG[j]->neighbors[i]->node->id<<endl;
        //                 //cout<<"id of the child "<<i<<" node id "<<taxaDAG[j]->neighbors[i]->node->id+1<<" "<<endl;
        //                 //cout<<"           neighbors[i="<<i<<"]->id="<<taxaDAG[j]->neighbors[i]->node->id<<endl;
        //             cout<<endl;
        //
        //          }
    }
    //    for(j=0;j<nvar;j++)
    //        cout<<j<<" "<<levelDAG[j]<<" "<<vec2[j]<<endl;
    
    int eq=0,step=0;
    //cout<<"Starting while..."<<endl;
    while(eq!=1){
        eq=1;
        step++;
        //         if(step==1 or step==2 or step==3)
        //        cout<<"-------STEP "<<step<<"-------"<<endl<<"j v1 v2"<<endl;
        for(j=0;j<nvar;j++){
            if(levelDAG[j]!=vec2[j])
                eq=0;
            //             if(step==1 or step==2 or step==3)
            //            cout<<j<<" "<<levelDAG[j]<<" "<<vec2[j]<<endl;
            levelDAG[j]=vec2[j];
        }
        for(j=0;j<nvar;j++){
            if(taxaDAG[j]->degree()>0){
                //cout<<"taxaDAG["<<j<<"]->neighbors[0]->node->id "<<taxaDAG[j]->neighbors[0]->node->id<<endl;
                vec2[j]=vec2[taxaDAG[j]->neighbors[0]->node->id]+1;
                for(i=1;i<taxaDAG[j]->degree();i++)
                    if(vec2[taxaDAG[j]->neighbors[i]->node->id]>=vec2[j])
                        vec2[j]=vec2[taxaDAG[j]->neighbors[i]->node->id]+1;
            }
        }
    }
    
    // For each predator the level corresponds to its longest food chain----------------------------
    if(verbose_mode == VB_MAX){
        cout<<"For each species its longest chain according to a food web"<<endl;
        for(j=0;j<nvar;j++)
            //if(findFoodWebID(j) != -1)
            //    cout<<dagNames[findFoodWebID(j)]<<"\t| "<<levelDAG[j]<<endl;
            //else
            cout<<*names[j]<<"\t| "<<levelDAG[j]<<endl;
    }
    //cout<<"Species - level"<<endl;
    //ofstream levelF;
    //levelF.open("Level",ios::app);
    //for(j=0;j<SpeciesNUM;j++)
    //    levelF<<j+1<<" "<<levelDAG[j]<<endl;
    //     for(i=0;i<tree.leafNum;i++)
    //     myfile<<"taxon id: "<<taxaTree[i]->id<<" | taxon name: "<<taxaTree[i]->name<<endl;
    //     myfile<<"root  id: "<<root->id<<" | root  name: "<<root->name<<endl;
    // // myfile.close();
    
    // The maximum level is the longest food chain of the food web ---------------------------------
    //     int maxlevel;
    //     maxlevel=0;
    //     for(i=0;i<SpeciesNUM;i++)
    //        if(maxlevel<levelDAG[i])
    //            maxlevel=levelDAG[i];
    
    // Decrease SpeciesNUM since you do not need to include the root to the Species anymore---------
    if(rooted)
        SpeciesNUM--;
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
