//
//  phyloTerraceSampling.cpp
//  iqtree
//
//  Created by Olga on 18.10.19.
//

#include "phyloTerraceSampling.h"
#include "alignment/alignment.h"
#include "alignment/superalignment.h"
#include "tree/iqtree.h"
#include "tree/phylotree.h"
#include "tree/phylosupertree.h"
#include "terrace.h"


void runTerraceSampling(Params &params){
    Alignment *alignment;
    IQTree *tree;
    
    alignment = new SuperAlignment(params);
    tree = new PhyloSuperTree((SuperAlignment*)alignment);
    
    ((PhyloTree*)tree)->readTreeFile(params.user_file);
    
    //((PhyloSuperTree*) tree)->mapTrees();

    stringstream tree_stream;
    tree->printTree(tree_stream, WT_SORT_TAXA);
    
    cout<<tree_stream.str()<<endl;
    
    
    //Terrace* terrace = tree->getTreeTerrace(tree_stream.str());
    vector<string> terraceSample;
    terraceSample = tree->getTerraceSample(params.terrace_sample_size,params.terrace_burnin);
    
/*terrace->burnin = params.terrace_burnin;
    
    if(terrace->isTrivial()){
        cout<<"A terrace is a trivial one or has disconnected components"<<endl;
    } else {
        cout<<"There are NNI-neighbours on the terrace. Starting the sampling"<<endl;
        vector<string> sample;
        sample = terrace->getSample(params.terrace_sample_size);
    
        vector<string>::iterator it;
        for(it=sample.begin();it<sample.end();it++){
            cout<<(*it)<<endl;
        }
    }
 */
}
