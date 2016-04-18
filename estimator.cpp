//
//  estimator.cpp
//  iqtree
//
//  Created by Olga on 21/03/16.
//
//

#include "estimator.h"
#include <string.h>

void estimatorAnalysis(Params* params, Alignment* alignment, IQTree* tree){
    cout<<"Starting estimator analysis...."<<endl;
// Output details ----------------------------------------------------------------------------------
// Pattern Info File
    string out_file = params->out_prefix;
    bool append = FALSE;
    
    if(append){
        out_file = "results.PatternInfo";   // Write results from multiple tree+alignment into one file // not a good idea, since you don't know in advance how many patterns do you have for a specific alignment
    }else{
        //out_file += ".PatternInfo";
        out_file += ".modALN";
    }
//--------------------------------------------------------------------------------------------------
    
    Alignment *aln_modified = new Alignment();
    aln_modified->createAlignmentPatternsOnly(alignment);
    aln_modified->printPhylip(out_file.c_str());
    
    
    
    //printPatternLhFreq(out_file.c_str(), (PhyloTree*)tree, NULL,append);
}




void printPatternLhFreq(const char*filename, PhyloTree *tree, double *ptn_lh,
                        bool append, const char *linename) {
    cout<<"Getting pattern information...."<<endl;
        int i;
        double *pattern_lh;
        if (!ptn_lh) {
            pattern_lh = new double[tree->getAlnNPattern()];
            tree->computeLikelihood();
            tree->computePatternLikelihood(pattern_lh);
        } else
            pattern_lh = ptn_lh;
        try {
            ofstream out;
            out.exceptions(ios::failbit | ios::badbit);
            if (append) {
                out.open(filename, ios::out | ios::app);
            } else {
                out.open(filename);
            }
            IntVector freq;
            tree->aln->getPatternFreq(freq);
            if (linename){
                out.width(10);
                out << left << linename;
            }
            int alnLEN = tree->getAlnNSite();
            int taxaNUM = tree->leafNum;
            for (i = 0; i < tree->getAlnNPattern(); i++){
                out<<" "<< taxaNUM <<" "<< alnLEN <<" "<<tree->getAlnNPattern()<<" "<<exp(pattern_lh[i])<<" "<<freq[i]<<" "<<(double)freq[i]/(double)alnLEN << endl;;
            }
            out.close();
            cout << "Pattern info was printed to " << filename << endl;
        } catch (ios::failure) {
            outError(ERR_WRITE_OUTPUT, filename);
        }
        
        if (!ptn_lh)
            delete[] pattern_lh;
}