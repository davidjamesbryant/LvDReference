/**
 * @File decompositionLileihood.cpp
 * @brief Routines for computing updating likelihoods decomposition trees.
 * @author David Bryant & Celine Scornavacca
 * @version 0.1
 */

#include "decompositionLikelihood.h"


/**
 Compute the tree likelihoods.
 
 It assumes that the ID of the leaf nodes in the tree corresponds to the
 corresponding row in seqs. This version has no compression.
 **/

Scalar  computeLikelihood(phylo<DecomNodeData>& decomTree, const SubstModel& model, const vector<sequence>& seqs, vector<Scalar>& siteL, Stopwatch& timer) {
    
    timer.stop();
    
    unsigned long nsites = seqs[1].size();
    siteL.resize(nsites);

    vector<phylo<DecomNodeData>::iterator> taxaPointers;

        
    Scalar logL = 0.0;
    timer.start();
    
      
    for(unsigned long  s = 0;s<nsites;s++) {
        for (auto p = decomTree.leftmost_leaf();!p.null();p=p.next_post()) {
            
                
            if (p.leaf() && p->isClade) {
                
//TODO: We could have the model remember the last value of t and related values.
                
                //cout<<"Leaf "<<p->id<<endl;
                base b = seqs[p->id][s];
                if (b<model.num_states()) {
                    //Unambiguous character
                    for (int i=0;i<model.num_states();i++)
                        p->partialVec(i) = model.Pij(i,b,p->length);
                } else if (b==bX) {
                    p->partialVec.fill(1.0);
//                    for (int i=0;i<model.num_states();i++)
//                        p->partialVec(i) = 1.0;
                } else {  //Ambiguous character
                    vector<Scalar> bases;
                    resolve_base(b,bases);
                    for (int i=0;i<model.num_states();i++) {
                        Scalar pb = 0.0;
                        for (int j=0;j<model.num_states();j++)
                            if (bases[j]==1)
                                pb+=+model.Pij(i,j,p->length);
                        p->partialVec(i) = pb;
                    }
                }
                p->exponent = 0;
            } else if (p.leaf() && !p->isClade) {
                //cout<<"Internal edge "<<endl;
                for (int i=0;i<model.num_states();i++)
                    for (int j=0;j<model.num_states();j++)
                        p->partialMat(i,j) = model.Pij(i,j,p->length);
                p->exponent = 0;
            } else {
                
                auto l = p.left();
                auto r = l.right();
                
                switch(p->mergeType){
                    case 1:
                        p->partialVec = l->partialVec.array() * r -> partialVec.array(); // .* multiplication on Eigen
                        break;
                    case 2:
                        //right child is the segment
                        for(int i=0;i<4;i++) {
                            p->partialMat(i,Eigen::all) = l->partialVec(i) * r->partialMat(i,Eigen::all);
                        }
                        break;
                    case 3:
                        for(int j=0;j<4;j++)
                            p->partialMat(Eigen::all,j) = l->partialMat(Eigen::all,j) * r->partialVec(j);
                        break;
                    case 4:
                        p->partialVec = l->partialMat * r->partialVec;
                        break;
                    case 5:
                        p->partialMat.noalias() = l->partialMat * r->partialMat;
                        break;
                    default:
                        cerr<<"Invalid merge index"<<endl;
                }
                
                //Deal with underflow
                p->exponent = l->exponent + r->exponent;
                
                Scalar maxval = 0;
                if (p->isClade)
                    maxval = p->partialVec.maxCoeff();
                else
                    maxval = p->partialMat.maxCoeff();
                
                while (0<maxval && maxval < UNDERFLOW_CUTOFF) {
                    maxval *= UNDERFLOW_MULTIPLIER;
                    if (p->isClade)
                        p->partialVec *=UNDERFLOW_MULTIPLIER;
                    else
                        p->partialMat *= UNDERFLOW_MULTIPLIER;
                    p->exponent-=UNDERFLOW_STEP;
                }
                
                
    
            }
        }
        
        Scalar Ls = 0.0;
      for(int i=0;i<model.num_states();i++)
          Ls += model.pi(i)*(decomTree.root()->partialVec(i));
      
      siteL[s] = log(Ls) + decomTree.root()->exponent;

      logL += siteL[s];
    
    }
    timer.stop();
    
    return logL;
}




/**
 Set up arrays for partial likelihoods and fill out transitions probabilities for the nodes in the decomposition tree corresponding to internal edges
 */
static void initialisePartials(phylo<DecomNodeData>& decomTree, int ntaxa, vector<phylo<DecomNodeData>::iterator>& taxaPointers, const SubstModel& model) {
    taxaPointers.resize(ntaxa);
    for(int i=0;i<ntaxa;i++)
        taxaPointers[i].set_null();

    for(auto p = decomTree.leftmost_leaf(); !p.null(); p=p.next_post()) {
        p->isDirty = false;
        if (p->isClade)
            p->partialVec.fill(0.0);
        else
            p->partialMat.fill(0.0);
        
        if (p.left().null()) {
            if (p->isClade) {
                taxaPointers[p->id] = p;
            } else {
//TODO add a class to model to return transition matrix.
                
                for(int i=0;i<4;i++)
                    for(int j=0;j<4;j++)
                        p->partialMat(i,j) = model.Pij(i,j,p->length);
            }
        }
    }
}

static void updatePartialsRecurse(phylo<DecomNodeData>::iterator p, phylo<DecomNodeData>& decomTree) {
    
    

        //Check that any dirty children have been updated
    for (auto q = p.left(); !q.null(); q = q.right()) {
        if (q->isDirty)
            updatePartialsRecurse(q,decomTree);
    }
    //Update partials on this node, depending on the merger type
    //Note - if this node is dirty then it has children.
    auto l = p.left();
    auto r = l.right();
    
    
    switch(p->mergeType){
        case 1:
            p->partialVec = l->partialVec.array() * r -> partialVec.array(); // .*
//            for(int i=0;i<4;i++)
//                p->partials[i][0] = l->partials[i][0] * r -> partials[i][0];
            break;
        case 2:
            for(int i=0;i<4;i++) {
                p->partialMat(i,Eigen::all) = l->partialVec(i) * r->partialMat(i,Eigen::all);
            }
            //right child is the segment
//            for(int i=0;i<4;i++)
//                for(int j=0;j<4;j++)
//                    p->partials[i][j] = l->partials[i][0] * r->partials[i][j];
            break;
        case 3:
            for(int j=0;j<4;j++)
                p->partialMat(Eigen::all,j) = l->partialMat(Eigen::all,j) * r->partialVec(j);

//            for(int i=0;i<4;i++)
//                for(int j=0;j<4;j++)
//                    p->partials[i][j] = l->partials[i][j] * r->partials[j][0];
            break;
        case 4:
            p->partialVec = l->partialMat * r->partialVec;
//            for(int i=0;i<4;i++) {
//                p->partials[i][0] = 0;
//                for(int j=0;j<4; j++)
//                    p->partials[i][0] += l->partials[i][j] * r->partials[j][0];
//            }
            break;
        case 5:
            //mul4x4_strassen_on_strassen_wrap(l->partialMat,r->partialMat,p->partialMat);
            p->partialMat.noalias() = l->partialMat * r->partialMat;
//            for(int i=0;i<4;i++)
//                for(int j=0;j<4; j++) {
//                    p->partialMat(i,j) = 0;
//                    Scalar pij = 0;
//                    for(int k=0;k<4; k++) {
//                        pij += l->partialMat(i,k)*r->partialMat(k,j);
//                    }
//                    p->partialMat(i,j) = pij;
//                }
            break;
        default:
            cerr<<"Invalid merge index"<<endl;
    }
    
    //Deal with underflow
    p->exponent = l->exponent + r->exponent;
    
    Scalar maxval = 0;
    if (p->isClade)
        maxval = p->partialVec.maxCoeff();
    else
        maxval = p->partialMat.maxCoeff();
    
    while (maxval > 0 && maxval < UNDERFLOW_CUTOFF) {
        maxval *= UNDERFLOW_MULTIPLIER;
        if (p->isClade)
            p->partialVec *=UNDERFLOW_MULTIPLIER;
        else
            p->partialMat *= UNDERFLOW_MULTIPLIER;
        p->exponent-=UNDERFLOW_STEP;
    }

    
    

    p->isDirty = false;
    
   /*  cout<<"End of update Partial "<<p->length<<" "<<p->id<<endl;
        debugPrintPartials(cout,decomTree);
    cout<<endl;
 */
}
/**
 Update the partial likelihoods for all the nodes ancestral to
 */
Scalar updatePartials(phylo<DecomNodeData>& decomTree, vector<phylo<DecomNodeData>::iterator>& taxaPointers, const SubstModel& model, const vector<pair<int, unsigned short>>& diffs) {
    //First mark the nodes which need updating
    for(unsigned int i=0;i<diffs.size();i++) {
        auto p = taxaPointers[diffs[i].first];
        if (p.null())
            continue;
        base b = diffs[i].second;
        if (b<model.num_states()) {
            //Unambiguous character
            for (int i=0;i<model.num_states();i++)
                p->partialVec(i) = model.Pij(i,b,p->length);
        } else if (b==bX) {
            p->partialVec.fill(1.0);
           //for (int i=0;i<model.num_states();i++)
                //p->partials[i][0] = 1.0;
        } else {
            vector<Scalar> bases;
            resolve_base(b,bases);
            for (int i=0;i<model.num_states();i++) {
                Scalar pb = 0.0;
                for (int j=0;j<model.num_states();j++)
                    if (bases[j]==1)
                        pb+=+model.Pij(i,j,p->length);
                p->partialVec(i) = pb;
            }
        }
        p->exponent = 0;
        p->isDirty = false; //This node is now updated
        p=p.par();
        
        //Mark vertices above - we can stop as soon as we meet one which is already marked
        while (p!=decomTree.header() && !p->isDirty) {
            p->isDirty=true;
            p=p.par();
        }
    }
    /* cout<<"Checking dirty nodes\n";
    debugPrintPartials(cout, decomTree);
    cout<<endl<<endl; */

    updatePartialsRecurse(decomTree.root(),decomTree);
        
    Scalar Ls = 0.0;
    for(int i=0;i<model.num_states();i++)
        Ls += model.pi(i)*(decomTree.root()->partialVec(i));
      
    return log(Ls) + decomTree.root()->exponent;
    
}


Scalar computeLikelihoodUsingUpdating(phylo<DecomNodeData>& decomTree, const SubstModel& model, const vector<sequence>& seqs, PatternSorter patternSorter, Stopwatch& timer) {
    
    timer.stop();
    
    vector<pair<Pattern, int>> patterns = compressPatterns(seqs,patternSorter);
    auto differences = computePatternDifferences(patterns);
    vector<phylo<DecomNodeData>::iterator> taxaPointers;
    initialisePartials(decomTree,seqs.size(),taxaPointers,model);

    
    timer.start();
    
    Scalar logL = 0.0;
    
    for (unsigned int i=0;i<patterns.size();i++) {
        Scalar patternlogL = updatePartials(decomTree,taxaPointers,model,differences[i]);
        logL += patternlogL * patterns[i].second;
       // if (i<5)
     //       cout<<"Pattern "<<i<<" logL = "<<patternlogL<<" freq = "<<patterns[i].second<<" exp = "<<decomTree.root()->exponent<<endl;
    }
    timer.stop();
    
    return logL;
}

void debugPrintPartials(ostream& os, phylo<DecomNodeData>& decomTree) {
    os<<"Partial likelihoods\n";
    for(auto p = decomTree.leftmost_leaf();!p.null();p=p.next_post()) {
        if (p->isClade)
            os<<"["<<p->partialVec<<"]";
        else
            os << "["<<p->partialMat<<"]";
        if (p.leaf() && p->isClade) {
            os<<" pendant edge, id = "<<p->id<<" length = "<<p->length;
        } else if (p.leaf() && !p->isClade) {
            os<<" internal edge, length = "<<p->length;
        }
        else {
            os<<" mergeType = "<<p->mergeType;
        }
        if (p->isDirty)
            os<<" dirty";
        else
            os<<" clean";
        os<<endl;
    }
}

/**
 Use the decomposition likelihood algorithm without updating to compute the likelihood for all the patterns.
 */
Scalar computeLikelihood(phylo<DecomNodeData>& decomTree, const SubstModel& model, const vector<pair<Pattern, int>>& patterns, vector<double>& patternL, Stopwatch& timer) {
    
    timer.stop();
    
    unsigned long npatterns = patterns.size();
    patternL.resize(npatterns);
    
    Scalar logL = 0.0;
    timer.start();

    for(unsigned long  s = 0;s<npatterns;s++) {
        Pattern site = patterns[s].first;
        
        for (auto p = decomTree.leftmost_leaf();!p.null();p=p.next_post()) {
            
                
            if (p.leaf() && p->isClade) {
                //cout<<"Leaf "<<p->id<<endl;
                base b = site[p->id];
                if (b<model.num_states()) {
                    //Unambiguous character
                    for (int i=0;i<model.num_states();i++)
                        p->partialVec(i) = model.Pij(i,b,p->length);
                } else if (b==bX) {
                    p->partialVec.fill(1.0);
                } else {  //Ambiguous character
                    vector<Scalar> bases;
                    resolve_base(b,bases);
                    for (int i=0;i<model.num_states();i++) {
                        Scalar pb = 0.0;
                        for (int j=0;j<model.num_states();j++)
                            if (bases[j]==1)
                                pb+=+model.Pij(i,j,p->length);
                        p->partialVec(i) = pb;
                    }
                }
                p->exponent = 0;
            } else if (p.leaf() && !p->isClade) {
                //cout<<"Internal edge "<<endl;
                for (int i=0;i<model.num_states();i++)
                    for (int j=0;j<model.num_states();j++)
                        p->partialMat(i,j) = model.Pij(i,j,p->length);
                p->exponent = 0;
            } else {
                
                auto l = p.left();
                auto r = l.right();
          
                switch(p->mergeType){
                    case 1:
                        p->partialVec = l->partialVec.array() * r -> partialVec.array();
                        break;
                    case 2:
                        //right child is the segment
                        for(int i=0;i<4;i++) {
                            p->partialMat(i,Eigen::all) = l->partialVec(i) * r->partialMat(i,Eigen::all);
                        }
                        break;
                    case 3:
                        for(int j=0;j<4;j++)
                            p->partialMat(Eigen::all,j) = l->partialMat(Eigen::all,j) * r->partialVec(j);
                        break;
                    case 4:
                        p->partialVec = l->partialMat * r->partialVec;
                        break;
                    case 5:
                        p->partialMat.noalias() = l->partialMat * r->partialMat;
                        break;
                    default:
                        cerr<<"Invalid merge index"<<endl;
                }
                
                //Deal with underflow
                p->exponent = l->exponent + r->exponent;
                
                Scalar maxval = 0;
                if (p->isClade)
                    maxval = p->partialVec.maxCoeff();
                else
                    maxval = p->partialMat.maxCoeff();
                
                while (maxval>0 && maxval < UNDERFLOW_CUTOFF) {
                    maxval *= UNDERFLOW_MULTIPLIER;
                    if (p->isClade)
                        p->partialVec *=UNDERFLOW_MULTIPLIER;
                    else
                        p->partialMat *= UNDERFLOW_MULTIPLIER;
                    p->exponent-=UNDERFLOW_STEP;
                }
                
                
    
            }
            /* //DEBUG
            cout<<"[";
            for (int i=0;i<4;i++)  {
                //DEBUG
                cout<<"[";
                for (int j=0;j<(int)p->partials[i].size();j++)
                    cout<<p->partials[i][j]<<" ";
                cout<<"]"<<endl;
            }
            cout<<"]"<<endl; */
        }
        
        Scalar Ls = 0.0;
      for(int i=0;i<model.num_states();i++)
          Ls += model.pi(i)*(decomTree.root()->partialVec(i));
      
      
        patternL[s] = log(Ls) + decomTree.root()->exponent;
     // if (s<5)
       //   cout<<s<<"\t"<<siteL[s]<<"\t"<<decomTree.root()->exponent<<endl;

      //debugPrintPartials(cout, decomTree);

      logL += patternL[s]*patterns[s].second;
    
    }
    timer.stop();
    

    return logL;
}
    
    
    
    


Scalar computeLikelihood(phylo<DecomNodeData>& decomTree, const SubstModel& model, const vector<pair<Pattern, int>>& patterns, const PatternDiffs& differences, vector<double>& patternL, Stopwatch& timer) {
    
    cerr<<"com[puting L"<<endl;
    
    timer.stop();
    int ntax = patterns[0].first.size();
    vector<phylo<DecomNodeData>::iterator> taxaPointers;
    initialisePartials(decomTree,ntax,taxaPointers,model);
    patternL.resize(patterns.size());
    
    timer.start();
    
    Scalar logL = 0.0;
    
    for (unsigned int i=0;i<patterns.size();i++) {
        cerr<<i<<endl;
        Scalar patternlogL = updatePartials(decomTree,taxaPointers,model,differences[i]);
        patternL[i] = patternlogL;
        logL += patternlogL * patterns[i].second;
    }
    timer.stop();
    
    return logL;
}

    
    
    


