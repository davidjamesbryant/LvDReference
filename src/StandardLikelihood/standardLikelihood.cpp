/**
 * @File standardLikelihood.cpp
 * @brief Routines for calculating likelihoods using standard algorithm
 * @author David Bryant & Celine Scornavacca
 * @version 0.1
 */

#include "standardLikelihood.h"

/**
 Compute the tree likelihoods.
 
 It assumes that the ID of the leaf nodes in the tree corresponds to the
 corresponding row in seqs. This version has no compression.
 **/


//TODO: Check the Underflow correction for overflow... 



Scalar  computeLikelihood(phylo<nodeData>& tree, const SubstModel& model, const vector<sequence>& seqs, vector<Scalar>& siteL, Stopwatch& timer) {
    
    timer.stop(); //We don't want to include the preprocessing time.
    
    unsigned long nsites = seqs[1].size();

    siteL.resize(nsites);

    typedef typename phylo<nodeData>::iterator ITERATOR;
    Scalar logL = 0.0;
    
    timer.start();
    
    for(unsigned long  s = 0;s<nsites;s++) {
        for(ITERATOR p = tree.leftmost_leaf();!p.null();p=p.next_post()) {
            
            p->partials = {0,0,0,0};
            if (p.leaf()) {
                base b = seqs[p->id][s];
                resolve_base(b,p->partials);
                p->exponent = 0;
            }
            else {
                //Assuming binary
                ITERATOR child1 = p.left();
                ITERATOR child2 = child1.right();
                Scalar t1 = child1->length;
                Scalar t2 = child2->length;
                Scalar partial_i;
                
                for (int i=0;i<4;i++) { //State here
                    Scalar p1 = 0, p2 = 0;
                    for(int j=0;j<4;j++)
                        p1 += model.Pij(i,j,t1)*child1->partials[j];
                    for(int j=0;j<4;j++)
                        p2 += model.Pij(i,j,t2)*child2->partials[j];
                    partial_i = p1*p2;
                    p->partials[i]=partial_i;
                }
                
                Scalar max_partial_i=*max_element(p->partials.begin(),p->partials.end());
                //Deal with underflow.
                p->exponent = child1->exponent + child2->exponent;
                while (max_partial_i > 0 && max_partial_i < UNDERFLOW_CUTOFF) {
                    max_partial_i *= UNDERFLOW_MULTIPLIER;
                    for (int i=0;i<4;i++) {
                        p->partials[i]*=UNDERFLOW_MULTIPLIER;
                    }
                    p->exponent-=UNDERFLOW_STEP;
                }
                
            }
            
        }
        
        Scalar Ls = 0.0;
        for(int i=0;i<4;i++)
            Ls += model.pi(i)*tree.root()->partials[i];
        
        siteL[s] = log(Ls) + tree.root()->exponent;
       //if (s<5)
         //  cout<<s<<"\t"<<siteL[s]<<"\t"<<tree.root()->exponent<<endl;
        logL += siteL[s];
        
    }
    timer.stop();

    
    return logL;
}


static void updatePartialsRecurse(phylo<nodeData>::iterator p, phylo<nodeData>& tree, const SubstModel& model) {

  using ITERATOR = phylo<nodeData>::iterator;
  
    //Check that any dirty children have been updated
    for (auto q = p.left(); !q.null(); q = q.right()) {
        if (q->isDirty)
            updatePartialsRecurse(q,tree,model);
    }
    //Update partials on this node, depending on the merger type
    //Note - if this node is dirty then it has children.
    
    ITERATOR child1 = p.left();
    ITERATOR child2 = child1.right();
    Scalar t1 = child1->length;
    Scalar t2 = child2->length;
    Scalar partial_i;
    for (int i=0;i<4;i++) { //State here
        Scalar p1 = 0, p2 = 0;
        for(int j=0;j<4;j++)
            p1 += model.Pij(i,j,t1)*child1->partials[j];
        for(int j=0;j<4;j++)
            p2 += model.Pij(i,j,t2)*child2->partials[j];
        partial_i = p1*p2;
        p->partials[i]=partial_i;
    }
    //Deal with underflow.
    Scalar max_partial_i = *max_element(p->partials.begin(),p->partials.end());
    p->exponent = child1->exponent + child2->exponent;
    while (max_partial_i > 0 && max_partial_i < UNDERFLOW_CUTOFF) {
        max_partial_i *= UNDERFLOW_MULTIPLIER;
        for (int i=0;i<4;i++) {
            p->partials[i]*=UNDERFLOW_MULTIPLIER;
        }
        p->exponent-=UNDERFLOW_STEP;
    }
    
    
    p->isDirty = false;
    
   /*  cout<<"End of update Partial "<<p->length<<" "<<p->id<<endl;   
        debugPrintPartials(cout,decomTree);
    cout<<endl;
 */
}


/**
 Update the partial likelihoods for all the nodes ancestral to taxa for which there is a change.
 This method assume that all notes have the isDirty flag set to false.
 */
static Scalar updatePartials(phylo<nodeData>& tree, vector<phylo<nodeData>::iterator>& taxaPointers, const SubstModel& model, vector<pair<int, unsigned short>>& diffs) {
    //First mark the nodes which need updating
    for(unsigned int i=0;i<diffs.size();i++) {
        auto p = taxaPointers[diffs[i].first];
        //if (p.null())
         //   continue; //Taxa not in tree, so skip.
        base b = diffs[i].second;
        
        resolve_base(b,p->partials);
        
        p->isDirty = false; //This node is now updated
        p=p.par();
        
        //Mark vertices above - we can stop as soon as we meet one which is already marked
        while (p!=tree.header() && !p->isDirty) {
            p->isDirty=true;
            p=p.par();
        }
    }

    updatePartialsRecurse(tree.root(),tree,model);
        
    Scalar Ls = 0.0;
    for(int i=0;i<4;i++)
        Ls += model.pi(i)*(tree.root()->partials[i]);
        
   
    
    
    return log(Ls) + tree.root()->exponent;
    
}




Scalar  computeLikelihoodUsingUpdating(phylo<nodeData>& tree, const SubstModel& model, const vector<sequence>& seqs, vector<Scalar>& siteL, PatternSorter patternSorter, Stopwatch& timer) {
 
    timer.stop();
    
   using ITERATOR = phylo<nodeData>::iterator;
    PatternVector patterns = compressPatterns(seqs,patternSorter);
    auto differences = computePatternDifferences(patterns);
    int ntax = seqs.size();
    vector<phylo<nodeData>::iterator> taxaPointers(ntax);
    for(int i=0; i<ntax;i++)
        taxaPointers[i].set_null();
        
    //Initialise partial arrays and store pointers to the taxa.
    for(ITERATOR p = tree.leftmost_leaf();!p.null();p=p.next_post()) {           
        p->partials = {0,0,0,0};
        p->exponent = 0;
        p->isDirty = false;
        if (p.leaf()) {
            taxaPointers[p->id]=p;
            //cout<<"p->id = "<<p->id<<endl;
        }
    }
    
    timer.start();
    Scalar logL = 0.0;
    
    
    
    for (unsigned int i=0;i<patterns.size();i++) {
        Scalar patternlogL = updatePartials(tree,taxaPointers,model,differences[i]);
        
        logL += patternlogL * patterns[i].second;
       // if (i<5)
       //     cout<<"Pattern "<<i<<" logL = "<<patternlogL<<" freq = "<<patterns[i].second<<endl;
    }
    timer.stop();
    
    return logL;
}
 
 
 
 
 
 
 
 
