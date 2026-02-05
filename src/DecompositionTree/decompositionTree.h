/**
 * @file decompositionTree.h
 * @brief Header file for routines for constructing and updating decomposition trees
 * @author David Bryant & Celine Scornavacca
 * @version 0.1
 */

#ifndef DECOMPTREE_H
#define DECOMPTREE_H

#include <queue>
#include <Eigen>
#include "phylib.h"
#include "../StandardLikelihood/standardLikelihood.h"

using namespace Phylib;

/*
We need different fields during the construction of the decomposition tree as when it is actually used. This header file just concentrates on the later. A decomposition tree is a special type of tree where the leaves correspond to edges in the original phylogeny. This is stored as an iterator for the tree.
 */


class DecomNodeData : public basic_newick{
public:
    bool external;  //Does this node correspond to a pendant edge in the phylogeny
    bool isClade;
    bool isDirty;
    
    //Store the partial likelihood information. We store both a vector and a matrix. This creates redundancy (we only need one),
    //but I figure that this allows compilation time improvements for minimal extra memory usage.
    //An alternative is to have a different class for clades and segments, both of which inherit from DecomNodeData
    
    Eigen::Matrix<Scalar, 4, 1>  partialVec;
    Eigen::Matrix<Scalar, 4, 4>  partialMat;
    
    int mergeType;
    int exponent;
    
    DecomNodeData() : Phylib::basic_newick() {
        external = false;
        exponent = 0;
        mergeType = 0;
    }
};

void constructPruningDecomp( phylo<basic_newick>& input, phylo<DecomNodeData>& dtree);
void constructDecompTree( phylo<basic_newick>& input, phylo<DecomNodeData>& decomTree, bool useGreedy);








Scalar  computeLikelihood(phylo<DecomNodeData>& decomTree, const SubstModel& model, const vector<sequence>& seqs, vector<Scalar>& siteL, Stopwatch& timer);
Scalar computeLikelihoodUsingUpdating(phylo<DecomNodeData>& decomTree, const SubstModel& model, const vector<sequence>& seqs, PatternSorter patternSorter, Stopwatch& timer);


Scalar computeLikelihood(phylo<DecomNodeData>& decomTree, const SubstModel& model, const vector<pair<Pattern, int>>& patterns, vector<double>& patternL, Stopwatch& timer);

Scalar computeLikelihood(phylo<DecomNodeData>& decomTree, const SubstModel& model, const vector<pair<Pattern, int>>& patterns, const PatternDiffs& differences, vector<double>& patternL, Stopwatch& timer);




#endif

