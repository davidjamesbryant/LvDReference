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

/**
 Like DecomNodeData but stores partial likelihoods for all sites/patterns at once.

 partialVec is 4 x nSites: column s holds the four state partial likelihoods for site s.
   Only allocated when isClade == true.

 partialMat is 4 x (4*nSites): the s-th 4x4 matrix occupies middleCols(4*s, 4).
   Only allocated when isClade == false.

 exponents is a per-site underflow correction vector (one entry per site).
 */
class DecomNodeDataAllSites : public basic_newick {
public:
    bool external;
    bool isClade;
    bool isDirty;

    int nSites;
    Eigen::Matrix<Scalar, 4, Eigen::Dynamic> partialVec;  // 4 x nSites        (clades only)
    Eigen::Matrix<Scalar, 4, Eigen::Dynamic> partialMat;  // 4 x (4*nSites)    (segments only)
    Eigen::VectorXi exponents;                             // one per site

    int mergeType;

    // Default constructor: does not allocate partial likelihood storage.
    DecomNodeDataAllSites() : basic_newick(), external(false), isClade(false),
                              isDirty(false), nSites(0), mergeType(0) {}

    // Allocating constructor: allocates partialVec (clade nodes) or partialMat (segment nodes).
    DecomNodeDataAllSites(bool isClade_, int nSites_)
        : basic_newick(), external(false), isClade(isClade_),
          isDirty(false), nSites(nSites_), mergeType(0) {
        if (isClade) {
            partialVec.resize(4, nSites);
            partialVec.setZero();
        } else {
            partialMat.resize(4, 4 * nSites);
            partialMat.setZero();
        }
        exponents.resize(nSites);
        exponents.setZero();
    }
};


void constructPruningDecomp( phylo<basic_newick>& input, phylo<DecomNodeData>& dtree);
void constructDecompTree( phylo<basic_newick>& input, phylo<DecomNodeData>& decomTree, bool useGreedy);









#endif

