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
#include "../Reference/PartialLikelihoodTensorReference.h"

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

 partialVec is 4 x size: column s holds the four state partial likelihoods for site s.
   Only allocated when isVectors == true.

 partialMat is 4 x (4*size): the s-th 4x4 matrix occupies middleCols(4*s, 4).
   Only allocated when isVectors == false.

 exponents is a per-site underflow correction vector (one entry per site).
 */
class DecomNodeDataAllSites : public basic_newick {
public:
    bool external;
    bool isClade;
    bool isDirty;
    int mergeType;
    int nSites;
    PartialLikelihoodTensorReference partialLikes;

//    Eigen::Matrix<Scalar, 4, Eigen::Dynamic> partialVec;  // 4 x size        (clades only)
//    Eigen::Matrix<Scalar, 4, Eigen::Dynamic> partialMat;  // 4 x (4*size)    (segments only)
    Eigen::VectorXi exponents;                             // one per site


    // Default constructor: does not allocate partial likelihood storage.
    DecomNodeDataAllSites() : basic_newick(), external(false), isClade(false),
                              isDirty(false), nSites(0), mergeType(0) {}

    // Conversion constructor from DecomNodeData: copies topology fields, leaves partials unallocated.
    // Required by Phylib::copy<DecomNodeData, DecomNodeDataAllSites>.
     DecomNodeDataAllSites(const DecomNodeData& d)
         : basic_newick(static_cast<const basic_newick&>(d)),
           external(d.external), isClade(d.isClade),
           isDirty(d.isDirty), nSites(0), mergeType(d.mergeType) {}

    // Copy constructor: allocates new partial Likelihoods
    DecomNodeDataAllSites(const DecomNodeDataAllSites& d) : external(d.external), isClade(d.isClade),
          isDirty(d.isDirty), nSites(d.nSites), mergeType(d.mergeType), partialLikes(d.partialLikes) {

        exponents.resize(nSites);
        exponents.setZero();
    }
};


void constructPruningDecomp( phylo<basic_newick>& input, phylo<DecomNodeData>& dtree);
void constructDecompTree( phylo<basic_newick>& input, phylo<DecomNodeData>& decomTree, bool useGreedy);









#endif

