
/**
 * @file decompositionLikelihood.h
 * @brief Header file for routines for constructing and updating likelihood calculations using decomposition trees
 * @author David Bryant & Celine Scornavacca
 * @version 0.1
 */

#ifndef DECOMPTREELIKE_H
#define DECOMPTREELIKE_H

#include <queue>
#include <Eigen>
#include "phylib.h"
#include "../StandardLikelihood/standardLikelihood.h"
#include "decompositionTree.h"


using namespace Phylib;

/*
We need different fields during the construction of the decomposition tree as when it is actually used. This header file just concentrates on the later. A decomposition tree is a special type of tree where the leaves correspond to edges in the original phylogeny. This is stored as an iterator for the tree.
 */

Scalar  computeLikelihood(phylo<DecomNodeData>& decomTree, const SubstModel& model, const vector<sequence>& seqs, vector<Scalar>& siteL, Stopwatch& timer);

Scalar computeLikelihoodUsingUpdating(phylo<DecomNodeData>& decomTree, const SubstModel& model, const vector<sequence>& seqs, PatternSorter patternSorter, Stopwatch& timer);

Scalar computeLikelihood(phylo<DecomNodeData>& decomTree, const SubstModel& model, const vector<pair<Pattern, int>>& patterns, vector<double>& patternL, Stopwatch& timer);

Scalar computeLikelihood(phylo<DecomNodeData>& decomTree, const SubstModel& model, const vector<pair<Pattern, int>>& patterns, const PatternDiffs& differences, vector<double>& patternL, Stopwatch& timer);




#endif

