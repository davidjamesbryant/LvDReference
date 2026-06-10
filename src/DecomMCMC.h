/*
 *
 *  DecomMCMC.h
 *
 *  Shared MCMC infrastructure for decomposition-tree branch-length sampling.
 *  Used by both RunSimulation and BranchLengthMCMC.
 *
 */

#pragma once

#include <vector>
#include "Phylib/phylib.h"
#include "standardLikelihood.h"
#include "DecompositionTree/decompositionTree.h"
#include "DecompositionTree/decompositionLikelihood.h"

using namespace Phylib;

struct MCMCResults {
    vector<double> logPosterior;  // log-posterior at each accepted/proposed state
    vector<double> iterTime;      // wall-clock seconds consumed by each iteration
    // DEBUG fields — comment out before production
    vector<int>    branch;        // index of branch proposed at each iteration
    vector<double> old_length;    // branch length before proposal
    vector<double> new_length;    // branch length after proposal
    vector<bool>   accepted;      // whether the proposal was accepted
    double         initialLogPost = 0.0;
    int            numBranches   = 0;
    int            numAccepted   = 0;
    int            numNegative   = 0;   // proposals rejected because newLen <= 0
};

// Assigns stable branch IDs: leaves keep ids 0..ntax-1; internal non-root
// nodes get ids -(ntax), -(ntax+1), ... so abs(p->id) indexes branches.
void assignBranchIDs(phylo<basic_newick>& tree);

// Core MH random-walk over all branch lengths on an already-built decomposition
// tree.  Caller is responsible for seeding the RNG before calling.
MCMCResults runDecomMCMC(
    phylo<DecomNodeDataAllSites>&    t,
    const SubstModel&                model,
    const vector<pair<Pattern,int>>& patterns,
    double prior_branch_rate,
    double proposal_width,
    int    num_iterations);
