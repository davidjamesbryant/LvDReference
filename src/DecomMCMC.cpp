/*
 *
 *  DecomMCMC.cpp
 *
 */

#include <chrono>
#include "DecomMCMC.h"

using namespace Phylib;
using namespace std;

void assignBranchIDs(phylo<basic_newick>& tree)
{
    int ntax = 0;
    for (auto p = tree.leftmost_leaf(); !p.null(); p = p.next_post())
        if (p.leaf()) ntax++;
    int nextId = ntax;
    for (auto p = tree.leftmost_leaf(); !p.null(); p = p.next_post())
        if (!p.leaf() && !p.root())
            p->id = -nextId++;
}

MCMCResults runDecomMCMC(
        phylo<DecomNodeDataAllSites>&    t,
        const SubstModel&                model,
        const vector<pair<Pattern,int>>& patterns,
        double prior_branch_rate,
        double proposal_width,
        int    num_iterations)
{
    Stopwatch sw;
    vector<double> patternL;
    Scalar logLik = computeLikelihood(t, model, patterns, patternL, sw);

    Scalar logPrior = 0.0;
    for (auto p = t.leftmost_leaf(); !p.null(); p = p.next_post())
        if (p.leaf())
            logPrior += log(prior_branch_rate) - prior_branch_rate * p->length;
    Scalar logPost = logLik + logPrior;

    unsigned int numBranches = 0;
    for (auto p = t.leftmost_leaf(); !p.null(); p = p.next_post())
        if (p.leaf()) numBranches++;
    vector<phylo<DecomNodeDataAllSites>::iterator> branches(numBranches);
    for (auto p = t.leftmost_leaf(); !p.null(); p = p.next_post())
        if (p.leaf())
            branches[abs(p->id)] = p;

    MCMCResults results;
    results.logPosterior.reserve(num_iterations);
    results.iterTime.reserve(num_iterations);
    results.branch.reserve(num_iterations);
    results.old_length.reserve(num_iterations);
    results.new_length.reserve(num_iterations);
    results.accepted.reserve(num_iterations);
    results.initialLogPost = logPost;
    results.numBranches    = numBranches;

    for (int iter = 0; iter < num_iterations; ++iter) {
        unsigned int branchIdx = random_num(numBranches);
        auto   p         = branches[branchIdx];
        double oldLen    = p->length;
        double newLen    = oldLen + randu(-proposal_width, proposal_width);

        bool accept = false;
        double likeTime = 0.0;
        if (newLen > 0.0) {
            auto t0 = chrono::steady_clock::now();
            Scalar newLogLik   = updateBranchLength(t, model, patterns, patternL, p, newLen);
            likeTime += chrono::duration<double>(chrono::steady_clock::now() - t0).count();

            Scalar newLogPrior = logPrior + prior_branch_rate * (oldLen - newLen);
            Scalar newLogPost  = newLogLik + newLogPrior;

            accept = (log(randu()) < newLogPost - logPost);

            if (accept) {
                logLik   = newLogLik;
                logPrior = newLogPrior;
                logPost  = newLogPost;
            } else {
                auto t1 = chrono::steady_clock::now();
                updateBranchLength(t, model, patterns, patternL, p, oldLen);  // restore
                likeTime += chrono::duration<double>(chrono::steady_clock::now() - t1).count();
            }
        } else {
            results.numNegative++;
        }

        if (accept) results.numAccepted++;

        results.logPosterior.push_back(logPost);
        results.iterTime.push_back(likeTime);
        results.branch.push_back(branchIdx);
        results.old_length.push_back(oldLen);
        results.new_length.push_back(newLen);
        results.accepted.push_back(accept);
    }

    return results;
}
