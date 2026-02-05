/**
 * @file standardLikelihood.h
 * @brief Header file for routines for calculating likelihoods using the standard algorithms
 * @author David Bryant & Celine Scornavacca
 * @version 0.1
 */

#ifndef STANDARDLIKE_H
#define STANDARDLIKE_H

#include "phylib.h"


using namespace Phylib;

/**
 Parameters for controlling underflow.
 A probability is stored as prob  = value * 10^(exponent)
 If value < UNDERFLOW_CUTOFF  we multiply value by e^UNDERFLOW_STEP and decrease
 exponent by UNDERFLOW_STEP
 */
static const Scalar UNDERFLOW_CUTOFF = 1e-20;
static const Scalar UNDERFLOW_STEP = 20;
static const Scalar UNDERFLOW_MULTIPLIER = exp(UNDERFLOW_STEP);


/**
 Declare the data stored at each node. In this case we store just the partial likleihoods
 for a single site (i.e. no compression). Next iteration would be to implement an
 interface and class for handling compressed patterns.
 **/

class nodeData :public basic_newick  {
public:
    nodeData() : basic_newick() {
        partials = {0,0,0,0};
        exponent = 0;
        isDirty = false;
    }
    nodeData(const basic_newick& node) : basic_newick(node) {
        partials = {0,0,0,0};
        exponent = 0;
        isDirty = false;
        
    }
    
    vector< Scalar > partials;
    int exponent; //Actual probability is partials x 10^exponent
    bool isDirty;
};


/*
Generic Substitution Model
*/

class SubstModel {
public:
	virtual Scalar pi(base i) const = 0;
	virtual Scalar Pij(base i, base j, Scalar t) const = 0;
	virtual int num_states() const = 0;
};


/*
Implementation of Jukes Cantor Model
*/
class JCModel: public SubstModel {
public:
	Scalar pi(base i) const override {return 0.25;}
	Scalar Pij(base i, base j, Scalar t) const override {
		if (i==j)
			return 0.25 + 0.75 * exp(-(4.0/3)*t);
		else
			return 0.25 - 0.25 * exp((-4.0/3)*t);
	}
	int num_states() const override {return 4;};
};

/*
 ChatGPT implementation of HKY model
 */

class HKY85ExplicitModel : public SubstModel {
public:
    // Pass any positive freqs; we normalize internally.
    HKY85ExplicitModel(double kappa, double piA, double piC, double piG, double piT)
      : kappa_(kappa)
    {
        double s = piA + piC + piG + piT;
        pis_ = { piA/s, piC/s, piG/s, piT/s };
        precompute_();
    }

    int    num_states() const override { return 4; }
    double pi(base b) const override { return pis_[static_cast<int>(b)]; }

    // t is a branch length in expected substitutions/site (rate = 1)
    double Pij(base bi, base bj, double t) const override {
        const int i = static_cast<int>(bi);
        const int j = static_cast<int>(bj);

        const bool iR = is_purine(i);
        const bool jR = is_purine(j);

        const double pi_i = pis_[i];
        const double pi_j = pis_[j];

        // Precomputed group constants
        const double piG   = iR ? piR_   : piY_;
        const double piNot = iR ? piY_   : piR_;
        const double gamma = iR ? gammaR_ : gammaY_;

        const double e1 = std::exp(-beta_ * t);           // exp(-t / mu0)
        const double e2 = std::exp(-gamma * beta_ * t);   // exp(-gamma_G * t / mu0)

        double P = 0.0;
        if (iR != jR) {
            // transversion (cross-group)
            P = pi_j * (1.0 - e1);
        } else if (i == j) {
            // stay in same state (diagonal)
            P = (pi_i * (piG + piNot * e1) + (piG - pi_i) * e2) / piG;
        } else {
            // transition (within-group, to j != i)
            P = (pi_j * (piG + piNot * e1) - pi_j * e2) / piG;
        }

        // numeric guard
        if (P < 0.0)      P = (P > -1e-12) ? 0.0 : P;
        else if (P > 1.0) P = (P < 1.0+1e-12) ? 1.0 : P;
        return P;
    }

private:
    static bool is_purine(int idx) { return idx == 0 /*A*/ || idx == 2 /*G*/; }

    void precompute_() {
        piR_ = pis_[0] + pis_[2];     // A + G
        piY_ = pis_[1] + pis_[3];     // C + T
        // average rate for the unscaled HKY off-diagonals q_ij ∝ (κ for transitions, 1 otherwise)*π_j
        mu0_ = 2.0 * ( piR_*piY_ + kappa_*(pis_[0]*pis_[2] + pis_[1]*pis_[3]) );
        beta_   = 1.0 / mu0_;                       // scale so t = expected subs/site
        gammaR_ = 1.0 + piR_ * (kappa_ - 1.0);
        gammaY_ = 1.0 + piY_ * (kappa_ - 1.0);
    }

    double kappa_;
    std::array<double,4> pis_;   // [A,C,G,T] sums to 1

    // cached constants
    double piR_{}, piY_{}, mu0_{}, beta_{}, gammaR_{}, gammaY_{};
};

/**
 Compute the tree likelihoods.
 
 It assumes that the ID of the leaf nodes in the tree corresponds to the
 corresponding row in seqs. This version has no compression.
 **/

Scalar  computeLikelihood(phylo<nodeData>& tree, const SubstModel& model, const vector<sequence>& seqs, vector<Scalar>& siteL, Stopwatch& timer);
Scalar  computeLikelihoodUsingUpdating(phylo<nodeData>& tree, const SubstModel& model, const vector<sequence>& seqs, vector<Scalar>& siteL, PatternSorter patternSorter, Stopwatch& timer);



#endif
