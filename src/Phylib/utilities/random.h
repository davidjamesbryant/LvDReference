/**
 *  Collection of random generation utilities (modernized).
 *
 *  Header-only; C++17.
 *  - Thread-local PRNG (mt19937_64) seeded from std::random_device
 *  - Distributions from <random>
 *
 *  Backward-compatible free functions retained:
 *    seed_random(), seed_random(int), random_num(x), randu(), randu(a,b),
 *    random_exp(mean), random_gaussian(), random_gamma(alpha, beta),
 *    random_discrete(p)
 *
 *  New class:
 *    Phylib::RandomEngine   // explicit engine instance with the same API
 *
 *  Author: David Bryant (modernization by ChatGPT)
 */

#ifndef RANDOM_H_INCLUDE
#define RANDOM_H_INCLUDE

#include <random>
#include <vector>
#include <numeric>
#include <algorithm>
#include <cmath>
#include <stdexcept>
#include "../global/stdIncludes.h"   // keep if your project expects this

namespace Phylib {

// ---------- Engine plumbing (thread-local) -----------------------------------

struct RandomEngine {
    using urbg_t = std::mt19937_64;

    // Thread-local engine seeded once by default
    static urbg_t& engine() {
        thread_local urbg_t eng{ seed_from_rd() };
        return eng;
    }

    // Reseed with non-deterministic seed
    static void reseed() {
        engine().seed(seed_from_rd());
    }
    // Reseed with fixed seed (deterministic)
    static void reseed(uint64_t seed) {
        engine().seed(seed);
    }

    // Uniform integer in [0, x-1]; returns 0 if x==0
    static unsigned int random_num(unsigned int x) {
        if (x == 0) return 0u;
        std::uniform_int_distribution<unsigned int> D(0u, x - 1);
        return D(engine());
    }

    // Uniform double in [0,1)
    static double randu() {
        // uniform_real_distribution is half-open by default
        static thread_local std::uniform_real_distribution<double> U(0.0, 1.0);
        return U(engine());
    }

    // Uniform double in [lower, upper)
    static double randu(double lower, double upper) {
        std::uniform_real_distribution<double> U(lower, upper);
        return U(engine());
    }

    // Exponential with given mean (rate = 1/mean)
    static double random_exp(double mean) {
        if (!(mean > 0.0)) throw std::invalid_argument("mean must be > 0");
        std::exponential_distribution<double> E(1.0 / mean);
        return E(engine());
    }

    // Standard normal N(0,1)
    static double random_gaussian() {
        static thread_local std::normal_distribution<double> N(0.0, 1.0);
        return N(engine());
    }

    // Gamma(shape=alpha, scale=beta). alpha>0, beta>0
    static double random_gamma(double alpha, double beta = 1.0) {
        if (!(alpha > 0.0 && beta > 0.0))
            throw std::invalid_argument("alpha, beta must be > 0");
        std::gamma_distribution<double> G(alpha, beta);
        return G(engine());
    }

    // Sample index i in {0..p.size()-1} with probabilities proportional to p[i].
    // (If you want 1..m indexing, add +1 at the callsite.)
    static std::size_t random_discrete(const std::vector<double>& p) {
        if (p.empty()) throw std::invalid_argument("p must be non-empty");
        // std::discrete_distribution tolerates unnormalized, even zero weights
        std::discrete_distribution<std::size_t> D(p.begin(), p.end());
        return D(engine());
    }

    
    
    
private:
    static urbg_t::result_type seed_from_rd() {
        std::random_device rd;
        // Combine multiple rd() calls to better fill 64 bits if rd is 32-bit
        std::seed_seq seq{ rd(), rd(), rd(), rd(), rd(), rd(), rd(), rd() };
        std::vector<urbg_t::result_type> seeds(1);
        seq.generate(seeds.begin(), seeds.end());
        return seeds[0];
    }
};

// ---------------- Backward-compatible free functions --------------------------

inline void seed_random() { RandomEngine::reseed(); }
inline void seed_random(int r_seed) { RandomEngine::reseed(static_cast<uint64_t>(r_seed)); }

inline unsigned int random_num(unsigned int x) {
    return RandomEngine::random_num(x);
}

// Uniform(0,1)
inline double randu() { return RandomEngine::randu(); }

// Uniform(lower, upper)
inline double randu(double lower, double upper) { return RandomEngine::randu(lower, upper); }

// Exponential(mean)
inline double random_exp(double mean) { return RandomEngine::random_exp(mean); }

// Standard normal
inline double random_gaussian() { return RandomEngine::random_gaussian(); }

// Gamma(alpha, beta)
inline double random_gamma(double alpha, double beta = 1.0) {
    return RandomEngine::random_gamma(alpha, beta);
}

// Discrete draw with weights p (returns index 0..p.size()-1)
inline std::size_t random_discrete(const std::vector<double>& p) {
    return RandomEngine::random_discrete(p);
}

//Shuffle elements in a range
template <class It>
inline void shuffle(It first, It last) {
    std::shuffle(first, last, RandomEngine::engine());
}

//Shuffle elements in a container
template <class Container>
inline void shuffle(Container& c) {
    using std::begin; using std::end;
    std::shuffle(begin(c), end(c), RandomEngine::engine());
}

//Returns random subset of size k from {0,1,...,n-1}}
//Not very efficient.
inline std::vector<bool> random_subset(int n, int k) {
    std::vector<bool> subset(n);
    std::fill(subset.begin(),subset.end(),false);
    int m=0; //Number of elements already in the set
    for(int i=0;i<n && m<k;i++)
        if (randu()*n<(k-m)) {
            subset[i]=true;
            m++;
        }
    return subset;
}


} // namespace Phylib

#endif // RANDOM_H_INCLUDE
