/*
 *
 *  LikelihoodAgreement.cpp
 *
 *  Checks the numerical agreement between the log-likelihoods computed using the
 *  pruning decomposition tree and the LvD decomposition tree. Each replicate simulates
 *  a uniform tree and an alignment, then evaluates the likelihood of the true tree
 *  with both decomposition trees on identical patterns.
 *
 */


#include <cmath>
#include <iomanip>
#include "Phylib/phylib.h"
#include "standardLikelihood.h"
#include "DecompositionTree/decompositionTree.h"
#include "DecompositionTree/decompositionLikelihood.h"
#include "DecomMCMC.h"

using namespace Phylib;


/**
 Parsing the program arguments
 */

class SimulationOptions
{
public:
    int min_ntax = 100;
    int max_ntax = 10000;
    int nsites = 1000;
    double root_to_tip = 0.0; //Expected number of mutations from root to a tip
    int num_replicates = 0;
    bool outputHeader = false;
};


static string Usage()
{
    string s;
    s = "LikelihoodAgreement [-h] -r <height> -R <replicates> [-m <min_ntax>] [-n <max_ntax>] [-s <nsites>]\n\n";
    s += "\t -h\t\t Output the header at top (must be the first argument)\n";
    s += "\t -r  <height>\t Expected number of mutations from root to a tip (required)\n";
    s += "\t -R  <number>\t Number of replicates (required)\n";
    s += "\t -m  <min_ntax>\t Minimum number of taxa (default: 100)\n";
    s += "\t -n  <max_ntax>\t Maximum number of taxa (default: 10000)\n";
    s += "\t -s  <nsites>\t Number of sites (default: 1000)\n";
    return s;
}

static bool is_int(const string& s, int& value)
{
    try
    {
        size_t pos;
        int v = std::stoi(s, &pos, 10);
        if (pos != s.size()) return false; // leftover chars
        value = v;
        return true;
    }
    catch (...)
    {
        return false; // invalid or out of range
    }
}

static bool is_double(const string& s, double& value)
{
    try
    {
        size_t pos;
        double v = std::stod(s, &pos);
        if (pos != s.size()) return false; // leftover chars
        value = v;
        return true;
    }
    catch (...)
    {
        return false; // invalid or out of range
    }
}

static bool parseArguments(int argc, char* argv[], SimulationOptions& options, string& errmsg)
{
    int first = 1;
    if (argc > 1 && string(argv[1]) == "-h")
    {
        options.outputHeader = true;
        first = 2;
    }

    bool have_r = false, have_R = false;

    for (int i = first; i < argc; i++)
    {
        string key = string(argv[i]);
        if (i + 1 >= argc)
        {
            errmsg = "ERROR: missing value for " + key + "\n\n" + Usage();
            return false;
        }
        string val = string(argv[++i]);

        bool ok;
        if (key == "-r")
            ok = have_r = is_double(val, options.root_to_tip);
        else if (key == "-R")
            ok = have_R = is_int(val, options.num_replicates);
        else if (key == "-m")
            ok = is_int(val, options.min_ntax);
        else if (key == "-n")
            ok = is_int(val, options.max_ntax);
        else if (key == "-s")
            ok = is_int(val, options.nsites);
        else
        {
            errmsg = "ERROR: unrecognised option " + key + "\n\n" + Usage();
            return false;
        }
        if (!ok)
        {
            errmsg = "ERROR: invalid value for " + key + "\n\n" + Usage();
            return false;
        }
    }

    if (!have_r || options.root_to_tip <= 0.0)
    {
        errmsg = "ERROR: -r <height> is required and must be positive\n\n" + Usage();
        return false;
    }
    if (!have_R || options.num_replicates < 1)
    {
        errmsg = "ERROR: -R <replicates> is required and must be positive\n\n" + Usage();
        return false;
    }
    if (options.min_ntax < 3)
    {
        errmsg = "ERROR: -m <min_ntax> must be at least 3\n\n" + Usage();
        return false;
    }
    if (options.min_ntax > options.max_ntax)
    {
        errmsg = "ERROR: -m <min_ntax> must not exceed -n <max_ntax>\n\n" + Usage();
        return false;
    }
    if (options.nsites < 1)
    {
        errmsg = "ERROR: -s <nsites> must be positive\n\n" + Usage();
        return false;
    }
    return true;
}


/**
 Code for simulating sequences on a tree (as in RunSimulation)
 */
class simNodeData : public basic_newick
{
public:
    simNodeData() : basic_newick()
    {
        state = 0;
    }

    simNodeData(const basic_newick& node) : basic_newick(node)
    {
        state = 0;
    }

    int state;
    vector<vector<double>> transitionMat;
};

static void simSequences(phylo<basic_newick>& T, const SubstModel& model, int nsites, vector<sequence>& alignment)
{
    phylo<simNodeData> simT;
    copy<basic_newick, simNodeData>(T, simT);

    for (auto p = simT.leftmost_leaf(); !p.root(); p = p.next_post())
    {
        p->transitionMat.resize(4);
        for (int i = 0; i < model.num_states(); i++)
        {
            p->transitionMat[i].resize(4);
            for (int j = 0; j < model.num_states(); j++)
                p->transitionMat[i][j] = model.Pij(i, j, p->length);
        }
    }
    vector<double> pi_vec(4);
    for (int i = 0; i < 4; i++)
        pi_vec[i] = model.pi(i);

    int ntax = alignment.size();
    for (int i = 0; i < ntax; i++)
        alignment[i].resize(nsites);

    for (int i = 0; i < nsites; i++)
    {
        for (auto p = simT.root(); !p.null(); p = p.next_pre())
        {
            if (p.root())
                p->state = random_discrete(pi_vec);
            else
                p->state = random_discrete(p->transitionMat[p.par()->state]);
            if (p.leaf())
                alignment[p->id][i] = p->state;
        }
    }
    simT.clear();
}


/**
 Log-likelihood of the tree, using a decomposition tree with the given layout.
 Uses the same likelihood routine as the MCMC in RunSimulation.
 */
static Scalar decompositionLogLikelihood(phylo<basic_newick>& tree, bool useLvD, const SubstModel& model,
                                         const vector<pair<Pattern, int>>& patterns)
{
    phylo<DecomNodeData> decomTree;
    if (useLvD)
        constructDecompTree(tree, decomTree, false);
    else
        constructPruningDecomp(tree, decomTree);

    phylo<DecomNodeDataAllSites> t;
    copy(decomTree, t);
    decomTree.clear();

    vector<double> patternL;
    Stopwatch timer;
    Scalar logL = computeLikelihood(t, model, patterns, patternL, timer);
    t.clear();
    return logL;
}


int main(int argc, char* argv[])
{
    SimulationOptions options;
    string err_msg;

    if (!parseArguments(argc, argv, options, err_msg))
    {
        cerr << err_msg << endl;
        exit(1);
    }

    if (options.outputHeader)
    {
        cout << "Replicate\t"
            << "ntax\t"
            << "root to tip\t"
            << "nsites\t"
            << "nPatterns\t"
            << "logL pruning\t"
            << "logL LvD\t"
            << "abs difference"
            << endl;
    }

    JCModel model;

    for (int rep = 0; rep < options.num_replicates; rep++)
    {
        //Simulate ntax so that log(ntax) is uniform between min_ntax and max_ntax (both reachable).
        int ntax = (int)lround(pow(10, log10(options.min_ntax)
                                       + randu() * (log10(options.max_ntax) - log10(options.min_ntax))));

        phylo<basic_newick> simT = sim_uniform_phylo<basic_newick>(ntax);
        scale(simT, options.root_to_tip / root_to_tip<basic_newick>(simT.root()));
        assignBranchIDs(simT);

        vector<sequence> alignment(ntax);
        simSequences(simT, model, options.nsites, alignment);
        vector<pair<Pattern, int>> patterns = sortPatterns(alignment);

        Scalar logLPruning = decompositionLogLikelihood(simT, false, model, patterns);
        Scalar logLLvD = decompositionLogLikelihood(simT, true, model, patterns);

        //Full precision, since the differences of interest are close to machine epsilon.
        cout << rep + 1 << "\t"
            << ntax << "\t"
            << options.root_to_tip << "\t"
            << options.nsites << "\t"
            << patterns.size() << "\t"
            << std::setprecision(17)
            << logLPruning << "\t"
            << logLLvD << "\t"
            << std::fabs(logLPruning - logLLvD)
            << std::setprecision(6)
            << endl;

        simT.clear();
    }
    return (0);
}
