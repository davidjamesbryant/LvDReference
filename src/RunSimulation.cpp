/*
 *
 *  Created by David Bryant on 28/03/25.
 *
 */


#include <unistd.h>
#include <chrono>
#include <list>
#include <charconv>
#include "Phylib/phylib.h"
#include "standardLikelihood.h"
#include "DecompositionTree/decompositionTree.h"
#include "DecompositionTree/decompositionLikelihood.h"

using namespace Phylib;


/**
 Parsing the program arguments
 */

class SimulationOptions
{
public:
    //Tree distribution parameters
    int ntax = 1000; //Number of taxa
    bool test_uniform = false; //Uniform on trees
    bool test_yule = false; //Yule = coalescent
    bool test_caterpillar = false; //Same as uniform
    bool test_beta_crit = false; //Beta critical
    double root_to_tip = 0.01; //Expected number of mutations from root to a tip

    //Sequence distribution parameters
    int nsites;

    //What to output
    bool runMCMC_PruningTree = true;
    bool runMCMC_DecomTree = true;
    bool outputHeader = true;
    bool outputTrees = false;
    int num_replicates;
    int num_iterations = 1000;
    double prior_branch_rate = 20.0;
    double proposal_width = 0.01;
};

std::ostream& operator<<(std::ostream& os, const SimulationOptions& opt)
{
    // Save current formatting flags so we can restore them at the end
    std::ios::fmtflags f = os.flags();
    os << std::boolalpha; // print bools as true/false

    os << "SimulationOptions {\n"
        << "  ntaxa = " << opt.ntax << "\n"
        << "  test_uniform = " << opt.test_uniform << "\n"
        << "  test_yule = " << opt.test_yule << "\n"
        << "  test_caterpillar = " << opt.test_caterpillar << "\n"
        << "  test_beta_crit = " << opt.test_beta_crit << "\n"
        << "  root_to_tip = " << opt.root_to_tip << "\n"
        << "  nsites = " << opt.nsites << "\n"
        << "  runMCMC_Pruning = " << opt.runMCMC_PruningTree << "\n"
        << "  runMCMC_LvD = " << opt.runMCMC_DecomTree << "\n"
        << "  outputHeader = " << opt.outputHeader << "\n"
        << "  num_replicates = " << opt.num_replicates << "\n"
        << "  num_iterations = " << opt.num_iterations << "\n"
        << "  prior_branch_rate = " << opt.prior_branch_rate << "\n"
        << "  proposal_width = " << opt.proposal_width << "\n"
        << "  ouputTrees = " << opt.outputTrees << "\n"
        << "}";

    os.flags(f); // restore flags
    return os;
}


string Usage()
{
    string s;
    s =
        "RunSimulation -BCUYclhp -n <ntax> -s <nsites> -r <height> -R <replicates> [-i <iterations>] [-P <rate>] [-w <proposal width>]\n\n";
    s += "\t -BCUYclhpt  (flags)\n";
    s += "\t\tB:\t Simulate from beta critical distribution\n";
    s += "\t\tC:\t Simulate from random caterpillar\n";
    s += "\t\tU:\t Simulate from uniform distribution\n";
    s += "\t\tY:\t Simulate from Yule distribution\n";
    s += "\t<At least one tree distribution must be selected>\n";
    s += "\t\tl:\t Output the timing for the MCMC with the LvD algorithm\n";
    s += "\t\th:\t Output the header at top\n";
    s += "\t\tp:\t Output the timing for the MCMC with the pruning algorithm\n";
    s += "\t\tt:\t Output trees to std err\n";
    s += "\t -n  <ntax>\t Number of taxa (required)\n";
    s += "\t -s  <nsites>\t Number of sites (required)\n";
    s += "\t -r  <height>\t Expected number of mutations from root to a tip (required)\n";
    s += "\t -R  <number>\t Number of replicates (required)\n";
    s += "\t -i  <number>\t Number of MCMC iterations (default: 1000)\n";
    s += "\t -P  <rate>\t Prior branch rate for exponential prior (default: 20)\n";
    s += "\t -w  <width>\t Half-width of uniform branch-length proposal (default: 0.01)\n";
    return s;
}

bool is_int(const std::string& s, int& value)
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


bool is_double(const std::string& s, double& value)
{
    try
    {
        value = std::atof(s.c_str());
        return true;
    }
    catch (...)
    {
        return false; // invalid or out of range
    }
}


bool parseArguments(int argc, char* argv[], SimulationOptions& options, string& errmsg)
{
    if (argc < 2)
    {
        errmsg = "ERROR: incorrect number of arguments\n\n" + Usage();
        return false;
    }

    // First argument must be the flags string
    string s = string(argv[1]);
    if (s[0] != '-')
    {
        errmsg = "ERROR: First argument must be the option string\n\n" + Usage();
        return false;
    }
    options.test_uniform = s.find('U') != string::npos;
    options.test_caterpillar = s.find('C') != string::npos;
    options.test_beta_crit = s.find('B') != string::npos;
    options.test_yule = s.find('Y') != string::npos;
    options.runMCMC_PruningTree = s.find('p') != string::npos;
    options.outputHeader = s.find('h') != string::npos;
    options.runMCMC_DecomTree = s.find('l') != string::npos;
    options.outputTrees = s.find('t') != string::npos;

    // Track which required args have been seen
    bool have_n = false, have_s = false, have_r = false, have_R = false;

    // Parse remaining key-value pairs
    for (int i = 2; i < argc; i++)
    {
        string key = string(argv[i]);
        if (i + 1 >= argc)
        {
            errmsg = "ERROR: missing value for " + key + "\n\n" + Usage();
            return false;
        }
        string val = string(argv[++i]);

        if (key == "-n")
        {
            if (!is_int(val, options.ntax))
            {
                errmsg = "ERROR: invalid value for -n\n\n" + Usage();
                return false;
            }
            have_n = true;
        }
        else if (key == "-s")
        {
            if (!is_int(val, options.nsites))
            {
                errmsg = "ERROR: invalid value for -s\n\n" + Usage();
                return false;
            }
            have_s = true;
        }
        else if (key == "-r")
        {
            if (!is_double(val, options.root_to_tip))
            {
                errmsg = "ERROR: invalid value for -r\n\n" + Usage();
                return false;
            }
            have_r = true;
        }
        else if (key == "-R")
        {
            if (!is_int(val, options.num_replicates))
            {
                errmsg = "ERROR: invalid value for -R\n\n" + Usage();
                return false;
            }
            have_R = true;
        }
        else if (key == "-i")
        {
            if (!is_int(val, options.num_iterations))
            {
                errmsg = "ERROR: invalid value for -i\n\n" + Usage();
                return false;
            }
        }
        else if (key == "-P")
        {
            if (!is_double(val, options.prior_branch_rate))
            {
                errmsg = "ERROR: invalid value for -P\n\n" + Usage();
                return false;
            }
        }
        else if (key == "-w")
        {
            if (!is_double(val, options.proposal_width) || options.proposal_width <= 0)
            {
                errmsg = "ERROR: invalid value for -w\n\n" + Usage();
                return false;
            }
        }
        else
        {
            errmsg = "ERROR: unrecognised option " + key + "\n\n" + Usage();
            return false;
        }
    }

    if (!have_n)
    {
        errmsg = "ERROR: -n <ntax> is required\n\n" + Usage();
        return false;
    }
    if (!have_s)
    {
        errmsg = "ERROR: -s <nsites> is required\n\n" + Usage();
        return false;
    }
    if (!have_r)
    {
        errmsg = "ERROR: -r <height> is required\n\n" + Usage();
        return false;
    }
    if (!have_R)
    {
        errmsg = "ERROR: -R <replicates> is required\n\n" + Usage();
        return false;
    }

    return true;
}


/**
 Code for simulating sequences on a tree
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

void simSequences(phylo<basic_newick>& T, const SubstModel& model, int nsites, vector<sequence>& alignment)
{
    //Set up the tree
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

    //Simulate characters
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
            {
                p->state = random_discrete(p->transitionMat[p.par()->state]);
            }
            if (p.leaf())
                alignment[p->id][i] = p->state;
        }
    }
    simT.clear();
}


// Assign stable branch IDs: leaves keep ids 0..ntax-1; internal non-root nodes
// get ids -(ntax), -(ntax+1), ... so that abs(p->id) indexes branches 0..numBranches-1.
static void assignBranchIDs(phylo<basic_newick>& tree)
{
    int ntax = 0;
    for (auto p = tree.leftmost_leaf(); !p.null(); p = p.next_post())
        if (p.leaf()) ntax++;
    int nextId = ntax;
    for (auto p = tree.leftmost_leaf(); !p.null(); p = p.next_post())
        if (!p.leaf() && !p.root())
            p->id = -nextId++;
}


struct MCMCResults {
    vector<double> logPosterior;  // log-posterior at each accepted/proposed state
    vector<double> iterTime;      // wall-clock seconds consumed by each iteration
    // DEBUG fields — comment out before production
    vector<int>    branch;        // index of branch proposed at each iteration
    vector<double> old_length;    // branch length before proposal
    vector<double> new_length;    // branch length after proposal
    vector<bool>   accepted;      // whether the proposal was accepted
    double         initialLogPost  = 0.0;
    int            numBranches    = 0;
    int            numAccepted    = 0;
    int            numNegative    = 0;   // proposals rejected because newLen <= 0
};


// Shared implementation for both decomposition-tree MCMC variants.
// useLvD=false → pruning layout; useLvD=true → LvD (greedy) layout.
static MCMCResults runDecomMCMC(
        phylo<DecomNodeDataAllSites>& t,
        const SubstModel&                model,
        const vector<pair<Pattern,int>>& patterns,
        const SimulationOptions&         options)
{

    Stopwatch sw;
    vector<double> patternL;
    Scalar logLik = computeLikelihood(t, model, patterns, patternL, sw);

    Scalar logPrior = 0.0;
    for (auto p = t.leftmost_leaf(); !p.null(); p = p.next_post())
        if (p.leaf())
            logPrior += log(options.prior_branch_rate) - options.prior_branch_rate * p->length;
    Scalar logPost = logLik + logPrior;

    unsigned int numBranches = 0;
    for (auto p = t.leftmost_leaf(); !p.null(); p = p.next_post())
        if (p.leaf()) numBranches++;
    vector<phylo<DecomNodeDataAllSites>::iterator> branches(numBranches);
    for (auto p = t.leftmost_leaf(); !p.null(); p = p.next_post())
        if (p.leaf())
            branches[abs(p->id)] = p;

    MCMCResults results;
    results.logPosterior.reserve(options.num_iterations);
    results.iterTime.reserve(options.num_iterations);
    results.branch.reserve(options.num_iterations);
    results.old_length.reserve(options.num_iterations);
    results.new_length.reserve(options.num_iterations);
    results.accepted.reserve(options.num_iterations);
    results.initialLogPost = logPost;
    results.numBranches    = numBranches;

    for (int iter = 0; iter < options.num_iterations; ++iter) {
        unsigned int    branchIdx = random_num(numBranches);
        auto   p         = branches[branchIdx];
        double oldLen    = p->length;
        double newLen    = oldLen + randu(-options.proposal_width, options.proposal_width);

        bool accept = false;
        double likeTime = 0.0;
        if (newLen > 0.0) {
            auto t0 = chrono::steady_clock::now();
            Scalar newLogLik   = updateBranchLength(t, model, patterns, patternL, p, newLen);
            likeTime += chrono::duration<double>(chrono::steady_clock::now() - t0).count();

            Scalar newLogPrior = logPrior + options.prior_branch_rate * (oldLen - newLen);
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




int main(int argc, char* argv[])
{
    SimulationOptions options;
    string err_msg;

    if (!parseArguments(argc, argv, options, err_msg))
    {
        cerr << err_msg << endl;
        exit(1);
    }

    std::string test_names[] = {"uniform", "yule", "caterpillar", "beta_crit"};
    bool run_test[] = {options.test_uniform, options.test_yule, options.test_caterpillar, options.test_beta_crit};

    ofstream treefile;

    seed_random(42);

    //Output header
    if (options.outputHeader)
    {
        cout << "Distribution\t"
            << "ntax\t"
            << "root to tip\t"
            << "total tree length\t"
            << "nsites\t"
            << "nPatterns\t"
            << "height pruning tree\t"
            << "av height pruning tree\t"
            << "height LvD tree\t"
            << "av height LvD tree\t";
        if (options.runMCMC_PruningTree)
            cout << "pruning time\t";
        if (options.runMCMC_DecomTree)
            cout << "decom time\t";
        if (options.runMCMC_PruningTree && options.runMCMC_DecomTree)
            cout << "ratio\t";
        cout << endl;
    }


    for (int testnum = 0; testnum < 4; testnum++)
    {
        if (!run_test[testnum])
            continue;

        int ntax = options.ntax;
        double root2tip = options.root_to_tip;
        int nsites = options.nsites;

        for (int rep = 0; rep < options.num_replicates; rep++)
        {
            //Simulate tree and characters
            phylo<basic_newick> simT;
            switch (testnum)
            {
            case 0: //uniform
                simT = sim_uniform_phylo<basic_newick>(ntax); break;
            case 1: //yule
                simT = coalescent<basic_newick>(ntax, 1.0); break;
            case 2: //Caterpillar
                simT = sim_caterpillar_phylo<basic_newick>(ntax); break;
            case 3: //Beta critical
                simT = sim_beta_critical_tree<basic_newick>(ntax, true); break;
            default: assert(false);
            };

            double currentRootToTip = root_to_tip<basic_newick>(simT.root());
            scale(simT, root2tip / currentRootToTip);
            double treeLength = computeLength(simT);
            assignBranchIDs(simT); //Numbers edges - uses negative numbers for internal edges.
            vector<string> taxa_names;
            for (int i = 0; i < ntax; i++)
                taxa_names.push_back("T" + to_string(i));

            vector<sequence> alignment(ntax);

            JCModel model;
            simSequences(simT, model, nsites, alignment);
            vector<pair<Pattern, int>> patterns = sortPatterns(alignment);
            long nPatterns = patterns.size();

            int pruningHeight = 0;
            double pruningTime = 0.0;
            double pruning_avheight = 0.0;

            if (options.runMCMC_PruningTree)
            {
                phylo<DecomNodeData> decomTree;
                constructPruningDecomp(simT, decomTree);
                pruningHeight = nodeHeight<DecomNodeData>(decomTree.root());
                pruning_avheight = averageNodeHeight<DecomNodeData>(decomTree.root());
                phylo<DecomNodeDataAllSites> pruningT;
                copy(decomTree, pruningT);

                seed_random(42);
                MCMCResults pruningResults = runDecomMCMC(pruningT, model, patterns, options);
                pruningTime = accumulate(pruningResults.iterTime.begin(), pruningResults.iterTime.end(), 0.0);
                decomTree.clear();
                pruningT.clear();
            }

            int decomHeight = 0;
            double decomTime = 0.0;
            double decom_avheight = 0.0;

            if (options.runMCMC_DecomTree)
            {
                phylo<DecomNodeData> decomTree;
                constructDecompTree(simT, decomTree, false);
                decomHeight = nodeHeight<DecomNodeData>(decomTree.root());
                decom_avheight = averageNodeHeight<DecomNodeData>(decomTree.root());
                phylo<DecomNodeDataAllSites> lvdT;
                copy(decomTree, lvdT);

                seed_random(42);
                MCMCResults decomResults = runDecomMCMC(lvdT, model, patterns, options);
                decomTime = accumulate(decomResults.iterTime.begin(), decomResults.iterTime.end(), 0.0);
                decomTree.clear();
                lvdT.clear();
            }


            //Output
            cout << test_names[testnum] << "\t"
                << ntax << "\t"
                << root2tip << "\t"
                << treeLength << "\t"
                << nsites << "\t"
                << nPatterns << "\t"
                << pruningHeight << "\t"
                << pruning_avheight << "\t"
                << decomHeight << "\t"
                << decom_avheight << "\t";

            if (options.runMCMC_PruningTree)
                cout << pruningTime << "\t";
            if (options.runMCMC_DecomTree)
                cout << decomTime << "\t";
            if (options.runMCMC_PruningTree && options.runMCMC_DecomTree)
                cout << (pruningTime / decomTime) << "\t";
            cout << endl;

            if (options.outputTrees)
            {
                print_newick<basic_newick>(cerr, simT, taxa_names);
                cerr << ";" << endl;
            }
            simT.clear();
        }
    }
    PartialLikelihoodTensor::printCounts(cerr);
    return (0);
}
