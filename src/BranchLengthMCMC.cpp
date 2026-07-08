/*
 *
 *  Created by David Bryant on 08/05/26.
 *
 *  MCMC analysis of branch lengths using multiple likelihood methods,
 *  comparing standard pruning, pruning-decomposition, and LvD-decomposition.
 */

#include <unistd.h>
#include <ctime>
#include <fstream>
#include <iomanip>
#include <numeric>
#include "Phylib/phylib.h"
#include "standardLikelihood.h"
#include "DecompositionTree/decompositionTree.h"
#include "DecompositionTree/decompositionLikelihood.h"
#include "DecomMCMC.h"

using namespace Phylib;


// ─────────────────────────────────────────────────────────────────────────────
// Options
// ─────────────────────────────────────────────────────────────────────────────

class MCMCOptions
{
public:
    double prior_branch_rate = 20.0; // Rate λ of exponential prior on branch lengths
    double proposal_width = 0.01; // Half-width of uniform proposal distribution
    int num_iterations = 1000; // Number of MCMC iterations

    bool run_standard = false; // Standard likelihood (NodeDataAllSites)
    bool run_pruning = false; // Decomposition with pruning tree
    bool run_lvd = false; // Decomposition with LvD tree
    bool output_header = true;

    string seq_filename;
    string tree_filename;
};

ostream& operator<<(ostream& os, const MCMCOptions& opt)
{
    ios::fmtflags f = os.flags();
    os << boolalpha;
    os << "MCMCOptions {\n"
        << "  prior_branch_rate     = " << opt.prior_branch_rate << "\n"
        << "  proposal_width = " << opt.proposal_width << "\n"
        << "  num_iterations = " << opt.num_iterations << "\n"
        << "  run_standard   = " << opt.run_standard << "\n"
        << "  run_pruning    = " << opt.run_pruning << "\n"
        << "  run_lvd        = " << opt.run_lvd << "\n"
        << "  seq_filename   = " << opt.seq_filename << "\n"
        << "  tree_filename  = " << opt.tree_filename << "\n"
        << "}";
    os.flags(f);
    return os;
}


// ─────────────────────────────────────────────────────────────────────────────
// Argument parsing
// ─────────────────────────────────────────────────────────────────────────────

static string Usage()
{
    string s;
    s = "BranchLengthMCMC [-SPLh] -p <prior> -u <proposal> -i <iterations> <seqfile> <treefile>\n\n";
    s += "\t Flags (combined as first argument):\n";
    s += "\t\t S:  Run MCMC using standard likelihood (NodeDataAllSites)\n";
    s += "\t\t P:  Run MCMC using decomposition tree with pruning layout\n";
    s += "\t\t L:  Run MCMC using decomposition tree with LvD layout\n";
    s += "\t\t h:  Output header row in output files\n";
    s += "\t <At least one of S, P, L must be selected>\n\n";
    s += "\t -p <prior>       Rate parameter λ for exponential prior on branch lengths (default 20)\n";
    s += "\t -u <proposal>    Half-width of uniform proposal distribution (default 0.01)\n";
    s += "\t -i <iterations>  Number of MCMC iterations (positive integer - default 1000)\n";
    s += "\t <seqfile>        PHYLIP-format sequence file\n";
    s += "\t <treefile>       Newick-format tree file\n";
    return s;
}

static bool is_int(const string& s, int& value)
{
    try
    {
        size_t pos;
        int v = stoi(s, &pos, 10);
        if (pos != s.size()) return false;
        value = v;
        return true;
    }
    catch (...)
    {
        return false;
    }
}

static bool is_double(const string& s, double& value)
{
    try
    {
        size_t pos;
        value = stod(s, &pos);
        return pos == s.size();
    }
    catch (...)
    {
        return false;
    }
}

static bool parseArguments(int argc, char* argv[], MCMCOptions& options, string& errmsg)
{
    if (argc != 10)
    {
        errmsg = "ERROR: incorrect number of arguments\n\n" + Usage();
        return false;
    }

    string flags = string(argv[1]);
    if (flags[0] != '-')
    {
        errmsg = "ERROR: first argument must be the flag string (e.g. -SPh)\n\n" + Usage();
        return false;
    }
    options.run_standard = flags.find('S') != string::npos;
    options.run_pruning = flags.find('P') != string::npos;
    options.run_lvd = flags.find('L') != string::npos;
    options.output_header = flags.find('h') != string::npos;

    if (!options.run_standard && !options.run_pruning && !options.run_lvd)
    {
        errmsg = "ERROR: at least one of S, P, L must be selected\n\n" + Usage();
        return false;
    }

    if (string(argv[2]) != "-p" || !is_double(string(argv[3]), options.prior_branch_rate) || options.prior_branch_rate
        <= 0)
    {
        errmsg = "ERROR reading prior rate: -p must be followed by a positive number\n\n" + Usage();
        return false;
    }

    if (string(argv[4]) != "-u" || !is_double(string(argv[5]), options.proposal_width) || options.proposal_width <= 0)
    {
        errmsg = "ERROR reading proposal width: -u must be followed by a positive number\n\n" + Usage();
        return false;
    }

    if (string(argv[6]) != "-i" || !is_int(string(argv[7]), options.num_iterations) || options.num_iterations <= 0)
    {
        errmsg = "ERROR reading iteration count: -i must be followed by a positive integer\n\n" + Usage();
        return false;
    }

    options.seq_filename = string(argv[8]);
    options.tree_filename = string(argv[9]);
    return true;
}

static string get_stem(const string& path)
{
    size_t slash = path.find_last_of("/\\");
    string filename = (slash == string::npos) ? path : path.substr(slash + 1);
    size_t dot = filename.find_last_of('.');
    return (dot != string::npos) ? filename.substr(0, dot) : filename;
}


// ─────────────────────────────────────────────────────────────────────────────
// MCMC functions
// Each function runs num_iterations of a Metropolis-Hastings random-walk on
// all branch lengths, using an exponential(prior_branch_rate) prior and a uniform
// proposal of half-width proposal_width on a single randomly chosen branch.
// ─────────────────────────────────────────────────────────────────────────────

static MCMCResults runStandardMCMC(
    phylo<basic_newick>& tree,
    const SubstModel& model,
    const vector<pair<Pattern, int>>& patterns,
    const MCMCOptions& options)
{
    // Fixed seed — all MCMC variants seed identically so proposals are directly comparable
    seed_random(42);

    // Build a NodeDataAllSites tree from the newick topology
    phylo<NodeDataAllSites> t;
    copy(tree, t);

    // Compute initial log-likelihood (also allocates and fills all partial arrays)
    Stopwatch sw;
    Scalar logLik = computeLikelihood(t, model, patterns, sw);

    // Initial log-prior: sum_branches [ log(rate) - rate * length ]
    Scalar logPrior = 0.0;
    for (auto p = t.leftmost_leaf(); !p.null(); p = p.next_post())
        if (!p.root())
            logPrior += log(options.prior_branch_rate) - options.prior_branch_rate * p->length;
    Scalar logPost = logLik + logPrior;


    int max_id = 0;
    for (auto p = t.leftmost_leaf(); !p.null(); p = p.next_post())
        if (!p.root())
            max_id = std::max(max_id, abs(p->id));
    vector<phylo<NodeDataAllSites>::iterator> branches(max_id+1);

    for (auto p = t.leftmost_leaf(); !p.null(); p = p.next_post())
        if (!p.root())
            branches[abs(p->id)] = p;

    //Remove null pointers - this can happen if not all taxa appear in the tree
    int numBranches=0;
    for (int i=0;i<=max_id; i++)
    {
        if (!branches[i].null())
        {
            branches[numBranches] = branches[i];
            numBranches++;
        }
    }
    branches.resize(numBranches);



    MCMCResults results;
    results.logPosterior.reserve(options.num_iterations);
    results.iterTime.reserve(options.num_iterations);
    results.branch.reserve(options.num_iterations);
    results.old_length.reserve(options.num_iterations);
    results.new_length.reserve(options.num_iterations);
    results.accepted.reserve(options.num_iterations);
    results.initialLogPost = logPost;
    results.numBranches = numBranches;

    for (int iter = 0; iter < options.num_iterations; ++iter)
    {
        int branchIdx = (int)random_num((unsigned int)numBranches);
        auto p = branches[branchIdx];
        double oldLen = p->length;
        double newLen = oldLen + randu(-options.proposal_width, options.proposal_width);

        bool accept = false;
        double likeTime = 0.0;
        if (newLen > 0.0)
        {
            auto t0 = chrono::steady_clock::now();
            Scalar newLogLik = updateBranchLength(t, model, patterns, p, newLen);
            likeTime += chrono::duration<double>(chrono::steady_clock::now() - t0).count();

            // Prior ratio for exponential: log p(new) - log p(old) = -rate*(new - old)
            Scalar newLogPrior = logPrior + options.prior_branch_rate * (oldLen - newLen);
            Scalar newLogPost = newLogLik + newLogPrior;

            // Metropolis-Hastings; proposal is symmetric so no Hastings correction needed
            accept = (log(randu()) < newLogPost - logPost);

            if (accept)
            {
                logLik = newLogLik;
                logPrior = newLogPrior;
                logPost = newLogPost;
            }
            else
            {
                auto t1 = chrono::steady_clock::now();
                updateBranchLength(t, model, patterns, p, oldLen); // restore
                likeTime += chrono::duration<double>(chrono::steady_clock::now() - t1).count();
            }
        }
        else
        {
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

static MCMCResults runDecomMCMC(
    phylo<basic_newick>& tree,
    const SubstModel& model,
    const vector<pair<Pattern, int>>& patterns,
    const MCMCOptions& options,
    int& treeHeight,
    double& avTreeHeight,
    bool useLvD)
{
    seed_random(42);

    phylo<DecomNodeData> decomBase;
    if (useLvD)
        constructDecompTree(tree, decomBase, false);
    else
        constructPruningDecomp(tree, decomBase);

    phylo<DecomNodeDataAllSites> t;
    copy(decomBase, t);

    treeHeight = nodeHeight<DecomNodeData>(decomBase.root());
    avTreeHeight = averageNodeHeight<DecomNodeData>(decomBase.root());

    return runDecomMCMC(t, model, patterns, options.prior_branch_rate, options.proposal_width, options.num_iterations);
}

// ─────────────────────────────────────────────────────────────────────────────
// main
// ─────────────────────────────────────────────────────────────────────────────

int main(int argc, char* argv[])
{
    MCMCOptions options;
    string errmsg;
    if (!parseArguments(argc, argv, options, errmsg))
    {
        cerr << errmsg << endl;
        exit(1);
    }

    if (options.output_header)
    {
        cout << "Alignment\t"
            << "nsites\t"
            << "nPatterns\t"
            << "nIterations\t"
            << "Treefile\t"
            << "ntax\t";
        if (options.run_standard)
            cout << "standard time\t";
        if (options.run_pruning)
        {
            cout << "height pruning tree\t"
                << "av_height pruning tree\t"
                << "pruning time\t";
        }
        if (options.run_lvd)
        {
            cout << "height LvD tree\t"
                << "av_height LvD tree\t"
                << "LvD time\t";
        }
        if (options.run_pruning && options.run_lvd)
            cout<<"ratio";
        cout << endl;
    }


    // Read sequence file
    vector<string> taxa_names;
    vector<sequence> alignment;
    {
        ifstream seqFile(options.seq_filename);
        if (!seqFile)
        {
            cerr << "Error opening sequence file: " << options.seq_filename << "\n";
            exit(2);
        }
        try
        {
            if (!read_phylip_seqs(seqFile, taxa_names, alignment))
            {
                cerr << "Error reading sequence file: " << options.seq_filename << "\n";
                exit(3);
            }
        }
        catch (const PhylibException& e)
        {
            cerr << "Error reading sequence file:\n" << e.getMessage() << "\n";
            exit(4);
        }
    }

    // Compress sites into patterns
    vector<pair<Pattern, int>> patterns = sortPatterns(alignment);

    // Read tree file (first tree only)
    phylo<basic_newick> newickTree;

    ifstream treeFile(options.tree_filename);
    if (!treeFile)
    {
        cerr << "Error opening tree file: " << options.tree_filename << "\n";
        exit(5);
    }
    try
    {
        read_newick(treeFile, newickTree, taxa_names, 0.0, true);
    }
    catch (const PhylibException& e)
    {
        cerr << "Error reading tree file:\n" << e.getMessage() << "\n";
        exit(6);
    }
    if (newickTree.empty())
    {
        cerr << "No tree found in tree file: " << options.tree_filename << "\n";
        exit(7);
    }

    if (!checkIfBinary(newickTree))
        resolveTree(newickTree);

    assignBranchIDs(newickTree);

    // Run selected MCMC analyses
    JCModel model;


    MCMCResults standardResults, pruningResults, lvdResults;
    int pruning_height, lvd_height;
    double pruning_avheight, lvd_avheight;
    pruning_height = lvd_height = 0;
    pruning_avheight = lvd_avheight = 0.0;


    if (options.run_standard)
        standardResults = runStandardMCMC(newickTree, model, patterns, options);
    if (options.run_pruning)
        pruningResults = runDecomMCMC(newickTree, model, patterns, options, pruning_height, pruning_avheight, false);
    if (options.run_lvd)
        lvdResults = runDecomMCMC(newickTree, model, patterns, options, lvd_height, lvd_avheight, true);

    // ─── Diagnostics ─────────────────────────────────────────────────────────
    cout << fixed << setprecision(4);
    cout << get_stem(options.seq_filename)<<"\t"   //seq name
    <<alignment[0].size()<<"\t"    //nsites
    <<patterns.size()<<"\t"     //npatterns
    <<options.num_iterations<<"\t"
    <<get_stem(options.tree_filename)<<"\t"
    <<taxa_names.size()<<"\t";
    if (options.run_standard)
         cout<<accumulate(standardResults.iterTime.begin(),standardResults.iterTime.end(), 0.0) <<"\t";
    double pruning_time = 0.0, lvd_time = 0.0;
    if (options.run_pruning)
        {
            cout << pruning_height<<"\t"<<pruning_avheight<<"\t";
            pruning_time = accumulate(pruningResults.iterTime.begin(),pruningResults.iterTime.end(), 0.0);
            cout<< pruning_time <<"\t";
        }
    if (options.run_lvd)
        {
            cout << lvd_height<<"\t"<<lvd_avheight<<"\t";
            lvd_time = accumulate(lvdResults.iterTime.begin(),lvdResults.iterTime.end(), 0.0);
            cout<< lvd_time <<"\t";
        }
    if (options.run_pruning && options.run_lvd)
        cout<< pruning_time / lvd_time;
    cout << endl;

    return 0;
}
