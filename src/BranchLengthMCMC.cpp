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
#include "Phylib/phylib.h"
#include "standardLikelihood.h"
#include "DecompositionTree/decompositionTree.h"

using namespace Phylib;


// ─────────────────────────────────────────────────────────────────────────────
// Options
// ─────────────────────────────────────────────────────────────────────────────

class MCMCOptions {
public:
    double prior_rate     = 1.0;    // Rate λ of exponential prior on branch lengths
    double proposal_width = 0.01;   // Half-width of uniform proposal distribution
    int    num_iterations = 1000;   // Number of MCMC iterations

    bool run_standard = false;   // Standard likelihood (NodeDataAllSites)
    bool run_pruning  = false;   // Decomposition with pruning tree
    bool run_lvd      = false;   // Decomposition with LvD tree
    bool output_header = true;

    string seq_filename;
    string tree_filename;
};

ostream& operator<<(ostream& os, const MCMCOptions& opt) {
    ios::fmtflags f = os.flags();
    os << boolalpha;
    os << "MCMCOptions {\n"
       << "  prior_rate     = " << opt.prior_rate     << "\n"
       << "  proposal_width = " << opt.proposal_width << "\n"
       << "  num_iterations = " << opt.num_iterations << "\n"
       << "  run_standard   = " << opt.run_standard   << "\n"
       << "  run_pruning    = " << opt.run_pruning    << "\n"
       << "  run_lvd        = " << opt.run_lvd        << "\n"
       << "  seq_filename   = " << opt.seq_filename   << "\n"
       << "  tree_filename  = " << opt.tree_filename  << "\n"
       << "}";
    os.flags(f);
    return os;
}


// ─────────────────────────────────────────────────────────────────────────────
// Argument parsing
// ─────────────────────────────────────────────────────────────────────────────

static string Usage() {
    string s;
    s  = "BranchLengthMCMC [-SPLh] -p <prior> -u <proposal> -i <iterations> <seqfile> <treefile>\n\n";
    s += "\t Flags (combined as first argument):\n";
    s += "\t\t S:  Run MCMC using standard likelihood (NodeDataAllSites)\n";
    s += "\t\t P:  Run MCMC using decomposition tree with pruning layout\n";
    s += "\t\t L:  Run MCMC using decomposition tree with LvD layout\n";
    s += "\t\t h:  Output header row in output files\n";
    s += "\t <At least one of S, P, L must be selected>\n\n";
    s += "\t -p <prior>       Rate parameter λ for exponential prior on branch lengths\n";
    s += "\t -u <proposal>    Half-width of uniform proposal distribution\n";
    s += "\t -i <iterations>  Number of MCMC iterations (positive integer)\n";
    s += "\t <seqfile>        PHYLIP-format sequence file\n";
    s += "\t <treefile>       Newick-format tree file\n";
    return s;
}

static bool is_int(const string& s, int& value) {
    try {
        size_t pos;
        int v = stoi(s, &pos, 10);
        if (pos != s.size()) return false;
        value = v;
        return true;
    } catch (...) {
        return false;
    }
}

static bool is_double(const string& s, double& value) {
    try {
        size_t pos;
        value = stod(s, &pos);
        return pos == s.size();
    } catch (...) {
        return false;
    }
}

static bool parseArguments(int argc, char* argv[], MCMCOptions& options, string& errmsg) {
    // Expected layout:
    //   argv[1]  flag string  e.g. -SPh
    //   argv[2]  -p   argv[3]  prior_rate
    //   argv[4]  -u   argv[5]  proposal_width
    //   argv[6]  -i   argv[7]  num_iterations
    //   argv[8]  seqfile
    //   argv[9]  treefile
    if (argc != 10) {
        errmsg = "ERROR: incorrect number of arguments\n\n" + Usage();
        return false;
    }

    string flags = string(argv[1]);
    if (flags[0] != '-') {
        errmsg = "ERROR: first argument must be the flag string (e.g. -SPh)\n\n" + Usage();
        return false;
    }
    options.run_standard  = flags.find('S') != string::npos;
    options.run_pruning   = flags.find('P') != string::npos;
    options.run_lvd       = flags.find('L') != string::npos;
    options.output_header = flags.find('h') != string::npos;

    if (!options.run_standard && !options.run_pruning && !options.run_lvd) {
        errmsg = "ERROR: at least one of S, P, L must be selected\n\n" + Usage();
        return false;
    }

    if (string(argv[2]) != "-p" || !is_double(string(argv[3]), options.prior_rate) || options.prior_rate <= 0) {
        errmsg = "ERROR reading prior rate: -p must be followed by a positive number\n\n" + Usage();
        return false;
    }

    if (string(argv[4]) != "-u" || !is_double(string(argv[5]), options.proposal_width) || options.proposal_width <= 0) {
        errmsg = "ERROR reading proposal width: -u must be followed by a positive number\n\n" + Usage();
        return false;
    }

    if (string(argv[6]) != "-i" || !is_int(string(argv[7]), options.num_iterations) || options.num_iterations <= 0) {
        errmsg = "ERROR reading iteration count: -i must be followed by a positive integer\n\n" + Usage();
        return false;
    }

    options.seq_filename  = string(argv[8]);
    options.tree_filename = string(argv[9]);
    return true;
}

static string get_stem(const string& path) {
    size_t slash = path.find_last_of("/\\");
    string filename = (slash == string::npos) ? path : path.substr(slash + 1);
    size_t dot = filename.find_last_of('.');
    return (dot != string::npos) ? filename.substr(0, dot) : filename;
}


// ─────────────────────────────────────────────────────────────────────────────
// Results
// ─────────────────────────────────────────────────────────────────────────────

struct MCMCResults {
    vector<double> logPosterior;  // log-posterior at each accepted/proposed state
    vector<double> iterTime;      // wall-clock seconds consumed by each iteration
};


// ─────────────────────────────────────────────────────────────────────────────
// MCMC stubs
// Each function runs num_iterations of a Metropolis-Hastings random-walk on
// all branch lengths, using an exponential(prior_rate) prior and a uniform
// proposal of half-width proposal_width on a single randomly chosen branch.
// ─────────────────────────────────────────────────────────────────────────────

static MCMCResults runStandardMCMC(
        phylo<basic_newick>&             tree,
        const SubstModel&                model,
        const vector<pair<Pattern,int>>& patterns,
        const MCMCOptions&               options)
{
    MCMCResults results;
    results.logPosterior.reserve(options.num_iterations);
    results.iterTime.reserve(options.num_iterations);
    // TODO: implement using NodeDataAllSites / updateBranchLength
    return results;
}

static MCMCResults runPruningMCMC(
        phylo<basic_newick>&             tree,
        const SubstModel&                model,
        const vector<pair<Pattern,int>>& patterns,
        const MCMCOptions&               options)
{
    MCMCResults results;
    results.logPosterior.reserve(options.num_iterations);
    results.iterTime.reserve(options.num_iterations);
    // TODO: implement using constructPruningDecomp / computeLikelihood(DecomNodeData)
    return results;
}

static MCMCResults runLvDMCMC(
        phylo<basic_newick>&             tree,
        const SubstModel&                model,
        const vector<pair<Pattern,int>>& patterns,
        const MCMCOptions&               options)
{
    MCMCResults results;
    results.logPosterior.reserve(options.num_iterations);
    results.iterTime.reserve(options.num_iterations);
    // TODO: implement using constructDecompTree / computeLikelihood(DecomNodeData)
    return results;
}


// ─────────────────────────────────────────────────────────────────────────────
// main
// ─────────────────────────────────────────────────────────────────────────────

int main(int argc, char* argv[]) {

    MCMCOptions options;
    string errmsg;
    if (!parseArguments(argc, argv, options, errmsg)) {
        cerr << errmsg << endl;
        exit(1);
    }

    // Read sequence file
    vector<string>   taxa_names;
    vector<sequence> alignment;
    {
        ifstream seqFile(options.seq_filename);
        if (!seqFile) {
            cerr << "Error opening sequence file: " << options.seq_filename << "\n";
            exit(2);
        }
        try {
            if (!read_phylip_seqs(seqFile, taxa_names, alignment)) {
                cerr << "Error reading sequence file: " << options.seq_filename << "\n";
                exit(3);
            }
        } catch (const PhylibException& e) {
            cerr << "Error reading sequence file:\n" << e.getMessage() << "\n";
            exit(4);
        }
    }

    // Compress sites into patterns (TSP-ordered to minimise inter-pattern distance)
    long tourLength = 0;
    vector<pair<Pattern,int>> patterns = tspPatternSort(alignment, tourLength);

    // Read tree file (first tree only)
    phylo<basic_newick> newickTree;
    {
        ifstream treeFile(options.tree_filename);
        if (!treeFile) {
            cerr << "Error opening tree file: " << options.tree_filename << "\n";
            exit(5);
        }
        try {
            read_newick(treeFile, newickTree, taxa_names, 0.0, true);
        } catch (const PhylibException& e) {
            cerr << "Error reading tree file:\n" << e.getMessage() << "\n";
            exit(6);
        }
        if (newickTree.empty()) {
            cerr << "No tree found in tree file: " << options.tree_filename << "\n";
            exit(7);
        }
    }
    if (!checkIfBinary(newickTree))
        resolveTree(newickTree);

    // Run selected MCMC analyses
    JCModel model;

    MCMCResults standardResults, pruningResults, lvdResults;
    if (options.run_standard)
        standardResults = runStandardMCMC(newickTree, model, patterns, options);
    if (options.run_pruning)
        pruningResults  = runPruningMCMC (newickTree, model, patterns, options);
    if (options.run_lvd)
        lvdResults      = runLvDMCMC    (newickTree, model, patterns, options);

    // ─── Output files ────────────────────────────────────────────────────────
    string stem = get_stem(options.seq_filename) + "_" + get_stem(options.tree_filename);

    // Posterior log-densities — one row per iteration, one column per active analysis
    {
        ofstream postFile(stem + "_posterior.tsv");
        if (options.output_header) {
            postFile << "iteration";
            if (options.run_standard) postFile << "\tstandard_logP";
            if (options.run_pruning)  postFile << "\tpruning_logP";
            if (options.run_lvd)      postFile << "\tlvd_logP";
            postFile << "\n";
        }
        for (int i = 0; i < options.num_iterations; i++) {
            postFile << (i + 1);
            if (options.run_standard) postFile << "\t" << standardResults.logPosterior[i];
            if (options.run_pruning)  postFile << "\t" << pruningResults.logPosterior[i];
            if (options.run_lvd)      postFile << "\t" << lvdResults.logPosterior[i];
            postFile << "\n";
        }
    }

    // Per-iteration wall-clock times — one row per iteration, one column per active analysis
    {
        ofstream timesFile(stem + "_times.tsv");
        if (options.output_header) {
            timesFile << "iteration";
            if (options.run_standard) timesFile << "\tstandard_time";
            if (options.run_pruning)  timesFile << "\tpruning_time";
            if (options.run_lvd)      timesFile << "\tlvd_time";
            timesFile << "\n";
        }
        for (int i = 0; i < options.num_iterations; i++) {
            timesFile << (i + 1);
            if (options.run_standard) timesFile << "\t" << standardResults.iterTime[i];
            if (options.run_pruning)  timesFile << "\t" << pruningResults.iterTime[i];
            if (options.run_lvd)      timesFile << "\t" << lvdResults.iterTime[i];
            timesFile << "\n";
        }
    }

    return 0;
}
