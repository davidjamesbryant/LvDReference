/*
 *
 *  Created by David Bryant on 28/03/25.
 *
 */


#include <unistd.h>
#include <ctime>
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

class SimulationOptions {
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
    bool outputPruningTime = true;
    bool outputLvDTime = true;
    bool outputCompressionCount = false;
    bool outputHeader = true;
    bool outputTrees = false;
    int num_replicates;
    int num_iterations = 1000;
    double prior_branch_rate = 20.0;
    double proposal_width    = 0.01;
};

std::ostream& operator<<(std::ostream& os, const SimulationOptions& opt) {
    // Save current formatting flags so we can restore them at the end
    std::ios::fmtflags f = os.flags();
    os << std::boolalpha;  // print bools as true/false
    
    os << "SimulationOptions {\n"
    << "  ntaxa = " << opt.ntax << "\n"
    << "  test_uniform = " << opt.test_uniform << "\n"
    << "  test_yule = " << opt.test_yule << "\n"
    << "  test_caterpillar = " << opt.test_caterpillar << "\n"
    << "  test_beta_crit = " << opt.test_beta_crit << "\n"
    << "  root_to_tip = " << opt.root_to_tip << "\n"
    << "  nsites = " << opt.nsites << "\n"
    << "  outputPruningTime = " << opt.outputPruningTime << "\n"
    << "  outputLvDTime = " << opt.outputLvDTime << "\n"
    << "  ouputCompressionCount = " << opt.outputCompressionCount << "\n"
    << "  outputHeader = " << opt.outputHeader << "\n"
    << "  num_replicates = " << opt.num_replicates << "\n"
    << "  num_iterations = " << opt.num_iterations << "\n"
    << "  prior_branch_rate = " << opt.prior_branch_rate << "\n"
    << "  proposal_width = " << opt.proposal_width << "\n"
    << "  ouputTrees = " << opt.outputTrees << "\n"
    << "}";
    
    os.flags(f);  // restore flags
    return os;
}


string Usage() {
    string s;
    s = "RunSimulation -BCUYclhp -n <ntax> -s <nsites> -r <height> -R <replicates> [-i <iterations>] [-P <rate>]\n\n";
    s += "\t -BCUYclhpt  (flags)\n";
    s += "\t\tB:\t Simulate from beta critical distribution\n";
    s += "\t\tC:\t Simulate from random caterpillar\n";
    s += "\t\tU:\t Simulate from uniform distribution\n";
    s += "\t\tY:\t Simulate from Yule distribution\n";
    s += "\t<At least one tree distribution must be selected>\n";
    s += "\t\tc:\t Output the total number of compression subpatterns\n";
    s += "\t\tl:\t Output the likelihood and timing for the LvD algorithm\n";
    s += "\t\th:\t Output the header at top\n";
    s += "\t\tp:\t Output the likelihood and timing for the pruning algorithm\n";
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

bool is_int(const std::string& s, int& value) {
    try {
        size_t pos;
        int v = std::stoi(s, &pos, 10);
        if (pos != s.size()) return false; // leftover chars
        value = v;
        return true;
    } catch (...) {
        return false;  // invalid or out of range
    }
}


bool is_double(const std::string& s, double& value) {
    try {
        value = std::atof(s.c_str());
        return true;
    } catch (...) {
        return false;  // invalid or out of range
    }
}


bool parseArguments(int argc, char* argv[], SimulationOptions& options, string& errmsg) {
    if (argc < 2) {
        errmsg = "ERROR: incorrect number of arguments\n\n" + Usage();
        return false;
    }

    // First argument must be the flags string
    string s = string(argv[1]);
    if (s[0] != '-') {
        errmsg = "ERROR: First argument must be the option string\n\n" + Usage();
        return false;
    }
    options.test_uniform            = s.find('U') != string::npos;
    options.test_caterpillar        = s.find('C') != string::npos;
    options.test_beta_crit          = s.find('B') != string::npos;
    options.test_yule               = s.find('Y') != string::npos;
    options.outputCompressionCount  = s.find('c') != string::npos;
    options.outputPruningTime       = s.find('p') != string::npos;
    options.outputHeader            = s.find('h') != string::npos;
    options.outputLvDTime           = s.find('l') != string::npos;
    options.outputTrees             = s.find('t') != string::npos;

    // Track which required args have been seen
    bool have_n = false, have_s = false, have_r = false, have_R = false;

    // Parse remaining key-value pairs
    for (int i = 2; i < argc; i++) {
        string key = string(argv[i]);
        if (i + 1 >= argc) {
            errmsg = "ERROR: missing value for " + key + "\n\n" + Usage();
            return false;
        }
        string val = string(argv[++i]);

        if (key == "-n") {
            if (!is_int(val, options.ntax))   { errmsg = "ERROR: invalid value for -n\n\n"  + Usage(); return false; }
            have_n = true;
        } else if (key == "-s") {
            if (!is_int(val, options.nsites))  { errmsg = "ERROR: invalid value for -s\n\n"  + Usage(); return false; }
            have_s = true;
        } else if (key == "-r") {
            if (!is_double(val, options.root_to_tip)) { errmsg = "ERROR: invalid value for -r\n\n" + Usage(); return false; }
            have_r = true;
        } else if (key == "-R") {
            if (!is_int(val, options.num_replicates)) { errmsg = "ERROR: invalid value for -R\n\n" + Usage(); return false; }
            have_R = true;
        } else if (key == "-i") {
            if (!is_int(val, options.num_iterations))    { errmsg = "ERROR: invalid value for -i\n\n" + Usage(); return false; }
        } else if (key == "-P") {
            if (!is_double(val, options.prior_branch_rate)) { errmsg = "ERROR: invalid value for -P\n\n" + Usage(); return false; }
        } else if (key == "-w") {
            if (!is_double(val, options.proposal_width) || options.proposal_width <= 0) { errmsg = "ERROR: invalid value for -u\n\n" + Usage(); return false; }
        } else {
            errmsg = "ERROR: unrecognised option " + key + "\n\n" + Usage();
            return false;
        }
    }

    if (!have_n) { errmsg = "ERROR: -n <ntax> is required\n\n"       + Usage(); return false; }
    if (!have_s) { errmsg = "ERROR: -s <nsites> is required\n\n"     + Usage(); return false; }
    if (!have_r) { errmsg = "ERROR: -r <height> is required\n\n"     + Usage(); return false; }
    if (!have_R) { errmsg = "ERROR: -R <replicates> is required\n\n" + Usage(); return false; }

    return true;
}


/**
 Code for simulating sequences on a tree
 */
class simNodeData :public basic_newick  {
public:
    simNodeData() : basic_newick() {
        state=0;
    }
    simNodeData(const basic_newick& node) : basic_newick(node) {
        state=0;
    }
    int state;
    vector< vector<double> > transitionMat;
};

void simSequences(phylo<basic_newick>& T, SubstModel& model, int nsites, vector<sequence>& alignment) {
    //Set up the tree
    phylo<simNodeData> simT;
    copy<basic_newick,simNodeData>(T,simT);
    
    for(auto p=simT.leftmost_leaf();!p.root();p=p.next_post()) {
        p->transitionMat.resize(4);
        for (int i=0;i<model.num_states();i++) {
            p->transitionMat[i].resize(4);
            for (int j=0;j<model.num_states();j++)
                p->transitionMat[i][j] = model.Pij(i,j,p->length);
        }
    }
    vector<double> pi_vec(4);
    for(int i=0;i<4;i++)
        pi_vec[i]=model.pi(i);
    
    //Simulate characters
    int ntax = alignment.size();
    for(int i=0;i<ntax;i++)
        alignment[i].resize(nsites);
    
    for(int i=0;i<nsites;i++) {
        for(auto p = simT.root(); !p.null();p=p.next_pre()) {
            if (p.root())
                p->state = random_discrete(pi_vec);
            else {
                p->state = random_discrete(p->transitionMat[p.par()->state]);
            }
            if (p.leaf())
                alignment[p->id][i] = p->state;
        }
    }
    simT.clear();
}

/**
 Given an alignment and a subset of the set of taxa, this calculates the number of distinct patterns for the alignment restricted to taxa in that subset. It does this using a lexicographic sort.
 */
long countPatterns( vector<sequence>& alignment,  vector<bool>& subset) {
    
    int ntax = alignment.size();
    
    //Count size of subset
    int n = std::count(subset.begin(), subset.end(), true);
    
    if (n==0)
        return 0; //This happens at unlabelled leaves
    
    vector<sequence> sub_alignment;
    for(int i=0;i<ntax;i++) {
        if (subset[i]) {
            sub_alignment.push_back(alignment[i]);
        }
    }
    auto patterns = compressPatterns(sub_alignment,lex);
    
    long count = patterns.size();
    patterns.clear();
    sub_alignment.clear();
    return count;
}

/**
 Compute the number of distinct patterns for the alignment restricted to each clade in the tree, summed over all clades.
 
 After being called, the parameter clade returns the set of taxa at or below p.
 */
long countAllPatterns( vector<sequence>& alignment, phylo<DecomNodeData>::iterator p, vector<bool>& clade) {
    long count = 0;
    int ntax = alignment.size();
    clade.resize(ntax);
    fill(clade.begin(),clade.end(),false);
    
    for(phylo<DecomNodeData>::iterator q = p.left();!q.null();q = q.right()) {
        vector<bool> subclade(ntax);
        count+= countAllPatterns(alignment,q,subclade);
        for(int i=0;i<ntax;i++)
            clade[i] = clade[i]||subclade[i];
    }
    if (p.leaf() && p->isClade)
        clade[p->id] = true;
    
    long cladeCount = countPatterns(alignment,clade);
    count+=cladeCount;
    p->meta_data = to_string(cladeCount)+","+to_string(count);
    
    
    
    return count;
}

/**
 Wrapper function which doesn't return clade
 */
long countAllPatterns(vector<sequence>& alignment, phylo<DecomNodeData>& tree) {
    vector<bool> clade;
    long count = countAllPatterns(alignment,tree.root(),clade);
    clade.clear();
    return count;
}

// Assign stable branch IDs: leaves keep ids 0..ntax-1; internal non-root nodes
// get ids -(ntax), -(ntax+1), ... so that abs(p->id) indexes branches 0..numBranches-1.
static void assignBranchIDs(phylo<basic_newick>& tree) {
    int ntax = 0;
    for (auto p = tree.leftmost_leaf(); !p.null(); p = p.next_post())
        if (p.leaf()) ntax++;
    int nextId = ntax;
    for (auto p = tree.leftmost_leaf(); !p.null(); p = p.next_post())
        if (!p.leaf() && !p.root())
            p->id = -nextId++;
}

double pathLengthToLeftmostLeaf(phylo<basic_newick>& tree) {
    double length = 0.0;
    auto p = tree.root().left();
    while(!p.null()) {
        length += p->length;
        p=p.left();
    }
    return length;
}


int main(int argc, char* argv[]) {
    
    SimulationOptions options;
    string err_msg;
    
    if (!parseArguments(argc,argv,options,err_msg)) {
        cerr<<err_msg<<endl;
        exit(1);
    }
    
    //cerr<<"Options:\n"<<options<<endl;
    
    std::string test_names[] = {"uniform", "yule", "caterpillar", "beta_crit"};
    bool run_test[] = {options.test_uniform, options.test_yule, options.test_caterpillar, options.test_beta_crit};
    
    JCModel model;
    
    ofstream treefile;
    
    //Output header
    if (options.outputHeader) {
        cout<<"Distribution\t"
        <<"ntax\t"
        <<"root to tip\t"
        <<"total tree length\t"
        <<"average dist\t"
        <<"nsites\t"
        <<"nPatterns\t"
        <<"preprocessing time\t"
        <<"height pruning tree\t"
        <<"height LvD tree\t";
        if (options.outputCompressionCount)
            cout<<"num subpatterns pruning\t"
            <<"num subpatterns LvD\t";
        if (options.outputPruningTime)
            cout<<"logL pruning\t"
            <<"pruning time\t";
        if (options.outputLvDTime)
            cout<<"logL LvD\t"
            <<"decom time\t";
        if (options.outputPruningTime&&options.outputLvDTime)
            cout<<"ratio\t";
        cout<<endl;
    }
    
    
    for (int testnum =0;testnum<4;testnum++) {
        if (!run_test[testnum])
            continue;
        
        
        int ntax = options.ntax;
        
        vector<sequence> alignment(ntax);
        
        double root2tip = options.root_to_tip;
        int nsites = options.nsites;
        for (int rep=0;rep<options.num_replicates;rep++) {
            phylo<basic_newick> simT;
            switch (testnum) {
                case 0:  //uniform
                    simT = sim_uniform_phylo<basic_newick>(ntax);
                    break;
                case 1: //yule
                    simT = coalescent<basic_newick>(ntax,1.0);
                    break;
                case 2: //Caterpillar
                    simT =  sim_caterpillar_phylo<basic_newick>(ntax);
                    break;
                case 3: //Beta critical
                    simT =  sim_beta_critical_tree<basic_newick>(ntax,true);
                    break;
            };
            double currentRootToTip = pathLengthToLeftmostLeaf(simT);
            scale(simT,root2tip/currentRootToTip);
            double treeLength = computeLength(simT);
            simSequences(simT,model,nsites,alignment);
            assignBranchIDs(simT);

            //OUTPUT
            vector<string> taxa_names;
            for(int i=0;i<ntax;i++)
                taxa_names.push_back("T"+to_string(i));

            Stopwatch stopwatch;

            long tourLength = 0;
            stopwatch.reset();
            vector<pair<Pattern, int>> patterns = tspPatternSort(alignment,tourLength);
            double siteSortingTime = stopwatch.get();

            auto differences = computePatternDifferences(patterns);
            long nPatterns = patterns.size();

            phylo<DecomNodeData> decomTree;
            constructPruningDecomp(simT, decomTree);
            int pruningHeight = nodeHeight<DecomNodeData>(decomTree.root());

            double compressionCountPruning = 0.0;
            if (options.outputCompressionCount)
                compressionCountPruning = (double)countAllPatterns(alignment,decomTree);

            vector<double> patternL;
            double Lpruning = 0.0, pruningTime = 0.0;
            if (options.outputPruningTime) {
                // Initial log-likelihood
                Stopwatch sw0;
                Lpruning = computeLikelihood(decomTree, model, patterns, differences, patternL, sw0);

                // MCMC on pruning decomp tree
                phylo<DecomNodeDataAllSites> pruningT;
                copy(decomTree, pruningT);
                {
                    Stopwatch sw;
                    vector<double> mcmcPatternL;
                    Scalar logLik = computeLikelihood(pruningT, model, patterns, mcmcPatternL, sw);

                    Scalar logPrior = 0.0;
                    for (auto p = pruningT.leftmost_leaf(); !p.null(); p = p.next_post())
                        if (p.leaf())
                            logPrior += log(options.prior_branch_rate) - options.prior_branch_rate * p->length;
                    Scalar logPost = logLik + logPrior;

                    int numBranches = 0;
                    for (auto p = pruningT.leftmost_leaf(); !p.null(); p = p.next_post())
                        if (p.leaf()) numBranches++;
                    vector<phylo<DecomNodeDataAllSites>::iterator> branches(numBranches);
                    for (auto p = pruningT.leftmost_leaf(); !p.null(); p = p.next_post())
                        if (p.leaf())
                            branches[abs(p->id)] = p;

                    seed_random(42);
                    for (int iter = 0; iter < options.num_iterations; ++iter) {
                        int branchIdx = (int)random_num((unsigned int)numBranches);
                        auto bp = branches[branchIdx];
                        double oldLen = bp->length;
                        double newLen = oldLen + randu(-options.proposal_width, options.proposal_width);
                        if (newLen > 0.0) {
                            auto t0 = chrono::steady_clock::now();
                            Scalar newLogLik = updateBranchLength(pruningT, model, patterns, mcmcPatternL, bp, newLen);
                            double elapsed = chrono::duration<double>(chrono::steady_clock::now() - t0).count();
                            Scalar newLogPrior = logPrior + options.prior_branch_rate * (oldLen - newLen);
                            Scalar newLogPost  = newLogLik + newLogPrior;
                            if (log(randu()) < newLogPost - logPost) {
                                logPrior = newLogPrior;
                                logPost  = newLogPost;
                            } else {
                                auto t1 = chrono::steady_clock::now();
                                updateBranchLength(pruningT, model, patterns, mcmcPatternL, bp, oldLen);
                                elapsed += chrono::duration<double>(chrono::steady_clock::now() - t1).count();
                            }
                            pruningTime += elapsed;
                        }
                    }
                }
            }

            decomTree.clear();
            constructDecompTree(simT, decomTree,false);
            int decomHeight = nodeHeight<DecomNodeData>(decomTree.root());

            double compressionCountDecom = 0.0;
            if (options.outputCompressionCount)
                compressionCountDecom = (double)countAllPatterns(alignment,decomTree);

            double Ldecom = 0.0, decomTime = 0.0;
            if (options.outputLvDTime) {
                // Initial log-likelihood
                Stopwatch sw0;
                Ldecom = computeLikelihood(decomTree, model, patterns, differences, patternL, sw0);

                // MCMC on LvD decomp tree
                phylo<DecomNodeDataAllSites> lvdT;
                copy(decomTree, lvdT);
                {
                    Stopwatch sw;
                    vector<double> mcmcPatternL;
                    Scalar logLik = computeLikelihood(lvdT, model, patterns, mcmcPatternL, sw);

                    Scalar logPrior = 0.0;
                    for (auto p = lvdT.leftmost_leaf(); !p.null(); p = p.next_post())
                        if (p.leaf())
                            logPrior += log(options.prior_branch_rate) - options.prior_branch_rate * p->length;
                    Scalar logPost = logLik + logPrior;

                    int numBranches = 0;
                    for (auto p = lvdT.leftmost_leaf(); !p.null(); p = p.next_post())
                        if (p.leaf()) numBranches++;
                    vector<phylo<DecomNodeDataAllSites>::iterator> branches(numBranches);
                    for (auto p = lvdT.leftmost_leaf(); !p.null(); p = p.next_post())
                        if (p.leaf())
                            branches[abs(p->id)] = p;

                    seed_random(42);
                    for (int iter = 0; iter < options.num_iterations; ++iter) {
                        int branchIdx = (int)random_num((unsigned int)numBranches);
                        auto bp = branches[branchIdx];
                        double oldLen = bp->length;
                        double newLen = oldLen + randu(-options.proposal_width, options.proposal_width);
                        if (newLen > 0.0) {
                            auto t0 = chrono::steady_clock::now();
                            Scalar newLogLik = updateBranchLength(lvdT, model, patterns, mcmcPatternL, bp, newLen);
                            double elapsed = chrono::duration<double>(chrono::steady_clock::now() - t0).count();
                            Scalar newLogPrior = logPrior + options.prior_branch_rate * (oldLen - newLen);
                            Scalar newLogPost  = newLogLik + newLogPrior;
                            if (log(randu()) < newLogPost - logPost) {
                                logPrior = newLogPrior;
                                logPost  = newLogPost;
                            } else {
                                auto t1 = chrono::steady_clock::now();
                                updateBranchLength(lvdT, model, patterns, mcmcPatternL, bp, oldLen);
                                elapsed += chrono::duration<double>(chrono::steady_clock::now() - t1).count();
                            }
                            decomTime += elapsed;
                        }
                    }
                }
            }
            
            
            
            //Output
            cout<<test_names[testnum]<<"\t"
            <<ntax<<"\t"
            <<root2tip<<"\t"
            <<treeLength<<"\t"
            <<(double)tourLength/(nsites-1)<<"\t"
            <<nsites<<"\t"
            <<nPatterns<<"\t"
            <<siteSortingTime<<"\t"
            <<pruningHeight<<"\t"
            <<decomHeight<<"\t";
            if (options.outputCompressionCount)
                cout<<compressionCountPruning<<"\t"
                <<compressionCountDecom<<"\t";
            if (options.outputPruningTime)
                cout<<Lpruning<<"\t"
                <<pruningTime<<"\t";
            if (options.outputLvDTime)
                cout<<Ldecom<<"\t"
                <<decomTime<<"\t";
            if (options.outputPruningTime&&options.outputLvDTime)
                cout<<(pruningTime/decomTime)<<"\t";
            cout<<endl;
            
            if (options.outputTrees) {
                print_newick<basic_newick>(cerr,simT,taxa_names);
                cerr<<";"<<endl;
            }
            
            decomTree.clear();
            simT.clear();
        }
    }
    return (0);
}

