/*
 *
 *  Created by David Bryant on 28/03/25.
 *
 */


#include <unistd.h>
#include <ctime>
#include <list>
#include <charconv>
#include "Phylib/phylib.h"
#include "standardLikelihood.h"
#include "DecompositionTree/decompositionTree.h"

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
    << "  ouputTrees = " << opt.outputTrees << "\n"
    << "}";
    
    os.flags(f);  // restore flags
    return os;
}


string Usage() {
    string s;
    s = "RunSimulation -BCUYclhp  -n <ntax> -s <nsites> -r <height> -R <number of replicates> [-t <filename>]\n\n";
    s +="\t -BCUYcdhp  (flags)\n";
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
    s += "\t -n  <ntax>\t Number of taxa\n";
    s += "\t -s  <nsites>\t Number of sites\n";
    s += "\t -r <height>\t Expected number of mutations from root to a tip\n";
    s += "\t -R <number>\t Number of replicates\n";
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
    if (argc != 10) {
        errmsg = "ERROR: incorrect number of arguments\n\n"+Usage();
        return false;
    }
    //Parse option string
    string s = string(argv[1]);
    if (s[0]!='-') {
        errmsg = "ERROR: First argument must be the option string\n\n"+Usage();
        return false;
    }
    options.test_uniform = s.find('U')!= std::string::npos;
    options.test_caterpillar = s.find('C')!= std::string::npos;
    options.test_beta_crit = s.find('B')!= std::string::npos;
    options.test_yule = s.find('Y')!= std::string::npos;
    
    options.outputCompressionCount = s.find('c')!= std::string::npos;
    options.outputPruningTime = s.find('p')!= std::string::npos;
    options.outputHeader = s.find('h')!= std::string::npos;
    options.outputLvDTime = s.find('l')!= std::string::npos;
    options.outputTrees = s.find('t')!= std::string::npos;
    
    
    //Number of taxa
    if (string(argv[2])!="-n" || !is_int(string(argv[3]),options.ntax)) {
        errmsg = "ERROR reading number of taxa\n\n"+Usage();
        return false;
    }
    
    //Number of sites
    if (string(argv[4])!="-s" || !is_int(string(argv[5]),options.nsites)) {
        errmsg = "ERROR reading number of sites\n\n"+Usage();
        return false;
    }
    
    //Number of sites
    if (string(argv[6])!="-r" || !is_double(string(argv[7]),options.root_to_tip)) {
        errmsg = "ERROR reading root to tip length\n\n"+Usage();
        return false;
    }
    
    //Number of replicates
    if (string(argv[8])!="-R" || !is_int(string(argv[9]),options.num_replicates)) {
        errmsg = "ERROR reading root to tip length\n\n"+Usage();
        return false;
    }
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
            
            double compressionCountPruning;
            if (options.outputCompressionCount)
                compressionCountPruning =  (double)countAllPatterns(alignment,decomTree);
            
            stopwatch.reset();
            vector<Scalar> patternL;
            double Lpruning, pruningTime;
            if (options.outputPruningTime) {
                Lpruning = computeLikelihood(decomTree, model, patterns, differences, patternL,stopwatch);
                pruningTime = stopwatch.get();
            }
            
            decomTree.clear();
            constructDecompTree(simT, decomTree,false);
            int decomHeight = nodeHeight<DecomNodeData>(decomTree.root());
            
            double compressionCountDecom;
            if (options.outputCompressionCount)
                compressionCountDecom =  (double)countAllPatterns(alignment,decomTree);
            
            stopwatch.reset();
            double Ldecom, decomTime;
            if (options.outputLvDTime) {
                Ldecom = computeLikelihood(decomTree, model, patterns, differences, patternL,stopwatch);
                decomTime = stopwatch.get();
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

