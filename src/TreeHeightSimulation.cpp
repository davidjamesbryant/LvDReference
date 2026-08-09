/*
 *
 *  Created by David Bryant on 28/03/25.
 *
 */


#include <unistd.h>
#include "Phylib/phylib.h"
#include "DecompositionTree/decompositionTree.h"

using namespace Phylib;


/**
 Parsing the program arguments
 */

class SimulationOptions
{
public:
    //Tree distribution parameters
    int min_ntax = 10;
    int max_ntax = 10000; //Number of taxa
    bool test_uniform = false; //Uniform on trees
    bool test_yule = false; //Yule = coalescent
    bool test_caterpillar = false; //Same as uniform
    bool test_beta_crit = false; //Beta critical


    //What to output

    bool outputHeader = true;

    int num_replicates;

};

std::ostream& operator<<(std::ostream& os, const SimulationOptions& opt)
{
    // Save current formatting flags so we can restore them at the end
    std::ios::fmtflags f = os.flags();
    os << std::boolalpha; // print bools as true/false

    os << "SimulationOptions {\n"
        << "  min_ntaxa = " << opt.min_ntax << "\n"
        << "  max_ntaxa = " << opt.max_ntax << "\n"
        << "  test_uniform = " << opt.test_uniform << "\n"
        << "  test_yule = " << opt.test_yule << "\n"
        << "  test_caterpillar = " << opt.test_caterpillar << "\n"
        << "  test_beta_crit = " << opt.test_beta_crit << "\n"
        << "  outputHeader = " << opt.outputHeader << "\n"
        << "  num_replicates = " << opt.num_replicates << "\n"
        << "}";

    os.flags(f); // restore flags
    return os;
}


string Usage()
{
    string s;
    s =
        "TreeHeightSimulation -BCUYh -m <min_ntax> -n <max_ntax> -R <replicates>\n\n";
    s += "\t -BCUYh  (flags)\n";
    s += "\t\tB:\t Simulate from beta critical distribution\n";
    s += "\t\tC:\t Simulate from random caterpillar\n";
    s += "\t\tU:\t Simulate from uniform distribution\n";
    s += "\t\tY:\t Simulate from Yule distribution\n";
    s += "\t<At least one tree distribution must be selected>\n";
    s += "\t\th:\t Output the header at top\n";
    s += "\t -m  <min_ntax>\t Minimum number of taxa (required) \n";
    s += "\t -n  <max_ntax>\t Maximum number of taxa (required)\n";
    s += "\t -R  <number>\t Number of replicates (required)\n";
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
    options.outputHeader = s.find('h') != string::npos;

    // Track which required args have been seen
    bool have_m = false, have_n = false, have_R = false;

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

        if (key == "-m")
        {
            if (!is_int(val, options.min_ntax))
            {
                errmsg = "ERROR: invalid value for -m\n\n" + Usage();
                return false;
            }
            have_m = true;
        }
        else if (key == "-n")
        {
            if (!is_int(val, options.max_ntax))
            {
                errmsg = "ERROR: invalid value for -n\n\n" + Usage();
                return false;
            }
            have_n = true;
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
        else
        {
            errmsg = "ERROR: unrecognised option " + key + "\n\n" + Usage();
            return false;
        }
    }

    if (!have_m)
    {
        errmsg = "ERROR: -m <min_ntax> is required\n\n" + Usage();
        return false;
    }
    if (!have_n)
    {
        errmsg = "ERROR: -n <max_ntax> is required\n\n" + Usage();
        return false;
    }
    if (!have_R)
    {
        errmsg = "ERROR: -R <replicates> is required\n\n" + Usage();
        return false;
    }
    if (options.min_ntax < 1)
    {
        errmsg = "ERROR: -m <min_ntax> must be at least 1\n\n" + Usage();
        return false;
    }
    if (options.min_ntax > options.max_ntax)
    {
        errmsg = "ERROR: -m <min_ntax> must not exceed -n <max_ntax>\n\n" + Usage();
        return false;
    }

    return true;
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

    //Output header
    if (options.outputHeader)
    {
        cout << "Distribution\t"
            << "ntax\t"
            << "height pruning tree\t"
            << "av height pruning tree\t"
            << "height LvD tree\t"
            << "av height LvD tree";
        cout << endl;
    }


    for (int testnum = 0; testnum < 4; testnum++)
    {
        if (!run_test[testnum])
            continue;

        int min_ntax = options.min_ntax;
        int max_ntax = options.max_ntax;

        for (int rep = 0; rep < options.num_replicates; rep++)
        {
            //Simulate tree
            //Simulate ntax so that log(ntax) is uniform between min_ntax and max_ntax.
            int ntax = floor(pow(10, log10(min_ntax) + randu() * (log10(max_ntax) - log10(min_ntax))));
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

            vector<string> taxa_names;
            for (int i = 0; i < ntax; i++)
                taxa_names.push_back("T" + to_string(i));

            phylo<DecomNodeData> decomTree;
            constructPruningDecomp(simT, decomTree);
            int pruningHeight = nodeHeight<DecomNodeData>(decomTree.root());
            double pruning_avheight = averageNodeHeight<DecomNodeData>(decomTree.root());
            decomTree.clear();

            constructDecompTree(simT, decomTree, false);
            int decomHeight = nodeHeight<DecomNodeData>(decomTree.root());
            double decom_avheight = averageNodeHeight<DecomNodeData>(decomTree.root());
            decomTree.clear();


            //Output
            cout << test_names[testnum] << "\t"
                << ntax << "\t"
                << pruningHeight << "\t"
                << pruning_avheight << "\t"
                << decomHeight << "\t"
                << decom_avheight << "\t";
            cout << endl;

            simT.clear();
        }
    }
    return (0);
}
