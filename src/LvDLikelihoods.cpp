/*
 *
 *  Created by David Bryant on 28/03/25.
 *
 */


//#DEFINE VERBOSE_SSR

#include <unistd.h>
#include <ctime>
#include <list>
#include "Phylib/phylib.h"
#include "standardLikelihood.h"
#include "DecompositionTree/decompositionTree.h"

using namespace Phylib;

static void printUsage(ostream& os, const string& appName) {
    os << appName << "\n\nCalculates likelihoods and running times for trees in a tree file\n";
    os << "Usage:\n\n" << appName << " [-h] <seqfile> <treefile> \n";
    os << endl;
}

string get_stem(const std::string &path) {
    // Find last directory separator
    size_t slash = path.find_last_of("/\\");
    std::string filename = (slash == std::string::npos)
                           ? path
                           : path.substr(slash + 1);

    // Remove final extension
    size_t dot = filename.find_last_of('.');
    if (dot != std::string::npos)
        return filename.substr(0, dot);

    return filename;
}


static bool handleInput(int argc, char* argv[], bool& printHeader, string& seqFilename, string& treeFilename)
{
    if (argc == 4) {
        if (string(argv[1])!="-h")
            return false;
        printHeader = true;
        seqFilename = string(argv[2]);
        treeFilename = string(argv[3]);
        return true;
    }
    if (argc==3) {
        printHeader = true;
        seqFilename = string(argv[1]);
        treeFilename = string(argv[2]);
        return true;
    }
    return false;
    
    if (argc <= 1) {
        // No filenames provided; still considered "parsed OK".
        return true;
    }
}



int main(int argc, char* argv[]) {
    
    string appName = "LvDLikelihoods";
    phylo<basic_newick> newickTree;
    vector<string> taxa_names;
    vector<sequence> alignment;

    ifstream is_tree;
    vector<string> filenames;
    bool printHeader;
    string seqFilename, treeFilename;
    
    bool OK = handleInput(argc, argv, printHeader, seqFilename, treeFilename);
    
    if (!OK)
        exit(1);
    
    if (printHeader) {
        cout<<"Alignment\t"
        <<"nsites\t"
        <<"nPatterns\t"
        <<"preprocessing time\t"
        <<"average site dist\t"
        <<"average pattern dist\t"
        <<"Treefile\t"
        <<"Number\t"
        <<"ntax\t"
        <<"root to tip\t"
        <<"total tree length\t"
        <<"height pruning tree\t"
        <<"height LvD tree\t"
        <<"height ratio\t"
        <<"logL pruning\t"
        <<"pruning time\t"
        <<"logL LvD\t"
        <<"decom time\t"
        <<"time ratio"
        <<endl;
    }

    //Read in the sequence file
    ifstream seqFile;
    seqFile.open(seqFilename.c_str());
    if (!seqFile) {
        cerr << "Error reading in sequence file " << seqFilename << "\n\n";
        printUsage(cerr, appName);
        exit(2);
    }
    
    try {
    if (!read_phylip_seqs(seqFile, taxa_names, alignment)) {
        cerr << "Problem reading sequence file";
        return (3);
    }
    } catch (const PhylibException& e) {
        cerr<<"Error reading sequence file:\n"<<e.getMessage()<<endl;
        exit(4);
    }
    seqFile.close();
    

    //Read in the tree file
    ifstream treeFile;
    treeFile.open(treeFilename.c_str());
    if (!treeFile) {
        cerr << "Error reading in tree file " << treeFilename << "\n\n";
        printUsage(cerr, appName);
        exit(5);
    }
    
    //Extract short filenames
    string shortSeqName =  get_stem(seqFilename);
    string shortTreeName = get_stem(treeFilename);
    
    int tree_num=0;
    bool done = false;
    JCModel model;
    Stopwatch stopwatch;
    
    //Preprocess sequences
    int nsites=alignment[0].size();
    long tourLength = 0;
    stopwatch.reset();
    vector<pair<Pattern, int>> patterns = tspPatternSort(alignment,tourLength);
    double siteSortingTime = stopwatch.get();
    
    auto differences = computePatternDifferences(patterns);
    long nPatterns = patterns.size();
    
    
    
    while (!done) {
        try {
            read_newick(treeFile, newickTree, taxa_names, 0.0, true);
        } catch (const PhylibException& e) {
            cerr<<"EXCEPTION: "<<e.getMessage()<<endl;
            done = true;
            exit(-1);
        }
        if (newickTree.empty()) {
            done = true;  // End of file
            break;
        }
        tree_num++;
        
        if (!checkIfBinary(newickTree)) {
            resolveTree(newickTree);
        }
        
        double rootToTip = root_to_tip<basic_newick>(newickTree.root());
        double treeLength = computeLength<basic_newick>(newickTree);
        
        
        int ntax = taxa_names.size();
        int orig_height = nodeHeight<basic_newick>(newickTree.root());
        Stopwatch stopwatch;
        
        cerr<<"Pruning tree"<<endl;
        
        phylo<DecomNodeData> decomTree;
        constructPruningDecomp(newickTree, decomTree);
        int pruning_height = nodeHeight<DecomNodeData>(decomTree.root());
        

        stopwatch.reset();
        vector<Scalar> patternL;
        double Lpruning = computeLikelihood(decomTree, model, patterns, differences, patternL,stopwatch);
        double pruningTime = stopwatch.get();
        

        decomTree.clear();
        constructDecompTree(newickTree, decomTree, false);
        
        check_topology<DecomNodeData>(decomTree);
        
        
        int decomp_height = nodeHeight<DecomNodeData>(decomTree.root());
        double height_ratio = ((double)pruning_height)/decomp_height;
        
        stopwatch.reset();
        double Ldecom = computeLikelihood(decomTree, model, patterns, differences, patternL,stopwatch);
        double decomTime = stopwatch.get();
        
        //Output
        cout<<shortSeqName<<"\t"
        <<nsites<<"\t"
        <<nPatterns<<"\t"
        <<siteSortingTime<<"\t"
        <<(double)tourLength/(nsites-1)<<"\t"
        <<(double)tourLength/(nPatterns-1)<<"\t"
        <<shortTreeName<<"\t"
        <<tree_num<<"\t"
        <<ntax<<"\t"
        <<rootToTip<<"\t"
        <<treeLength<<"\t"
        <<pruning_height<<"\t"
        <<decomp_height<<"\t"
        <<height_ratio<<"\t"
        <<Lpruning<<"\t"
        <<pruningTime<<"\t"
        <<Ldecom<<"\t"
        <<decomTime<<"\t"
        <<((double)pruningTime/decomTime)
        <<endl;
    }
    
    treeFile.close();
    
    return(0);
}
