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
    os << appName << "\n\nCalculates heights of trees and decomposition trees for trees in a tree file\n";
    os << "Usage:\n\n" << appName << " [-h] <treefile> <treefile> .. <treefile>\n";
    os <<  "\t -h\t output header";
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

static bool handleInput(int argc, char* argv[], bool& printHeader, vector<string>& filenames)
{
    filenames.clear();
    if (argc <= 1) {
        // No filenames provided; still considered "parsed OK".
        return true;
    }
    
    int firstFile = 1;
    if (string(argv[1])=="-h") {
        printHeader = true;
        firstFile++;
    } else
        printHeader = false;
        

    filenames.reserve(static_cast<size_t>(argc - 1));
    for (int i = firstFile; i < argc; ++i) {
        if (argv[i] && argv[i][0] != '\0') {
            filenames.emplace_back(argv[i]);
        }
    }
    return true;
}



int main(int argc, char* argv[]) {
    
    phylo<basic_newick> newickTree;
    vector<string> taxa_names;

    ifstream is_tree;
    vector<string> filenames;
    bool printHeader;
    
    bool OK = handleInput(argc, argv, printHeader, filenames);
    if (!OK)
        exit(1);
    
    if (printHeader) {
        cout<<"Filename\t"
        <<"Number\t"
        <<"ntax\t"
        <<"root to tip\t"
        <<"total tree length\t"
        <<"height pruning tree\t"
        <<"height LvD tree\t"
        <<"height ratio\t"
        <<endl;
    }
    //Print header
    
    for(size_t filenum =0; filenum<filenames.size();filenum++) {
        is_tree.open(filenames[filenum].c_str());
        if (!is_tree) {
            cerr << "Error reading in tree file " << filenames[filenum] << "\n\n";
            continue;
        }
        
        //Extract short filename
        
        string shortFileName =  get_stem(filenames[filenum]);
        
        
        int tree_num=0;
        bool done = false;
        
        while (!done) {
            taxa_names.clear();
            try {
                read_newick(is_tree, newickTree, taxa_names, 0.0, false);
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
            
            //print_newick<basic_newick>(cout, newickTree, taxa_names, true, false);
            if (!checkIfBinary(newickTree)) {
                resolveTree(newickTree);
            }
            
            double rootToTip = root_to_tip<basic_newick>(newickTree.root());
            double treeLength = computeLength<basic_newick>(newickTree);
            
            
            int ntax = taxa_names.size();
            int orig_height = nodeHeight<basic_newick>(newickTree.root());
            
            
            phylo<DecomNodeData> decomTree;
            constructPruningDecomp(newickTree, decomTree);
            int pruning_height = nodeHeight<DecomNodeData>(decomTree.root());
            constructDecompTree(newickTree, decomTree, false);
            int decomp_height = nodeHeight<DecomNodeData>(decomTree.root());
            double ratio = ((double)pruning_height)/decomp_height;
            
            //Output
            cout<<shortFileName<<"\t"
            <<tree_num<<"\t"
            <<ntax<<"\t"
            <<rootToTip<<"\t"
            <<treeLength<<"\t"
            <<pruning_height<<"\t"
            <<decomp_height<<"\t"
            <<ratio<<"\t"
            <<endl;
        }
        is_tree.close();
    }
    
    return(0);
}
