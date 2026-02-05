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
    os << appName << "\n\nReads in a sequence and a file of trees, then outputs a tidied up version of the treefile with comments to stderr\n";
    os << "Usage:\n\n" << appName << " <seqfile> <treefile> \n";
    os << endl;
}

string get_stem(const string &path) {
    // Find last directory separator
    size_t slash = path.find_last_of("/\\");
    string filename = (slash == string::npos)
                           ? path
                           : path.substr(slash + 1);

    // Remove final extension
    size_t dot = filename.find_last_of('.');
    if (dot != string::npos)
        return filename.substr(0, dot);

    return filename;
}


static bool handleInput(int argc, char* argv[], bool& printHeader, string& seqFilename, string& treeFilename)
{
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

// Helper: remove punctuation, keep everything else as-is
string removePunctuation(const string& s) {
    string out;
    out.reserve(s.size());
    for (unsigned char ch : s) {
        if (!ispunct(ch)) {
            out.push_back(ch);
        }
    }
    return out;
}

// Helper: lowercase version (ASCII)
string toLowerAscii(const string& s) {
    string out;
    out.reserve(s.size());
    for (unsigned char ch : s) {
        out.push_back(static_cast<char>(tolower(ch)));
    }
    return out;
}

// "Hamming-like" distance that allows different lengths
size_t hammingLikeDistance(const string& a, const string& b) {
    const size_t n = a.size();
    const size_t m = b.size();
    const size_t common = (n < m ? n : m);
    size_t dist = 0;

    for (size_t i = 0; i < common; ++i) {
        if (a[i] != b[i]) {
            ++dist;
        }
    }
    dist += (n > m ? n - m : m - n); // penalty for length difference
    return dist;
}

enum MatchType {
    Exact = 0,
    PunctInsensitive = 1,
    CaseInsensitive = 2,
    Hamming = 3
};


// Return index of best match in candidates, or -1 if candidates is empty.
int findClosestMatchIndex(const string& query,
                          const vector<string>& candidates,
                          MatchType& bestType) {
    if (candidates.empty()) return -1;

    

    const string queryNoPunct = removePunctuation(query);
    const string queryLower   = toLowerAscii(query);

    bestType = Hamming;
    int bestDist = numeric_limits<int>::max();
    int bestIdx = -1;

    for (int i = 0; i < candidates.size(); ++i) {
        const string& cand = candidates[i];

        // 1. Exact match
        if (cand == query) {
            // This is the best possible; you can return immediately
            bestType = Hamming;
            return i;
        }

        // 2. Differs only in punctuation
        string candNoPunct = removePunctuation(cand);
        if (candNoPunct == queryNoPunct) {
            if (bestType > PunctInsensitive) {
                bestType = PunctInsensitive;
                bestDist = 0;
                bestIdx = i;
            }
            continue; // no need to check case/Hamming
        }

        // 3. Differs only by case
        string candLower = toLowerAscii(cand);
        if (candLower == queryLower) {
            if (bestType > CaseInsensitive) {
                bestType = CaseInsensitive;
                bestDist = 0;
                bestIdx = i;
            }
            continue; // no need to check Hamming
        }

        // 4. Fallback: Hamming-like distance
        int dist = hammingLikeDistance(query, cand);
        if (bestType > Hamming || (bestType == Hamming && dist < bestDist)) {
            bestType = Hamming;
            bestDist = dist;
            bestIdx = i;
        }
    }

    return bestIdx;
}

int main(int argc, char* argv[]) {
    
    string appName = "TidyTree";
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
    
    while (!done) {
        vector<string> taxaInTree;
        try {
            read_newick(treeFile, newickTree, taxaInTree, 0.0, false);
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
        
        //Clear out internal taxon labels
        for(auto p = newickTree.root();!p.null();p=p.next_pre()) {
            if (p->id >= 0 && !p.leaf()) {
                cerr<<"Removed internal label "+taxaInTree[p->id]<<endl;
                p->id = -1;
            }
        }
        
        //Match tree taxa names with sequence taxa names
        string MatchNames[] = {"Exact","Punctuation","Upper-lower case","Approx match"};

        cerr<<"Tree number "<<tree_num<<" Matching tree taxa with sequence taxa\ntaxon in tree\tsequence\tMismatch type"<<endl;
        for(auto p = newickTree.root();!p.null();p=p.next_pre()) {
            if (p.leaf()) {
                string taxon_tree_name = taxaInTree[p->id];
                MatchType matchType;
                int idx = findClosestMatchIndex(taxon_tree_name, taxa_names, matchType);
                if (idx<0) {
                    cerr<<"Cannot match taxon  "<<taxon_tree_name<<endl;
                    exit(1);
                }
                string closest_taxon = taxa_names[idx];
                if (closest_taxon!=taxon_tree_name)
                    cerr<<taxon_tree_name<<"\t"<<closest_taxon<<"\t"<<MatchNames[(int)matchType]<<endl;
                p->id = idx;
            }
        }
        
        //Output tree to stdout, followed by ";". Doesn't output metadata
       
        print_newick(cout,newickTree,taxa_names,true,false);
        cout<<";"<<endl;
    }
    
    treeFile.close();
    
    return(0);
}
