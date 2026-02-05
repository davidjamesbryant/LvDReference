/**
 * @file phylip_seq.cpp
 * @brief Code file for routines for reading and writing phylip sequence files
 * @author David Bryant
 * @version 1.0
 */


#include"phylip_seq.h"

namespace Phylib {

using namespace std;
/*****************************
 Sequence data
 *****************************/
/*
 typedef unsigned short int  base;
 const int bA = 0;
 const int bC = 1;
 const int bG = 2;
 const int bT = 3;
 const int bR = 4; //Purine
 const int bY = 5; //Pyrimadine
 const int bS = 6; //GC
 const int bW = 7; //AT
 const int bK = 8; //GT
 const int bM = 9; //AC
 const int bD = 10; //AGT
 const int bV = 11; //ACG
 const int bH = 12; //ACT
 const int bB = 13; //CGT
 const int bX = 14; //Missing or gap
 
 
 typedef vector<base> sequence;
 typedef vector<base> pattern;
 */


/********************
 write_name(ostream& os, const string& name)
 
 writes out first 10 letters of name to os. If name is shorter than 10 letters
 then prints out spaces to pad out to 10 letters.
 
 We use the C++ string type to handle names (rather than char[11]) so that
 in future, we won't need to limit names to 10 characters.
 
 UPDATE: just use whole string
 
 ********************/
static void write_name(ostream& os, const string& name) {
//    unsigned int i;
//    for (i=0;i<10;i++) {
//        if (i<name.length())
//            os<<name[i];
//        else
//            os<<" ";
//    }
    os<<name;
}

/********************
 read_name(istream& os, string& name)
 
 reads in a name of length at most 10 letters. If the name is
 more than 10 letters long then only the first 10 letters are
 returned.
 ********************/
static void read_name(istream& os, string& name) {
    //char s[11]; /* 10 characters plus terminating 0 */
    //cin.width(10);
    //os>>s;
    //cerr<<"Reading "<<s<<endl;
    //name = string(s);
    os>>name;
    //cerr<<"Reading "<<name<<endl;
}


/* Printb
 Prints out a single base to ostream os */
static void printb(ostream& os, base b) {
    char ch='X';
    
    switch (b) {
        case bA: ch='A'; break;
        case bC: ch='C';break;
        case bG: ch='G';break;
        case bT: ch='T';break;
            
        case bR: ch='R';break;
        case bY: ch='Y';break;
        case bS: ch='S'; break;
        case bW: ch='W';break;
        case bK: ch='K';break;
        case bM: ch='M'; break;
        
        case bD: ch='D'; break;
        case bV: ch='V'; break;
        case bH: ch='H'; break;
        case bB: ch='B'; break;
            
        case bX: ch='?';break;
    }
    os<<ch;
}

/* Output operator for a sequence */
ostream& operator<<(ostream& os, const sequence& seq) {
    unsigned int i;
    for (i=0;i<seq.size();i++)
        printb(os,seq[i]);
    return os;
}

/* Return array of 1s and 0s for resolving ambiguity codes*/
void resolve_base(base b, vector<Scalar>& bases) {
    if (b==bA)
        bases =  {1,0,0,0};
    else if (b==bC)
        bases =  {0,1,0,0};
    else if (b==bG)
        bases =  {0,0,1,0};
    else if (b==bT)
        bases =  {0,0,0,1};
    
    else if (b==bR)
        bases =  {1,0,1,0};
    else if (b==bY)
        bases =  {0,1,0,1};
    else if (b==bW)
        bases =  {1,0,0,1};
    else if (b==bS)
        bases =  {0,1,1,0};
    else if (b==bK)
        bases =  {0,0,1,1};
    else if (b==bM)
        bases =  {1,1,0,0};
    
    else if (b==bD)
        bases =  {1,0,1,1};
    else if (b==bV)
        bases =  {1,1,1,0};
    else if (b==bH)
        bases =  {1,1,0,1};
    else if (b==bB)
        bases =  {0,1,1,1};
    
    else //(b==bX)
        bases =  {1,1,1,1};
}


/* convert a string to a sequence */
static void string_to_sequence(const string& s, sequence& seq) {
    unsigned int i;
    base b;
    
    seq.clear();
    seq.resize(s.length());
    for (i=0;i<s.length();i++) {
        switch (toupper(s[i])) {
            case 'A': b=bA; break;
            case 'C': b=bC; break;
            case 'G': b=bG; break;
            case 'T': case 'U': b=bT; break;
            case 'R': b=bR; break;
            case 'Y': b=bY; break;
            case 'S': b=bS; break;
            case 'W': b=bW; break;
            case 'K': b=bK; break;
            case 'M': b=bM; break;
                
            case 'D': b=bD; break;
            case 'V': b=bV; break;
            case 'H': b=bH; break;
            case 'B': b=bB; break;
                
            case '?': case 'N': case 'X': case '-': b = bX; break;
            default:
                throw PhylibException((string)"Encountered the char '"+s[i]+"' when reading sequence");
        }
        seq[i]=b;
    }
}


/********
 read_phylip_seqs
 
 Crude parser of the sequence file.
 takes a possibly empty list of sequence names as input.
 For each name+sequence read in, it checks if name is already in
 taxa_names. If it is, the seq is inserted into the corresponding entry of seqs.
 If it is not, then name is added to the end of taxa_names and the sequence is
 added to the end of seqs.
 *********/
bool read_phylip_seqs(istream& is, vector<string>& taxa_names, vector<sequence>& seqs) {
    string s;
    int use_old_names;
    int num_taxa;
    unsigned int num_bases;
    unsigned int index;
    sequence seq;
    
    use_old_names = (!taxa_names.empty());
    
    //read in number of taxa and num bases
    is>>num_taxa; 
    is>>num_bases;
    
    index=0;
    while (is) {
        read_name(is,s);
        if (!is) break; //test for eof (only generated after a read)
        //find if the name is already read. If it is, replace seq. If not, add a new seq.
        //this is a safety check for when we read the tree before the sequence.
        //if this is taking a lot of time, we could make a map type for taxa_names (log n search)
        if (use_old_names) {
            for(index=0;index<taxa_names.size() && taxa_names[index]!=s;index++)
                if(taxa_names[index]==s)
                    break;
            if(index==taxa_names.size()) {
                throw PhylibException("Error:- label"+s+" in sequence file is not the name of an existing taxon");
                return false;
            }
        }
        else
            taxa_names.push_back(s);
        
        is>>s;
        string_to_sequence(s,seq);
        
        //cerr<<"Read sequence of length "<<seq.size()<<endl;
        
        if(use_old_names)
            seqs[index]=seq;
        else
            seqs.push_back(seq);
        
        if (seq.size()!=num_bases) {
            throw PhylibException("Sequence "+taxa_names.back()+" is the wrong length");
            return false;
        }
        
        if ((int)seqs.size()>=num_taxa)
        	break;	
    }
    return true;
}

bool read_phylip_seqs(const string& filename, vector<string>& taxa_names, vector<sequence>& seqs) {
    ifstream the_file(filename.c_str());
    if (the_file)
        return read_phylip_seqs(the_file,taxa_names, seqs);
    else
        throw PhylibException("Could not open file: "+filename);
    
    return false;
}


/****************
 write_phylip_seqs(ostream& os, const vector<string>& taxa_names, const vector<sequence>& seqs);
 
 outputs the sequences.
 ****************/
bool write_phylip_seqs(ostream& os, const vector<string>& taxa_names, const vector<sequence>& seqs) {
    
    unsigned int i;
    os<<taxa_names.size()<<"\t"<<seqs[0].size()<<"\n"; /* output name and sequence length */
    
    if (taxa_names.size()!=seqs.size())
        throw PhylibException("There is a different number of taxa than sequences");
    for (i=0;i<seqs.size();i++) {
        write_name(os,taxa_names[i]);
        os<<seqs[i]<<"\n";
    }
    return true;
}

bool write_phylip_seqs(const string& filename, const vector<string>& taxa_names, const vector<sequence>& seqs) {
    ofstream the_file(filename.c_str());
    if (the_file)
        return write_phylip_seqs(the_file,taxa_names, seqs);
    else
        throw PhylibException("Could not write to file: "+filename);
    
    return false;
}


/**********
 Routines for sorting and compressing site patterns
 ***********/

/*
 getUncompressedSitePatterns
 
 This functiojn takes a vector of sequences and returns a vector of site patterns, without eliminating
 duplicates etc.
 */
static PatternVector getUncompressedSitePatterns(const vector<sequence>& sequences) {
    size_t sequenceLength = sequences[0].size();
    
    PatternVector patterns;
    
    // Extract columns (site patterns)
    for (size_t col = 0; col < sequenceLength; ++col) {
        Pattern pattern;
        for (const auto& seq : sequences) {
            pattern.push_back(seq[col]);
        }
        patterns.emplace_back(pattern, 1);  // Initialize frequency count to 1
    }
    
    return patterns;
}


// Define a custom comparator for patterns. This is lexicographic. Returns true if a
//comes before b.
struct LexicographicalComparator {
    bool operator()(const Pattern& a, const Pattern& b) const {
        for (size_t i = 0; i < a.size() && i < b.size(); ++i) {
            if (a[i] != b[i]) {
                return a[i] < b[i];  // Normal order
            }
        }
        return false;
    }
};

// Define a custom comparator for patterns. This is lexicographic, but the ordering
// of the bases is reversed for every second digit. This should reduce the number of
//differences between adjacent patterns.
struct ZigzagComparator {
    bool operator()(const Pattern& a, const Pattern& b) const {
        for (size_t i = 0; i < a.size() && i < b.size(); ++i) {
            // Reverse the order of bases for every second digit
            if (i % 2 == 1) {
                if (a[i] != b[i]) {
                    return a[i] > b[i];  // Reverse order
                }
            } else {
                if (a[i] != b[i]) {
                    return a[i] < b[i];  // Normal order
                }
            }
        }
        return false;
    }
};


/**
 Take a vector of sequences and returns a vector of patterns with frequencies.
 By default, patterns are sorted lexicographically. If zigzag is true then
 they are sorted lexicogrqaphically but the ordering of even taxa is reversed.
 **/

static PatternVector sortPatterns(const vector<sequence>& sequences,bool zigzag=false) {
    size_t sequenceLength = sequences[0].size();
    
    //Copy the patterns onto a map with comparator depending on whether we zigzag
    // the rows or not.
    
    map<Pattern, int,ZigzagComparator> patternCountsZig;
    map<Pattern, int,LexicographicalComparator> patternCountsLex;
    
    for (size_t col = 0; col < sequenceLength; ++col) {
        Pattern pattern;
        for (const auto& seq : sequences) {
            pattern.push_back(seq[col]);
        }
        if (zigzag)
            patternCountsZig[pattern]++;  // Increment frequency count
        else
            patternCountsLex[pattern]++;
    }
    if (zigzag)
        return {patternCountsZig.begin(), patternCountsZig.end()};
    else
        return {patternCountsLex.begin(), patternCountsLex.end()};
}



/*************
 Use a TSP heuristic to sort the patterns.
 **************/

// Edge structure for Prim's algorithm
struct Edge {
    int from, to;
    double weight;
    bool operator<(const Edge& other) const {
        return weight < other.weight;
    }
};

/**
 Returns number of differences between Pattern a and Pattern b
 **/

static int HammingDist(const Pattern& a, const Pattern& b) {
    size_t n = a.size();
    int diffs = 0;
    for(size_t i=0;i<n;i++)
        if (a[i]!=b[i])
            diffs++;
    return diffs;
}


// Compute MST for patterns using Prim's algorithm
static vector<vector<int>> computeMST(const PatternVector& patterns) {
    int n=patterns.size();
    
    vector<bool> inTree(n, false);
    vector<double> minDist(n, numeric_limits<double>::max());
    vector<int> parent(n, -1);
    vector<vector<int>> tree(n); // adjacency list for MST
    
    minDist[0] = 0;
    
    for (int i = 0; i < n; ++i) {
        int u = -1;
        for (int v = 0; v < n; ++v) {
            if (!inTree[v] && (u == -1 || minDist[v] < minDist[u]))
                u = v;
        }
        
        inTree[u] = true;
        if (parent[u] != -1) {
            tree[u].push_back(parent[u]);
            tree[parent[u]].push_back(u);
        }
        
        for (int v = 0; v < n; ++v) {
            double dist = HammingDist(patterns[u].first,patterns[v].first);
            if (!inTree[v] && dist < minDist[v]) {
                minDist[v] = dist;
                parent[v] = u;
            }
        }
    }
    
    return tree;
}

// DFS traversal of MST to generate Euler tour
static void dfs(int u, const vector<vector<int>>& tree, vector<bool>& visited, vector<int>& tour) {
    visited[u] = true;
    tour.push_back(u);
    for (int v : tree[u]) {
        if (!visited[v]) {
            dfs(v, tree, visited, tour);
            tour.push_back(u); // backtrack
        }
    }
}

// Remove duplicates to form TSP path
static vector<int> tspApproximation(PatternVector patterns) {
    int n=patterns.size();
    auto tree = computeMST(patterns);
//    cerr<<"Computed MST"<<endl;
    
    vector<bool> visited(n, false);
    vector<int> eulerTour;
    dfs(0, tree, visited, eulerTour);
    
    
    
    // Shortcut repeated vertices to form TSP path
    vector<bool> seen(n, false);
    vector<int> tspPath;
    for (int v : eulerTour) {
        if (!seen[v]) {
            tspPath.push_back(v);
            seen[v] = true;
        }
    }
    
//    cerr<<"Tour is: "<<endl;
//    for(int i=0;i<n;i++)
//        cerr<<" "<<tspPath[i];
//    cerr<<endl;
    
        
        
    return tspPath;
}

PatternVector tspPatternSort(const vector<sequence>& sequences, long& length) {
    cerr<<"Sorting sites"<<endl;
    
    PatternVector lexSorted = sortPatterns(sequences,false);
    int nPatterns = lexSorted.size();
    vector<int> tspPath = tspApproximation(lexSorted);
    PatternVector tourPatterns(nPatterns);
    for(size_t i=0;i<tspPath.size();i++)
        tourPatterns[i] = lexSorted[tspPath[i]];
   
    length = 0;
    for(int i =1;i<tourPatterns.size();i++)
        length+=HammingDist(tourPatterns[i-1].first,tourPatterns[i].first);
    
    cerr<<"Finished Sorting sites"<<endl;

    return tourPatterns;
}


/****
 Extract patterns from the alignment and then sort them according to one of three methods:
 uncompressed: no compression
 lex: lexicographically
 zigzag: lexicographically, but with bases in alternate taxa treated in reverse order
 tsp: sorted using a 2-approximate tsp heuristic.
 ***/
PatternVector compressPatterns(const vector<sequence>& sequences, PatternSorter method) {
    long length;
    switch (method) {
        case uncompressed:
            return getUncompressedSitePatterns(sequences);
        case lex:
            return sortPatterns(sequences,false);
        case zigzag:
            return sortPatterns(sequences,true);
        case tsp:
            return tspPatternSort(sequences,length);
    };
}


/****
 Take a vector of Patterns and return, for each pattern, a vector of changes from the
 previous one. We define all bases in the first pattern to have changed.
 This means that we can update likelihoods only looking at the taxa which have changed
 from one site to the next.
 ****/
PatternDiffs computePatternDifferences(const vector<pair<Pattern, int>>& patterns) {
    PatternDiffs differences;
    
    //First pattern - differences in all positions.
    vector<pair<int, unsigned short>> first_diff;
    const auto& firstPattern = patterns[0].first;
    for (size_t j = 0; j < firstPattern.size(); ++j) {
        first_diff.emplace_back(j, firstPattern[j]);  // Store only position and new value
    }
    differences.push_back(first_diff);
    
    /* // Output first pattern
     cout << "First pattern: ";
     for (const auto& p : first_diff) {
     cout << "(" << p.first << ", " << p.second << ") ";
     }
     cout << endl;
     */
    
    
    for (size_t i = 1; i < patterns.size(); i++) {
        const auto& prevPattern = patterns[i - 1].first;
        const auto& currPattern = patterns[i].first;
        vector<pair<int, unsigned short>> diff;
        
        for (size_t j = 0; j < prevPattern.size(); ++j) {
            if (prevPattern[j] != currPattern[j]) {
                diff.emplace_back(j, currPattern[j]);  // Store only position and new value
            }
        }
        
        differences.push_back(diff);
    }
    
    return differences;
}

}

