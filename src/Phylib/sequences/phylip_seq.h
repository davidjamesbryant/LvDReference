/**
 * @file phylip_seq.h
 * @brief Header file for routines for reading and writing phylip sequence files
 * @author David Bryant
 * @version 1.0
 */


#ifndef PHYLIP_SEQ
#define PHYLIP_SEQ



#include"../global/stdIncludes.h"
#include"../utilities/phylibException.h"
#include"phylib.h" // Include phylib.h first to define Scalar

namespace Phylib {
using namespace std;

/*****************************
 Sequence data
 *****************************/
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

//TODO Do all ambiguities.




enum PatternSorter {uncompressed,lex,zigzag,tsp};

typedef vector<base> sequence;
typedef vector<base> Pattern;
typedef vector<pair<Pattern,int> > PatternVector;
typedef vector<vector<pair<int, unsigned short>>> PatternDiffs;


/* Output operator for a sequence */
ostream& operator<<(ostream& os, const sequence& seq);

/* Return array of 1s and 0s for resolving ambiguity codes*/
void resolve_base(base b, vector<Scalar>& bases);

/********
 read_phylip_seqs
 *********/
bool read_phylip_seqs(istream& is, vector<string>& taxa_names, vector<sequence>& seqs);

bool read_phylip_seqs(const string& filename, vector<string>& taxa_names, vector<sequence>& seqs);

/****************
 write_phylip_seqs(ostream& os, const vector<string>& taxa_names, const vector<sequence>& seqs);
 
 outputs the sequences.
 ****************/
bool write_phylip_seqs(ostream& os, const vector<string>& taxa_names, const vector<sequence>& seqs);

bool write_phylip_seqs(const string& filename, const vector<string>& taxa_names, const vector<sequence>& seqs);


/*************
 Extract a vector of patterns with their frequencies.
 ****/

PatternVector compressPatterns(const vector<sequence>& sequences, PatternSorter method);

/**
 Extract patterns and sort them (approximately) to minimise the tour length.
 */
PatternVector tspPatternSort(const vector<sequence>& sequences, long& length);


/**
 Compute the differences between consecutive patterns.
 */
PatternDiffs computePatternDifferences(const PatternVector& patterns);


}
#endif
