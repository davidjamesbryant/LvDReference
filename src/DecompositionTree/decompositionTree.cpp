/**
 * @File decompositionTree.cpp
 * @brief Routines for constructing and updating decomposition trees. Uses Eigen to speed up array calculations
 * @author David Bryant & Celine Scornavacca
 * @version 0.1
 */

#include "decompositionTree.h"

static const int NULLPTR = -1;



//TODO: When merging we shouldn't need to check that the components have null parents unless there is a bug with the pointer handling. Perhaps run a test for this.
class DecomVertexData: public Phylib::basic_newick {
    
public:
    int bndry0,bndry1,bndry2; //bndry0 - segment with this node as a bud. //bndry1, bndry2 segments which this node as a root.
    int num;
    int height; //Max length, in number of edges, of a path from this node to a descendent.
    
    DecomVertexData(): Phylib::basic_newick() {
        bndry0 = bndry1 = bndry2 = NULLPTR;
        num=-1;
    }
    DecomVertexData(const basic_newick& node): basic_newick(node){
        bndry0 = bndry1 = bndry2 = NULLPTR;
        num=-1;
    }
    
};


class DecomTableData {
    
public:

    DecomTableData() {
        lchild = rchild = par = NULLPTR;
        mergeType = 0;
        root.set_null();
        bud.set_null();
        length = -1;
        id = -1;
    }
    
    int lchild, rchild, par; //S1 and S2 in the paper.
    phylo<DecomVertexData>::iterator root, bud;
    int mergeType;
    
    //Store here any data fields to preserve from the original phylogeny
    Scalar length; //Fields from basic_nexus
    int id;
    
    bool is_clade() const {return bud.null();}
};


static int countTaxa(const phylo<basic_newick>& input) {
    int n=0;
    //typedef typename phylo<basic_newick>::iterator ITERATOR;
    using ITERATOR = phylo<basic_newick>::iterator;

    for(ITERATOR p = input.root();!p.null();p=p.next_pre()) {
        if (p.leaf())
            n++;
    }
    return n;
}

/*******
DEBUGGING CODE*/

static void printDecomPhylo(ostream& os, phylo<DecomVertexData>::iterator subtree, int indent) {
    typedef typename phylo<DecomVertexData>::iterator ITERATOR;

    for(int k=0;k<indent;k++)
        os<<"    ";
    os<<subtree->num<<": bndry0 = "<<subtree->bndry0<<" bndry1 = "<<subtree->bndry1<<" bndry2 = "<<subtree->bndry2<<" id = "<<subtree->id<<" length = "<<subtree->length<<endl;
    
    for(ITERATOR p=subtree.left(); !p.null(); p=p.right())
        printDecomPhylo(os,p,indent+1);
}

static void printDecomPhylo(ostream& os, phylo<DecomVertexData>& tree) {
	printDecomPhylo(os, tree.root(), 0);

}

static void printDecomTable(ostream& os,  vector<DecomTableData>& dtable) {


	os<<"\tlchild\trchild\tpar\tlength\tid\troot\tbud\tmergeType\n";
	for(int i=0;i<(int)dtable.size();i++) {
		os<<i<<"\t"<<dtable[i].lchild<<"\t"<<dtable[i].rchild<<"\t"<<dtable[i].par<<"\t";
		os<<dtable[i].length<<"\t"<<dtable[i].id<<"\t";
		if (dtable[i].root.null())
			os<<"NULL\t";
		else
			os<<(dtable[i].root)->num<<"\t";
		if (dtable[i].bud.null())
			os<<"NULL\t";
		else
			os<<(dtable[i].bud)->num<<"\t";
		os<<dtable[i].mergeType;
		os<<"\n";
	
	}
}


/*
 Utility class which stores a prospective merger. S1 and S2 are components, labelled to correspond to figure in the
 paper, v is the boundary node they share in common
 
 */
class Merger {
public:
    int S1,S2;
    phylo<DecomVertexData>::iterator v;
    int mergeType;
    Merger(int new_S1, int new_S2, phylo<DecomVertexData>::iterator new_v, int newType): S1(new_S1), S2(new_S2), v(new_v), mergeType(newType) {}
    Merger(): S1(0), S2(0), mergeType(-1) {v.set_null();}
    Merger(const Merger& M) {
        S1 = M.S1; S2 = M.S2; v = M.v; mergeType = M.mergeType;
    }
};


/*
 Finds the first element in the list which has a valid merger with this one.
 Valid means that the merger gives a segment AND i>j.
 
 */
static Merger* getMinValidNeighbor(int i, vector<DecomTableData>& dtable) {
    Merger* minMerge = nullptr;
    int min_j=i,j;
    phylo<DecomVertexData>::iterator v;
    
    //First check if component i can be merged with the component above it.
    v = dtable[i].root;
    j = v->bndry0; //Find component j v as a bud.
    if (j!=NULLPTR && j<min_j && dtable[j].par==NULLPTR) {
        if (dtable[i].is_clade() && v->bndry2!=NULLPTR) 
            minMerge = new Merger(j,i,v,3);
        else if(dtable[i].is_clade() && v->bndry2==NULLPTR) 
            minMerge = new Merger(j,i,v,4);
        else if (v->bndry2==NULLPTR)
            minMerge = new Merger(j,i,v,5);
        if (minMerge!=nullptr)
            min_j = j;
    }

    //Now check if component i can be merged with its sibling
    j = v->bndry1;
    if (i==j)
        j = v->bndry2;
    if (j!=NULLPTR && j<min_j && dtable[j].par==NULLPTR) {
        if (dtable[i].is_clade()) {
            if (dtable[j].is_clade())
                minMerge = new Merger(i,j,v,1);
            else
                minMerge = new Merger(i,j,v,2);
            min_j = j;
        } else if (dtable[j].is_clade()) {
            minMerge = new Merger(j,i,v,2);
            min_j = j;
        }
    }

    //Now check if component i can be merged with the component below it.
    v = dtable[i].bud;
    if (!v.null()) {
        if (v->bndry2==NULLPTR) {
            j = v->bndry1;
            if (j!=NULLPTR && j<min_j && dtable[j].par==NULLPTR) {
                if (dtable[j].is_clade())
                    minMerge = new Merger(i,j,v,4);
                else
                    minMerge = new Merger(i,j,v,5);
                min_j = j;
            }
        } else {
            j = v->bndry1;
            if (j!=NULLPTR && j<min_j && dtable[j].par==NULLPTR && dtable[j].is_clade()) {
                minMerge = new Merger(i,j,v,3);
                min_j = j;
            } 
            j = v->bndry2;
            if (j!=NULLPTR && j<min_j && dtable[j].par==NULLPTR && dtable[j].is_clade()) {
                minMerge = new Merger(i,j,v,3);
                min_j = j;
            }
        }
    }

    return minMerge;
}


/*
 Get a list of valid mergtes which can be made with this one.
 Valid means that the merger gives a segment AND i>j.
 
 Attempts to avoid a type V merger, unless that is the only option.
 */
static Merger* getGreedyValidNeighbor(int i, vector<DecomTableData>& dtable) {
    
    int j;
    list<Merger> validMergers;
    
    phylo<DecomVertexData>::iterator v;
    
    //First check if component i can be merged with the component above it.
    v = dtable[i].root;
    j = v->bndry0; //Find component j v as a bud.
    if (j!=NULLPTR && dtable[j].par==NULLPTR) {
        if (dtable[i].is_clade() && v->bndry2!=NULLPTR)
            validMergers.push_back(Merger(j,i,v,3));
        else if(dtable[i].is_clade() && v->bndry2==NULLPTR)
            validMergers.push_back(Merger(j,i,v,4));
        else if (v->bndry2==NULLPTR)
            validMergers.push_back(Merger(j,i,v,5));
    }

    //Now check if component i can be merged with its sibling
    j = v->bndry1;
    if (i==j)
        j = v->bndry2;
    if (j!=NULLPTR  && dtable[j].par==NULLPTR) {
        if (dtable[i].is_clade()) {
            if (dtable[j].is_clade())
                validMergers.push_back(Merger(i,j,v,1));
            else
                validMergers.push_back(Merger(i,j,v,2));
        } else if (dtable[j].is_clade()) {
            validMergers.push_back(Merger(j,i,v,2));
        }
    }

    //Now check if component i can be merged with the component below it.
    v = dtable[i].bud;
    if (!v.null()) {
        if (v->bndry2==NULLPTR) {
            j = v->bndry1;
            if (j!=NULLPTR  && dtable[j].par==NULLPTR) {
                if (dtable[j].is_clade())
                    validMergers.push_back(Merger(i,j,v,4));
                else
                    validMergers.push_back(Merger(i,j,v,5));
            }
        } else {
            j = v->bndry1;
            if (j!=NULLPTR  && dtable[j].par==NULLPTR && dtable[j].is_clade()) {
                validMergers.push_back(Merger(i,j,v,3));
            }
            j = v->bndry2;
            if (j!=NULLPTR  && dtable[j].par==NULLPTR && dtable[j].is_clade()) {
                validMergers.push_back(Merger(i,j,v,3));
            }
        }
    }

    //We now have a list of valid mergers. Identify the non-V type merge with minimum index, otherwise identify the V type with minimum index.
    int min_notV_j = -1;
    int min_V_j = -1;
    Merger minMerge, minVMerge;
    
    for (Merger M : validMergers) {
        int index = min(M.S1,M.S2);
        if (M.mergeType==5) {
            if (min_V_j < 0 || index<min_V_j) {
                min_V_j = index;
                minVMerge = M;
            }
        } else {
            if (min_notV_j < 0 || index<min_notV_j) {
                min_notV_j = index;
                minMerge = M;
            }
        }
    }
    if (min_notV_j >=0)
        return new Merger(minMerge);
    else if (min_V_j)
        return new Merger(minVMerge);
    else
        return nullptr;
}


/*
 Identify what type of merger this is (according to diagram in paper). Note we are assuming that S1 and S2 correspond to the components in that diagram.
 */
static int getMergeType(int S1, int S2,vector<DecomTableData>& dtable) {
    if (dtable[S1].is_clade()) {
        if (dtable[S2].is_clade())
            return 1;
        else
            return 2;
    }
    if (dtable[S1].bud->bndry2 != NULLPTR)
        return 3;
    if (dtable[S2].is_clade())
        return 4;
    return 5;
}

static void buildDecomTable(phylo<DecomVertexData>& tree, int n, vector<DecomTableData>& dtable, bool useGreedy, bool debug = false) {
    typedef typename phylo<DecomVertexData>::iterator ITERATOR;
    
    dtable.resize(4*n-5);
    int i = 0, j = n; //Index for leaves and for internal
    
    //Number the nodes. This is mainly for debugging only.
    for(ITERATOR v = tree.leftmost_leaf();!v.null();v = v.next_post()) {
        if (v.leaf())
            v->num = i++;
        else
            v->num = j++;
    }

    i=0; j=n;
    //Fill out the components corresponding to edges.
    for(ITERATOR v = tree.leftmost_leaf();!v.root();v = v.next_post()) {
        if (v.leaf()) {
            dtable[i].root = v.par();
            if (dtable[i].root -> bndry1 == NULLPTR)
                dtable[i].root -> bndry1 = i;
            else
                dtable[i].root -> bndry2 = i;
            dtable[i].id = v->id;
            //cout<<"**"<<v->id<<endl;
            dtable[i].length = v->length;
            dtable[i].mergeType = 0;
            i++;
        } else {
            dtable[j].root = v.par();
            if (dtable[j].root -> bndry1 == NULLPTR)
                dtable[j].root -> bndry1 = j;
            else
                dtable[j].root -> bndry2 = j;
            dtable[j].bud = v;
            v -> bndry0 = j;
            dtable[j].length = v->length;
            dtable[i].mergeType = 0;
            j++;
        }
    }
    
    if (debug) {
      printDecomPhylo(cout, tree); cout<<endl; printDecomTable(cout,dtable); //DEBUG
    }
    
    /*Now the mergers. We loop through the components, adding any newly formed components to the end of the list. Initially
     we just have the components corresponding to pendant edges, and then internal edges.
     
     Define the 'layer' of a component by letting all edges have layer one, and the layer of a new components formed by a merger
     equal to one more than the maximum of the layers for components it is formed from. Eqauivalently, the layer is the height
     in the decomposition tree.
     
     The way we loop through components means that we consider a maximal set of mergers for each layer, and we don't start going on to the
     next layer until a maximal set of mergers has been attempted in the previous layer. This is important for the log(n) height bound, and can
     no doubt be described way better than it has been here.
    */
    
    int k = j;
    
    for(int i=0;i<4*n-5;i++) {
        Merger* merge;
        if (useGreedy)
            merge = getGreedyValidNeighbor(i,dtable);
        else
            merge = getMinValidNeighbor(i,dtable);
        
        if (merge==nullptr)
            continue;
        
        if (debug) {
          cout<<"Merging S1 = "<<merge->S1<<" S2 = "<<merge->S2<<" to give "<<k; //DEBUG
          cout<<" with merge type "<<merge->mergeType<<endl; //DEBUG
        }
        
        int S1 = merge->S1, S2 = merge->S2, m = merge->mergeType;

        //DEBUG
        if (m!=getMergeType(S1,S2,dtable)) {
            cerr<<"Merge type mismatch "<<m<<" "<<getMergeType(S1,S2,dtable)<<endl;
            exit(1);
        }   


        auto u = dtable[S1].root;
        auto v = dtable[S1].bud;
        auto w = dtable[S2].bud;
        dtable[k].mergeType = m;
        dtable[k].lchild = S1;
        dtable[k].rchild = S2;
        dtable[S1].par = k;
        dtable[S2].par = k;
        
        switch(m) {
            case 1:
                dtable[k].root = u;
                u->bndry1 = k;
                u->bndry2 = NULLPTR;
                break;
            case 2:
                dtable[k].root = u; 
                dtable[k].bud = w;
                u->bndry1 = k;
                u->bndry2 = NULLPTR;
                w->bndry0 = k;
                break;
            case 3:
                dtable[k].root = u;
                if (u->bndry1==S1)
                    u->bndry1 = k;
                else
                    u->bndry2 = k;
                dtable[k].bud = v;
                v->bndry0 = k;
                if (v->bndry1==S2)
                    v->bndry1 = v->bndry2;
                v->bndry2 = NULLPTR;
                break;
            case 4:
                dtable[k].root = u;
                if (u->bndry1==S1)
                    u->bndry1 = k;
                else
                    u->bndry2 = k;
                v->bndry0 = v->bndry1 = v->bndry2 = NULLPTR;
                break;
            case 5:
                dtable[k].root = u;
                if (u->bndry1==S1)
                    u->bndry1 = k;
                else
                    u->bndry2 = k;
                v->bndry0 = v->bndry1 = v->bndry2 = NULLPTR;
                dtable[k].bud = w;
                w->bndry0 = k;
                break;
            default:
                cerr<<"Invalid merge type"<<endl;
                exit(1);
        }

       
            
        if (debug&&false) {
          printDecomPhylo(cout, tree); cout<<endl; printDecomTable(cout,dtable); cout<<"=============\n";  //DEBUG 
        }
         
        k++;
    }
}


/*
 Recursively extract the phylo object from the table used to construct the Decomposition tree.
 */
static void  ExtractDecomTree( vector<DecomTableData>& dtable, phylo<DecomNodeData>& dtree, int rootNode) {
    //cout<<"rootNode = "<<rootNode<<endl;
    
    
  typedef typename phylo<DecomNodeData>::iterator ITERATOR;
    
    dtree.clear();
    
    ITERATOR newNode = dtree.insert_child(dtree.header());
    if (dtable[rootNode].lchild != NULLPTR ) {
        phylo<DecomNodeData> subtreeL, subtreeR;
        ExtractDecomTree(dtable,subtreeL,dtable[rootNode].lchild);
        ExtractDecomTree(dtable,subtreeR,dtable[rootNode].rchild);
        ITERATOR child = dtree.graft_child(newNode,subtreeL);
        dtree.graft_sibling(child,subtreeR);
        assert(!newNode.left().null() && !newNode.left().right().null());
    } else {
        newNode->length = dtable[rootNode].length; //Single edge component
        newNode->id = dtable[rootNode].id;
    }
    newNode->mergeType = dtable[rootNode].mergeType;
    newNode->isClade = dtable[rootNode].is_clade();
}


static void  ExtractDecomTree( vector<DecomTableData>& dtable, phylo<DecomNodeData>& dtree) {
    ExtractDecomTree(dtable,dtree,dtable.size()-1);
}

/**
Construct the decomposition tree for a given phylogeny
//TODO Params
**/

void constructDecompTree( phylo<basic_newick>& input, phylo<DecomNodeData>& decomTree, bool useGreedy)  {
    
    
    //Copy the phylogeny to add the extra fields
    Phylib::phylo<DecomVertexData> phylo;
    int n = countTaxa(input);
    copy(input,phylo);
    
    vector<DecomTableData> dtable;
    buildDecomTable(phylo,n,dtable,useGreedy,false);
    
//    cout<<"**** Decomposition table ***\n";
//    printDecomTable(cout,dtable);
//    cout<<endl;


    ExtractDecomTree(dtable,decomTree);
    
    dtable.clear();
    
}


/**
The following routines set up a decomposition tree which is the same as the pruning 
algorithm. This is used just for comparisons.
**/


/*
 Recursively extract the phylo object from the table used to construct the Decomposition tree.
 
 Returns the decomposition tree with top node corresponding to the clade rooted at p 
 together with the edge immediately above p.
 */
static void  ExtractPruningDecomTree( phylo<basic_newick>::iterator p, phylo<DecomNodeData>& dtree) {
    //cout<<"rootNode = "<<rootNode<<endl;
    
    
    typedef typename phylo<DecomNodeData>::iterator ITERATOR;
    dtree.clear();
    ITERATOR newNode = dtree.insert_child(dtree.header());
    newNode->isClade = true;

    if (p.left().null())  {
        //p is a leaf. We just want to return the component containing the edge connecting
        //p to its parent.
        newNode->length = p->length;
        newNode->id = p->id;
        //cerr<<"\t"<<p->id<<endl;
    } else {
    
        //Return  (edge,(subtreeL,subtreeR))
        //This adds three components to the decomposition tree:
        //   newNode   ((subtreeL,subtreeR),edge)
        //   pNode     (subtreeL,subtreeR)
        //   edgeNode  edge

        phylo<DecomNodeData> subtreeL, subtreeR;
        ExtractPruningDecomTree(p.left(), subtreeL);
        ExtractPruningDecomTree(p.left().right(),subtreeR);
        
        ITERATOR edgeNode = dtree.insert_child(newNode);
        edgeNode->length = p->length;
        edgeNode->isClade = false;


        ITERATOR pNode = dtree.insert_sibling(edgeNode);
        ITERATOR nodeL = dtree.graft_child(pNode,subtreeL);
        dtree.graft_sibling(nodeL,subtreeR);
        pNode->mergeType=1; //Merge two clades.
        pNode->isClade = true;

        newNode->mergeType = 4; //S1=edge, S2 = clade with root = bud of S1.
        assert(!newNode.left().right().null());
    }
}

void constructPruningDecomp( phylo<basic_newick>& input, phylo<DecomNodeData>& dtree)  {
    
    typedef typename phylo<DecomNodeData>::iterator ITERATOR;

    dtree.clear();
    ITERATOR rootNode = dtree.insert_child(dtree.header());
    phylo<DecomNodeData> subtreeL, subtreeR;
    ExtractPruningDecomTree(input.root().left(), subtreeL);
    ExtractPruningDecomTree(input.root().left().right(), subtreeR);
    ITERATOR nodeL = dtree.graft_child(rootNode,subtreeL);
    dtree.graft_sibling(nodeL,subtreeR);
    rootNode->mergeType = 1; //Merge of two clades.
    rootNode->isClade = true;
    
    cerr<<"Constructed pruning decom tree"<<endl;
}
      


