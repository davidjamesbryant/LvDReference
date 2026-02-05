/**
 * @file phylo_util.h
 * @brief Header file various utility routines on phylogenies
 * @author David Bryant
 * @version 0.1
 */



#ifndef PHYLO_UTIL_INCLUDE
#define PHYLO_UTIL_INCLUDE

#include"../global/stdIncludes.h"
#include"../utilities/phylibException.h"
#include"../utilities/random.h"
#include"phylo.h"


namespace Phylib {
    template <typename T> int count_rooted_cherries(const phylo<T>& tree) {
        typedef typename phylo<T>::const_iterator ITERATOR;
        int count = 0;
        for(ITERATOR p = tree.root();!p.null();p=p.next_pre()) {
            if (!p.leaf()&& p.left().leaf()) {
                ITERATOR q = p.left();
                if (!q.right().null() && q.right().leaf() && q.right().right().null())
                    count++;
            }
        }
        return count;
    }
    
    template <typename T> int count_leaves(const phylo<T>& tree) {
        typedef typename phylo<T>::const_iterator ITERATOR;
        int count = 0;
        for(ITERATOR p = tree.root();!p.null();p=p.next_pre()) {
            if (p.leaf())
                count++;
        }
        return count;
    }


    template <typename T> int count_unrooted_cherries(const phylo<T>& tree) {
        typedef typename phylo<T>::const_iterator ITERATOR;
        int count = count_rooted_cherries(tree);
        //One more cherry can occur if the tree is rooted along a pendant branch in a cheery.
        ITERATOR root = tree.root();
        
        int nc = 0, nlc = 0; //Number children and number leaf children of root
        for(ITERATOR p = root.left();!p.null();p=p.right()) {
            nc++;
            if (p.leaf())
                nlc++;
        }
        
        if (nc==2 && nlc==1) {
            //Check to see if nonleaf child also has exactly one leaf child.
            
            ITERATOR child = root.left();
            if (child.leaf())
                child=child.right();
            
            nlc = 0;
            for(ITERATOR p = child.left();!p.null();p=p.right()) {
                if (p.leaf())
                    nlc++;
            }
            if (nlc==1)
                count++;
        }
        return count;
        
    }
    
    //TODO: Make 'from' const.
    /**
     Copies the topology and the id fields only.
     **/
    template<typename A, typename B> void copy_topology_and_taxa(const phylo<A>& from, phylo<B>& to) {
        copy_topology(from,to);
        typename phylo< B>::iterator q = to.root();
        for(typename phylo< A>::const_iterator p=from.root();!p.null();p=p.next_pre()) {
            q->id = p->id;
            q = q.next_pre();
        }
    }
    
    
    /**
     Copies the topology and the data fields. It assumes that type B has an assignment operator which can copy from type a.
     **/
    template<typename A, typename B> void copy(phylo<A>& from, phylo<B>& to) {
        copy_topology(from,to);
        typename phylo< B>::iterator q = to.root();
        for(typename phylo< A>::const_iterator p=from.root();!p.null();p=p.next_pre()) {
            *q = *p;
            q = q.next_pre();
        }
    }
    
    
    
    /**
     * ReRoot the topology.
     * Makes the node newRoot the root of the new tree. All iterators will remain valid, so a pointer to the old root
     * will remain pointing to that node (though it will no longer be the root). This is useful for eliminating
     * degree two vertices that may be created.
     */
    template<typename T> void reRoot(phylo<T>& tree, typename phylo<T>::iterator newRoot) {
        typedef typename phylo<T>::iterator ITERATOR;
        
        if (newRoot.par()!=tree.header()) {
            reRoot(tree,newRoot.par()); //Recurse. After this step, the root is immediately above this node.
            ITERATOR oldRoot = newRoot.par();
            oldRoot->length = newRoot->length;
            newRoot->length = 0.0;
            tree.graft_child(tree.header(),tree,newRoot); //Detach the tree rooted at the new root and make it a new root
            tree.graft_child(newRoot,tree,oldRoot); //Move the old tree to make it a subtree of the new root
        }
        
    }
    
    
    /**
     Reverses the order of the children of the given node.
     **/
    template<typename T> void reverseChildren(phylo<T>& tree, typename phylo<T>::iterator node) {
        
        if (node.leaf())
            return;
        
        typename phylo<T>::iterator me=node.left();
        while(!me.right().null())
            tree.graft_child(node,tree,me.right());
    }
    
    
    
    
    
    /**
     computeHeights
     Computes the heights for each node.
     The height is defined to be the length of the path to the left most leaf descendent.
     
     Use checkClock to test if the tree is actually clock-like (ultrametric).
     
     Assumes that the note (type T) has a double called height, which is filled with the height values.
     **/
    template<typename T> void computeHeights(phylo<T>& tree) {
        typedef typename phylo<T>::iterator ITERATOR;
        
        for(ITERATOR p = tree.leftmost_leaf();!p.null();p=p.next_post()) {
            if (p.leaf())
                p->height = 0.0;
            else
                p->height = (p.left()->height + p.left()->length);
        }
    }
    
    /**
    Rescale the lengths of the branches, given a scaling value lambda
     **/
    template<typename T> void scale(phylo<T>& tree, double lambda) {
        
        for(auto p = tree.leftmost_leaf();!p.null();p=p.next_post())
            p->length *= lambda;
    }
    
    /**
    Returns the sum of branch lengths in a tree
     **/
    template<typename T> double computeLength(phylo<T>& tree) {
        double total = 0.0;
        for(auto p = tree.leftmost_leaf();!p.null();p=p.next_post())
            total+= p->length;
        return total;
    }
    
    /**
     Computes the lengths of the branches, given the values for the heights.
     **/
    template<typename T> void computeLengthsFromHeights(phylo<T>& tree) {
        typedef typename phylo<T>::iterator ITERATOR;
        
        for(ITERATOR p = tree.leftmost_leaf();!p.null();p=p.next_post()) {
            if (p.leaf())
                p->length = p.par()->height;
            else
                p->length = (p.par()->height - p->height);
        }
    }
    
    /**
     CheckClock
     
     Computes the heights for each node, checking that the tree is clock-like.
     The height is defined to be the length of the path to the left most leaf descendent.
     
     Assumes that the datatype has a double field called height, which is filled with these values.
     
     EPSILON is the precision to which we evalute the heights.
     **/
    template<typename T> bool checkClock(phylo<T>& tree, double EPSILON = 1e-5) {
        typedef typename phylo<T>::iterator ITERATOR;
        
        computeHeights(tree);
        
        for(ITERATOR p = tree.root();!p.null(); p=p.next_pre()) {
            if (!p.leaf()) {
                for(ITERATOR q = p.left().right();!q.null();q=q.right())
                    if (std::abs(p->height - (q->height + q->length))>EPSILON)
                        return false;
            }
        }
        return true;
    }


 
 
     /**
     * makeUltrametric
     * 
     * Make the tree ultrametric 
     */
    template<typename T> void makeUltrametric(typename phylo<T>::iterator p, double max, double dist2root=0.0) {
		if (p.leaf())
            p->length +=max-dist2root;
		else{
			 for(auto q = p.left();!q.null();q=q.right()) 
			 	makeUltrametric<T>(q,max,dist2root + q->length);
		}    
    }
 
 
      /**
     * computeTimes
     * 
     * Compute computeTimes and store them in height
     */
    template<typename T> void computeTimes(typename phylo<T>::iterator p,double time) {        

        if(!p.root()){
        	p->height = p->length + time;
        	time += p->length ;
        }	
        else{
        	p->height =0;
        } 
        for(auto q = p.left();!q.null();q=q.right()) {
           computeTimes<T>(q, time) ;
           
        }       
    }
       
     /**
     * maxLeafPath
     * 
     * Compute the longest path between the root and any leaf.
     */
    template<typename T> double maxLeafPath(typename phylo<T>::iterator p) {        

        double max_height = 0;
        for(auto q = p.left();!q.null();q=q.right()) {
            double h = maxLeafPath<T>(q) ;
            if (h > max_height)
               max_height = h;
        }
        if(!p.root())
        	max_height += p->length;

        return max_height;
                
    }
    
      
    /**
     * nodeHeight
     * 
     * Compute the height of the tree (in number of edges from the node to the furthest descendent leaf).
     * Does not require a 'height' field in the node
     */
    template<typename T> int nodeHeight(typename phylo<T>::iterator p) {        
        if (p.leaf())
            return 0;
        else {
            int max_height = 0;
            for(auto q = p.left();!q.null();q=q.right()) {
                int h = nodeHeight<T>(q) + 1;
                if (h > max_height)
                    max_height = h;
            }
            return max_height;
        }        
    }

/**
 * root_to_tip
 *
 * Compute the height of the tree (in total path length to furthest descendent leaf).
 * Does not require a 'height' field in the node
 */
template<typename T> double root_to_tip(typename phylo<T>::iterator p) {
    if (p.leaf())
        return 0;
    else {
        double max_height = 0;
        for(auto q = p.left();!q.null();q=q.right()) {
            double h = root_to_tip<T>(q) + q->length;
            if (h > max_height)
                max_height = h;
        }
        return max_height;
    }
}


    /**
     * summedNodeHeight
     *
     * Compute the  number of edges from the root to a leaf summed over all leaves
     */
    template<typename T> int summedNodeHeight(typename phylo<T>::iterator p) {
        int nleaves;
        return summedNodeHeight<T>(p,nleaves);
    }
    
    template<typename T> int summedNodeHeight(typename phylo<T>::iterator p, int& nleaves) {
        if (p.leaf()) {
            nleaves = 1;
            return 0;
        }
        else {
            int height_sum = 0;
            nleaves = 0;
            for(auto q = p.left();!q.null();q=q.right()) {
                int subtreeLeaves;
                auto childSum = summedNodeHeight<T>(q,subtreeLeaves);
                height_sum += childSum + subtreeLeaves;
                nleaves += subtreeLeaves;
            }
            return height_sum;
        }
    }

    /**
          checkIfBinary
     
            Checks that every non-leaf node has exactly two children.
     */
    template<typename T> bool checkIfBinary(phylo<T>& tree) {
        typedef typename phylo<T>::iterator ITERATOR;
        bool isBinary = true;
        
        for (ITERATOR p = tree.root();!p.null();p=p.next_pre()) {
            if (!p.left().null()) {
                if (p.left().right().null())
                    isBinary = false;
                else if (!p.left().right().right().null())
                    isBinary = false;
            }
        }
        return isBinary;
    }



    /**
      * resolveTree
      *
      * Any nodes with more than three children are resolved arbitrarily (actually, using a caterpillar *tree...
    *
      */
    template<typename T> int resolveTree(phylo<T>& tree) {
        typedef typename phylo<T>::iterator ITERATOR;
        
        ITERATOR p = tree.root();
        int maxNumChildren = 0;
        
        while(!p.null()) {
            int numc = 0;
            //Count the number of children
            for (ITERATOR q = p.left();!q.null();q=q.right())
                numc++;
            maxNumChildren = max(maxNumChildren,numc);
            
            if (numc == 1) { //Surpress degree one node (ignores branch length)
                ITERATOR q = p.left();
                tree.contract(p);
                p = q;
                continue;
            }
            while (numc > 2) {
                //Move first two children to be children of a new vertex
                ITERATOR q = tree.insert_child(p);
                tree.graft_child(q, tree, p.left().right());
                tree.graft_sibling(q.left(), tree, p.left().right());
                numc--;
            }
            p = p.next_pre();
        }
        return maxNumChildren;
    }

    /* Constructs the upgma tree given a distance matrix dist. The id values of the leaves are
     set to 0...ntax-1 and and internal nodes get id 2. The branch lengths are given by the algorithm.
     Hence it is assumed that T contains the fields of basic_newick as well as a height field
     */
    //test
    template<typename T> void upgma(phylo<T>& tree, const vector<vector<double> >& dist) {
        uint ntax = static_cast<uint>(dist.size());
        tree.clear();
        
        //Initialise a copy of the distance matrix.
        vector<vector<double> > D(dist.size());
        for(uint i=0;i<ntax;i++) {
            D[i].resize(ntax);
            std::copy(dist[i].begin(),dist[i].end(),D[i].begin());
        }
        
        typedef typename phylo<T>::iterator ITERATOR;
        
        //Create the ntax subtrees containing the leaves.
        vector< phylo<T> > subtrees(ntax);
        for(uint i=0;i<ntax;i++) {
            ITERATOR p = subtrees[i].insert_child(subtrees[i].header());
            p->id = i;
            p->height = 0.0;
        }
        vector<uint> n(ntax); //Number of taxa in each cluster
        fill(n.begin(),n.end(),1);
        
        //Begin the loops
        for(uint r=ntax;r>1;r--) {
            //Choose two to amalgamate.
            double min_d;
            int min_i,min_j;
            min_i = -1; //Indicates that we haven't found a pair yet.
            for(uint i=0;i<ntax-1;i++) {
                if (n[i]==0)
                    continue;
                for(uint j=i+1;j<ntax;j++) {
                    if (n[j]==0)
                        continue;
                    double dij = dist[i][j]/((double)n[i]*n[j]);
                    if (min_i < 0 || min_d > dij) {
                        min_i = i; min_j = j; min_d = dij;
                    }
                }
            }
            
            phylo<T> newNode;
            ITERATOR newRoot = newNode.insert_child(newNode.header());
            newRoot->id = -1;
            newNode.graft_child(newRoot,subtrees[min_j],subtrees[min_j].root());
            newNode.graft_child(newRoot,subtrees[min_i],subtrees[min_i].root());
            newRoot->height = std::max(min_d/2.0, newRoot.left()->height);
            newRoot->height = std::max(newRoot->height, newRoot.left().right()->height);
            
            n[min_i] += n[min_j];
            n[min_j] = 0;
            subtrees[min_j].clear();
            subtrees[min_i].clear();
            subtrees[min_i].graft_child(subtrees[min_i].header(),newNode,newNode.root());
        }
        
        //At this point, subtres[0] contains the root of the upgma tree.
        computeLengthsFromHeights(subtrees[0]);
        tree = subtrees[0];
        subtrees[0].clear();
    }
    
       
    
}

#endif
    
