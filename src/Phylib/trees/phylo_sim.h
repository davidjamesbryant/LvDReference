/**
 * phylo_sim.h
 *
 * Utilities for generating random phylogenies.
 *
 */

#ifndef PHYLO_SIM_H_
#define PHYLO_SIM_H_

#include"../global/stdIncludes.h"
#include"../utilities/phylibException.h"
#include "phylo.h"
#include "phylo_util.h"
#include <random>
#include"../utilities/random.h"

namespace Phylib {

template<typename T>  phylo<T>  sim_uniform_phylo(int ntax) {
    if (ntax<=0)
        throw PhylibException("Simulating tree with 0 or negative taxa!");
    
    phylo<T> tree;
    
    //Initial tree has one leaf.
    auto p = tree.insert_child(tree.header());
    p->id = 0;
    
    
    int nnodes = 1;
    for (int i = 1; i<ntax; i++) {
        typename phylo<T>::iterator p = tree.root(); //position to insert next leaf
        int pos = std::rand()%nnodes; //Number of steps from root for this position
        for (int j=0; j<pos; j++) {
            p = p.next_pre();
        }
        //Insert node above p.
        tree.insert_sibling(p); //add sibling
        tree.graft_child(p.right(), tree, p); //move node to below sibling
        //Insert new leaf as sibling of p.
        tree.insert_sibling(p);
        p.right()->id = i;
        nnodes+=2;
    }
    
    for(auto p=tree.root().next_pre();!p.null();p=p.next_pre())
        p->length = random_exp(1.0);

    return tree;
}




template<typename T> phylo<T> sim_caterpillar_phylo(int ntax) {
    
    typedef typename phylo<T>::iterator ITERATOR;
    phylo<T> tree;
    
    vector<int> leaves(ntax);
    for(int i=0;i<ntax;i++)
        leaves[i] = i;
    shuffle(leaves);
    //std::random_shuffle(leaves.begin(),leaves.end());
    //Initial tree has one leaf.
    tree.insert_child(tree.header());
    tree.root()->id = leaves[0];
    ITERATOR bottom = tree.root(); //Bottom root in the tree.
    
    for(int i=1;i<ntax;i++) {
        //Insert node below bottom
        ITERATOR newnode1 = tree.insert_child(bottom);
        ITERATOR newnode2 = tree.insert_sibling(newnode1);
        newnode1->id = bottom->id;
        bottom->id = -1;
        newnode2->id = leaves[i];
        bottom = newnode2;
    }
    
    for(auto p=tree.root().next_pre();!p.null();p=p.next_pre())
        p->length = random_exp(1.0);
    
    return tree;
}




/* Constructs a tree with the coalescent distribution, with given birth rate lambda.
 This version works in O(n) time. It requires T to inherit basic_newick.
 */

template<typename T> phylo<T> coalescent(uint ntax, const double birthRate) {
    phylo<T> tree;
    
    typedef typename phylo<T>::iterator ITERATOR;
    
    tree.clear();
    tree.insert_child(tree.header()); //Inserts the root node.
    
    //Add ntax leaves branching off the root.
    //We store a list of heights of the current maximum nodes.
    list<double> heights;
    double height = 0.0;
    
    for(uint i=0;i<ntax;i++) {
        ITERATOR p = tree.insert_child(tree.root());
        p->length = 0;
        p->id = ntax - i-1; //So that ids go from 0 to ntax - 1
        heights.push_front(0.0);
    }
    
    
    for (int r = ntax;r>=2;r--) {
        
        height += random_exp(birthRate/r);
        
        //Choose two to amalgamate. First choose a pair of numbers from 0...r-1 such that i \neq j.
        //Then choose the iterators pointing to the corresponding nodes and height entries.
        int i = random_num(r);
        int j = random_num(r-1);
        if (j==i)
            j++;
        
        ITERATOR x,y;
        list<double>::iterator bx,by;
        
        
        x = tree.root().left();
        bx = heights.begin();
        int k=0;
        while(k<min(i,j)) {
            x = x.right();
            bx++;
            k++;
        }
        
        k++;
        y = x.right();
        by = bx; by++;
        while (k<max(i,j)) {
            y = y.right();
            by++;
            k++;
        }
        
        //The length of these branches is the  height minus their current height.
        x->length = height - (*bx);
        y->length = height - (*by);
        
        //Add a new node z which will be the parent of x and y.
        if (r>2) {
            ITERATOR z = tree.insert_sibling(y);
            z->length = 0.0;
            *by = height;
            heights.erase(bx);
            tree.graft_child(z,tree,y);
            tree.graft_child(z,tree,x);
        }
    }
    heights.clear();
    return tree;
}

inline int randomSplit(int n) {
    vector<double> p(n+1);
    p[0] = p[n] = 0.0;
    for(int i=1;i<n;i++)
        p[i] = 1.0/(i*(n-i));
    return random_discrete(p);
}

inline double harmonic(int n) {
    double hn = 0.0;
    for(int i=1;i<=n;i++)
        hn+=1.0/i;
    return hn;
}

template<typename T> phylo<T> sim_beta_critical_recurse(uint ntax) {
    phylo<T> tree;
    auto root = tree.insert_child(tree.header());
    if (ntax == 1)
        return tree;
    int m = randomSplit(ntax);
    
    phylo<T> tree1 = sim_beta_critical_recurse<T>(m);
    phylo<T> tree2 = sim_beta_critical_recurse<T>(ntax-m);
    
    auto u1 = tree.graft_child(root,tree1);
    auto u2 = tree.graft_sibling(u1,tree2);
    u1 -> length = random_exp(1.0/harmonic(m));
    u2 -> length = random_exp(1.0/harmonic(ntax-m));
    
    return tree;
}


template<typename T>  phylo<T> sim_beta_critical_tree(uint ntax, bool make_ultrametric = false) {
    phylo<T> tree = sim_beta_critical_recurse<T>(ntax);
    
    //Shuffle leaf ids.
    vector<int> leaves(ntax);
    for(int i=0;i<ntax;i++)
        leaves[i] = i;
    shuffle(leaves);
    
    int i=0;
    //Add taxon labels
    for(auto p = tree.leftmost_leaf();!p.null();p=p.next_post()) {
        if (p.leaf())
            p->id = leaves[i++];
    }
    if (make_ultrametric) {
        double height = maxLeafPath<T>(tree.root());
        makeUltrametric<T>(tree.root(),height);
    }
    
    return tree;
}





}

#endif /*PHYLO_SIM_H_*/




