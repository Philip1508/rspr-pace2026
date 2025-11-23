/*******************************************************************************
rspr.h

Calculate approximate and exact Subtree Prune and Regraft (rSPR)
distances and the associated maximum agreement forests (MAFs) between pairs
of rooted binary trees.
Supports arbitrary labels. See the
README for more information.

Copyright 2009-2014 Chris Whidden
whidden@cs.dal.ca
http://kiwi.cs.dal.ca/Software/RSPR
April 29, 2014
Version 1.2.2

This file is part of rspr.

rspr is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

rspr is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with rspr.  If not, see <http://www.gnu.org/licenses/>.

*******************************************************************************/

#define RSPR
//#define DEBUG 1
//#define DEBUG_CONTRACTED 1
//#define DEBUG_APPROX 1
//#define DEBUG_CLUSTERS 1
//#define DEBUG_SYNC 1
//#define DEBUG_UNDO 1
//#define DEBUG_DEPTHS 1
//#define DEBUG_CASE_COUNTER 1
//#define MULT_PICK_LARGEST_GROUP 1

#pragma once

#include <cstdio>
#include <cstdlib>
#include <string>
#include <cstring>
#include <iostream>
#include <sstream>
#include <climits>
#include <vector>
#include <map>
#include <set>
#include <list>
#include <algorithm>
#include "../DataStructures/Forest.h"
#include "../DataStructures/ClusterForest.h"
#include "../UndecipheredComponents/LCA.h"
#include "../DataStructures/ClusterInstance.h"
#include "../DataStructures/SiblingPair.h"
#include "../DataStructures/UndoMachine.h"

#include "Utility/rSprUtility.h"
#include "Algorithm/Binary/rSprAlgorithm3Approx.h"
#include "Algorithm/Binary/rSprAlgorithmBB.h"
using namespace std;



enum RELAXATION {STRICT, NEGATIVE_RELAXED, ALL_RELAXED};

const string whitespaces = " \t\f\v\n\r";


// ----------------------------------
// FUNCTION PROTOTYPES
// ----------------------------------
// note: not using undo
int rSPR_worse_3_mult_approx_hlpr(Forest *T1, Forest *T2, list<Node *> *singletons, list<Node *> *sibling_pairs, Forest **F1, Forest **F2, bool save_forests);
int rSPR_worse_3_mult_approx(Forest *T1, Forest *T2);
int rSPR_worse_3_mult_approx(Forest *T1, Forest *T2, bool sync);

int rSPR_branch_and_bound_mult_range(Forest *T1, Forest *T2, int start_k);
int rSPR_branch_and_bound_mult_range(Forest *T1, Forest *T2, int start_k, int end_k);
int rSPR_branch_and_bound_mult(Forest *T1, Forest *T2, int k);
int rSPR_branch_and_bound_mult_hlpr(Forest *T1, Forest *T2, int k, list<Node*> *sibling_groups, list<Node*> *singletons, Node *protected_stack, list<pair<Forest,Forest>> *AFs, int* num_ties);



int rSPR_3_approx_hlpr(Forest *T1, Forest *T2, list<Node *> *singletons,
		list<Node *> *sibling_pairs);
int rSPR_3_approx(Forest *T1, Forest *T2);
int rSPR_worse_3_approx_hlpr(Forest *T1, Forest *T2, list<Node *> *singletons, list<Node *> *sibling_pairs, Forest **F1, Forest **F2, bool save_forests);
int rSPR_worse_3_approx(Forest *T1, Forest *T2);
int rSPR_worse_3_approx(Forest *T1, Forest *T2, bool sync);
int rSPR_worse_3_approx(Node *subtree, Forest *T1, Forest *T2);
int rSPR_worse_3_approx(Node *subtree, Forest *T1, Forest *T2, bool sync);
int rSPR_worse_3_approx_binary_hlpr(Forest *T1, Forest *T2, list<Node *> *singletons, list<Node *> *sibling_pairs, Forest **F1, Forest **F2, bool save_forests);
int rSPR_worse_3_approx_binary(Forest *T1, Forest *T2, bool sync);
int rSPR_worse_3_approx_binary(Forest *T1, Forest *T2);



int rSPR_branch_and_bound(Forest *T1, Forest *T2);
int rSPR_branch_and_bound(Forest *T1, Forest *T2, int k);
int rSPR_branch_and_bound(Forest *T1, Forest *T2, int k,
		map<string, int> *label_map,
		map<int, string> *reverse_label_map);

int rSPR_branch_and_bound_range(Forest *T1, Forest *T2, int end_k);
int rSPR_branch_and_bound_range(Forest *T1, Forest *T2, int start_k,
		int end_k);
int rSPR_branch_and_bound_hlpr(Forest *T1, Forest *T2, int k,
		set<SiblingPair> *sibling_pairs, list<Node *> *singletons, bool cut_b_only,
		list<pair<Forest,Forest> > *AFs, list<Node *> *protected_stack,
		int *num_ties);
int rSPR_branch_and_bound_hlpr(Forest *T1, Forest *T2, int k,
		set<SiblingPair> *sibling_pairs, list<Node *> *singletons, bool cut_b_only,
		list<pair<Forest,Forest> > *AFs, list<Node *> *protected_stack,
		int *num_ties, Node *prev_T1_a, Node *prev_T1_c);



int rSPR_total_approx_distance(Node *T1, vector<Node *> &gene_trees);
int rSPR_total_approx_distance(Node *T1, vector<Node *> &gene_trees,
		int threshold);
int rSPR_total_distance(Node *T1, vector<Node *> &gene_trees);
int rSPR_total_distance(Node *T1, vector<Node *> &gene_trees,
		vector<int> *original_scores);
void rSPR_pairwise_distance(Node *T1, vector<Node *> &gene_trees);
void rSPR_pairwise_distance(Node *T1, vector<Node *> &gene_trees, bool approx);
void rSPR_pairwise_distance(Node *T1, vector<Node *> &gene_trees, int start, int end);
void rSPR_pairwise_distance(Node *T1, vector<Node *> &gene_trees, int start, int end, bool approx);
void rSPR_pairwise_distance(Node *T1, vector<Node *> &gene_trees, int max_spr);
void rSPR_pairwise_distance(Node *T1, vector<Node *> &gene_trees, int max_spr, int start, int end);
void rSPR_pairwise_distance_unrooted(Node *T1, vector<Node *> &gene_trees);
void rSPR_pairwise_distance_unrooted(Node *T1, vector<Node *> &gene_trees, int start, int end);
void rSPR_pairwise_distance_unrooted(Node *T1, vector<Node *> &gene_trees, bool approx);
void rSPR_pairwise_distance_unrooted(Node *T1, vector<Node *> &gene_trees, int start, int end, bool approx);
int rf_total_distance(Node *T1, vector<Node *> &gene_trees);
int rf_total_distance_unrooted(Node *T1, vector<Node *> &gene_trees);
void rf_pairwise_distance(Node *T1, vector<Node *> &gene_trees);
void rf_pairwise_distance(Node *T1, vector<Node *> &gene_trees, int start, int end);
void rf_pairwise_distance_unrooted(Node *T1, vector<Node *> &gene_trees);
void rf_pairwise_distance_unrooted(Node *T1, vector<Node *> &gene_trees, int start, int end);
int rSPR_total_distance_unrooted(Node *T1, vector<Node *> &gene_trees, int threshold);
int rSPR_total_distance_unrooted(Node *T1, vector<Node *> &gene_trees, int threshold, vector<int> *original_scores);



int rSPR_branch_and_bound_simple_clustering(Node *T1, Node *T2, Forest **out_F1, Forest **out_F2);
int rSPR_branch_and_bound_simple_clustering(Node *T1, Node *T2, bool verbose, map<string, int> *label_map, map<int, string> *reverse_label_map);
int rSPR_branch_and_bound_simple_clustering(Node *T1, Node *T2, bool verbose, map<string, int> *label_map, map<int, string> *reverse_label_map, int min_k, int max_k);
int rSPR_branch_and_bound_simple_clustering(Node *T1, Node *T2);
int rSPR_branch_and_bound_simple_clustering(Node *T1, Node *T2, bool verbose);
int rSPR_branch_and_bound_simple_clustering(Node *T1, Node *T2, bool verbose, int min_k, int max_k);
int rSPR_branch_and_bound_simple_clustering(Node *T1, Node *T2, bool verbose, map<string, int> *label_map, map<int, string> *reverse_label_map, int min_k, int max_k, Forest **out_F1, Forest **out_F2);


void reduction_leaf_mult(Forest *T1, Forest *T2);
void reduction_leaf(Forest *T1, Forest *T2);
void reduction_leaf(Forest *T1, Forest *T2, UndoMachine *um);




Node *find_subtree_of_approx_distance(Node *n, Forest *F1, Forest *F2, int target_size);
Node *find_best_root(Node *T1, Node *T2);
double find_best_root_acc(Node *T1, Node *T2);
void find_best_root_hlpr(Node *T2, int pre_separator, int group_1_total,
		int group_2_total, Node **best_root, double *best_root_b_acc);
void find_best_root_hlpr(Node *n, int pre_separator, int group_1_total,
		int group_2_total, Node **best_root, double *best_root_b_acc,
		int *p_group_1_descendants, int *p_group_2_descendants, int *num_ties);
int rf_distance(Node *T1, Node *T2);
int count_differing_bipartitions(Node *n);
bool contains_bipartition(Node *n, int pre_start, int pre_end,
		int group_1_total, int group_2_total, int *p_group_1_descendants,
		int *p_group_2_descendants);
void modify_bipartition_support(Node *T1, Node *T2, enum RELAXATION relaxed);
void modify_bipartition_support(Node *n, Forest *F1, Forest *F2,
		Node *T1, Node *T2, vector<int> *F1_descendant_counts, enum RELAXATION);
void modify_bipartition_support(Forest *F1, Forest *F2, Node *n1);
bool is_nonbranching(Forest *T1, Forest *T2, Node *T1_a, Node *T1_c, Node *T2_a, Node *T2_c);
bool outgroup_root(Node *T, set<string, StringCompare> outgroup);
bool outgroup_root(Node *n, vector<int> &num_in, vector<int> &num_out);
bool outgroup_reroot(Node *n, vector<int> &num_in, vector<int> &num_out);
void count_in_out(Node *n, vector<int> &num_in, vector<int> &num_out,
		set<string, StringCompare> &outgroup);
void randomize_tree_with_spr(Node* T1, Node* T2, int count);
/*Joel's part*/


int rSPR_branch_and_bound_simple_clustering(Forest *T1, Forest *T2, bool verbose, map<string, int> *label_map, map<int, string> *reverse_label_map);
int rSPR_branch_and_bound_simple_clustering(Forest *T1, Forest *T2);
int rSPR_branch_and_bound_simple_clustering(Forest *T1, Forest *T2, bool verbose);


int rSPR_total_distance(Forest *T1, vector<Node *> &gene_trees);

// ----------------------------------
// ENF OF FUNCTION PROTOTYPES
// ----------------------------------



bool BB = false;
bool APPROX_CHECK_COMPONENT = false;
bool APPROX_REVERSE_CUT_ONE_B = false;
bool APPROX_REVERSE_CUT_ONE_B_2 = false;
bool APPROX_CUT_ONE_B = false;
bool APPROX_CUT_TWO_B = false;
bool APPROX_CUT_TWO_B_ROOT = false;
bool APPROX_EDGE_PROTECTION = false;
bool CUT_ONE_B = false;
bool REVERSE_CUT_ONE_B = false;
bool REVERSE_CUT_ONE_B_2 = false;
bool REVERSE_CUT_ONE_B_3 = false;
bool CUT_TWO_B = false;
bool CUT_TWO_B_ROOT = false;
bool CUT_ALL_B = false;
bool CUT_AC_SEPARATE_COMPONENTS = false;
bool CUT_ONE_AB = false;
bool CLUSTER_REDUCTION = false;
bool PREFER_RHO = false;
bool MAIN_CALL = true;
bool MEMOIZE = false;
bool MULTIFURCATING = false;
bool MULT_4_BRANCH = false;
bool USE_CASE_7 = true;
bool ALL_MAFS = false;
int NUM_CLUSTERS = 0;
int MAX_CLUSTERS = -1;
bool UNROOTED_MIN_APPROX = false;
bool VERBOSE = false;
bool CLAMP = false;
int MAX_SPR = 1000;
int CLUSTER_MAX_SPR = MAX_SPR;
int MIN_SPR = 0;
bool FIND_RATE = false;
bool EDGE_PROTECTION = false;
bool EDGE_PROTECTION_TWO_B = false;
bool ABORT_AT_FIRST_SOLUTION = false;
bool PREORDER_SIBLING_PAIRS = false;
bool DEEPEST_ORDER = false;
bool DEEPEST_PROTECTED_ORDER = false;
bool NEAR_PREORDER_SIBLING_PAIRS = false;
bool LEAF_REDUCTION = false;
bool LEAF_REDUCTION2 = false;
bool SPLIT_APPROX = false;
bool IN_SPLIT_APPROX = false;
int SPLIT_APPROX_THRESHOLD = 25;
float INITIAL_TREE_FRACTION = 0.4;
bool COUNT_LOSSES = false;
bool CUT_LOST = false;
bool CHECK_MERGE_DEPTH = false;
bool check_all_pairs = true;
bool PREFER_NONBRANCHING = false;
int CLUSTER_TUNE = -1;
int SIMPLE_UNROOTED_LEAF = 0;
bool SHOW_CLUSTERS = false;



class ProblemSolution {
public:
string T1;
string T2;
int k;

ProblemSolution(Forest *t1, Forest *t2, int new_k) {
	T1 = t1->str();
	T2 = t2->str();
	k = new_k;
}
	};

	map<string, ProblemSolution> memoized_clusters = map<string, ProblemSolution>();

/*******************************************************************************
	RSPR WORSE_3_MULT_APPROX
*******************************************************************************/

/* rSPR_worse_3_mult_approx
 * Calculate an approximate maximum agreement forest and SPR distance for two multifurcating trees
 * RETURN At most 3 times the rSPR distance
 * NOTE: destructive. The computed forests replace T1 and T2.
 * T1 and T2 can be  multifurcating forests.
 */
int rSPR_worse_3_mult_approx(Forest *T1, Forest *T2) {
	return rSPR_worse_3_mult_approx(T1, T2, true);
}

int rSPR_worse_3_mult_approx(Forest *T1, Forest *T2, bool sync) {
	// match up nodes of T1 and T2
	if (sync) {
if (!sync_twins(T1, T2))
	return 0;
	}

	if (LEAF_REDUCTION) {
	  reduction_leaf_mult(T1, T2);
	}
//	cout << "T1: "; T1->print_components();
//	cout << "T2: "; T2->print_components();
	// find sibling pairs of T1
	list<Node *> *sibling_groups = T1->find_sibling_groups();


	/*
	list<Node *>::iterator c;
	list<Node *>::iterator ch;
	
	for (c = sibling_groups->begin(); c != sibling_groups->end(); c++) {
	  cout << "Sibling group parent: " << (*c)->str() << endl;
		list<Node *> children = (*c)->get_children();
		for (ch = children.begin(); ch != children.end(); ch++) {
		  cout << "\tChild: " << (*ch)->str() <<endl;
		}
		list<list<Node *>> *identical_sibling_groups = (*c)->find_identical_sibling_groups();
		if (identical_sibling_groups->size() > 0) {
		  cout << "Printing identical sibling group:" << endl;
		  for (auto g = identical_sibling_groups->begin(); g != identical_sibling_groups->end(); g++) {
			cout << "Sibling identical group:" << endl;
			for (auto ge = g->begin(); ge != g->end(); ge++) {
			  cout << "\tchild:" << (*ge)->get_name() << endl;
			}			  
		  }
		}
		else {		  
		  cout << "No identical sibling groups" <<endl;
		}
	}	
	//return 0;
	*/
	
	
	// find singletons of T2
	list<Node *> singletons = T2->find_singletons();
	//list<pair<Forest,Forest> > AFs = list<pair<Forest,Forest> >();

	Forest *F1;
	Forest *F2;

	T2->max_preorder = T2->components[0]->get_max_preorder_number(0);//preorder_number(0);
	int ans = rSPR_worse_3_mult_approx_hlpr(T1, T2, &singletons, sibling_groups, &F1, &F2, true);

	F1->swap(T1);
	F2->swap(T2);
	sync_twins(T1,T2);


	delete sibling_groups;
	delete F1;
	delete F2;
	return ans;
}


// rSPR_worse_3_mult_approx recursive helper function
int rSPR_worse_3_mult_approx_hlpr(Forest *T1, Forest *T2, list<Node *> *singletons, list<Node *> *sibling_groups, Forest **F1, Forest **F2, bool save_forests) {

  int num_cut = 0;
  Node* previous_group = NULL;
  while(!singletons->empty() || !sibling_groups->empty()) {

    // Case 1 - Remove singletons
    while(!singletons->empty()) {

      Node *T2_a = singletons->back();
      singletons->pop_back();
      #ifdef DEBUG_APPROX
      cout << "Handling singleton: " << T2_a->str() << endl;
      #endif
      // find twin in T1
      Node *T1_a = T2_a->get_twin();
      // if this is in the first component of T_2 then
      // it is not really a singleton.
      // TODO: problem when we cluster and have a singleton as the
      //		first comp of T2
      //    NEED TO MODIFY CUTTING?
      // 		HERE AND IN BB?
      Node *T1_a_p = T1_a->parent();
      if (T1_a_p == NULL)
	continue;


      if (T2_a == T2->get_component(0)){
	if (!T1->contains_rho()) {
	  T1->add_rho();
	  T2->add_rho();
	  num_cut++;
	}
      }

	//continue;

      bool is_sibling_group = T1_a_p->is_sibling_group();
      // cut the edge above T1_a
      T1_a->cut_parent();
      if (!T1_a->is_leaf()) {
	T1_a_p->decrement_non_leaf_children();//although would this ever be a non leaf?
      }
      T1->add_component(T1_a);

      //only contract if one node
      if (T1_a_p->get_children().size() == 1) {
	if (is_sibling_group) {
	  sibling_groups->remove(T1_a_p);
          #ifdef DEBUG_APPROX
	  cout << "Removed ";
	  for (list<Node*>::iterator i = T1_a_p->get_children().begin(); i != T1_a_p->get_children().end(); i++) {
	    cout << (*i)->str();
	  }
	  cout  << " from sibling groups" << endl;
	  #endif
	}
	Node *possible_previous_sibling = T1_a_p->get_children().front();
	bool was_sibling_group = possible_previous_sibling->is_sibling_group();
	Node *node = T1_a_p->contract(true);
	if (node != NULL) {
	  node->recalculate_non_leaf_children(); //can we tell what this would be instead of recalculating?

	  if (node->is_sibling_group()) {
	    if (was_sibling_group) {
	      list<Node*>::iterator i = find(sibling_groups->begin(), sibling_groups->end(), possible_previous_sibling);
	      *i = node;
	    }
	    else{
	      sibling_groups->push_front(node);
	    }
	  }
	}
      }
    }//!singletons->empty()


    if(!sibling_groups->empty()) {
      //Get the first group that has identical sibling groups, otherwise default to the group on the back
      list<Node*>::iterator i = sibling_groups->end();
      i--;
      Node *T1_sibling_group = sibling_groups->back();
      list<list<Node*>> identical_sibling_groups = list<list<Node*>>();
      T1_sibling_group->find_identical_sibling_groups(&identical_sibling_groups);
      for (; i != sibling_groups->begin(); i-- ){
	(*i)->find_identical_sibling_groups(&identical_sibling_groups);
	if (identical_sibling_groups.size() > 0) {
	  T1_sibling_group = (*i);
	  break;
	}
      }

      #ifdef DEBUG_APPROX
      cout << "F2: ";
      T2->print_components();
      cout << endl;
      cout << "F1: ";
      T1->print_components();
      cout << endl;
      #endif

      /*
	 Case where a subset of the group have the same parent both in T1 and T2
	 Step 5 in paper
      */
      // Case 2 - Contract identical sibling pair
      if (identical_sibling_groups.size() > 0) {
	list<list<Node *>>::iterator i;
	for (i = identical_sibling_groups.begin(); i != identical_sibling_groups.end(); i++) {
	  list<Node *> T2_group = (*i);
	  Node *T2_p = T2_group.front()->parent();
	  #ifdef DEBUG_APPROX
	  cout << "Contracting T1... " << endl;
	  #endif
	  Node *T1_group_new = T1_sibling_group->contract_twin_group(&T2_group);
	  #ifdef DEBUG_APPROX
	  cout << "Contracting T2... " << endl;
	  #endif
	  Node *T2_group_new = T2_p->contract_sibling_group(&T2_group);

	  T1_group_new->set_twin(T2_group_new);
	  T2_group_new->set_twin(T1_group_new);

	  // check if T2_p is a singleton after the contraction
	  if (T2_p->is_singleton() && T2_p != T2->get_component(0)) {
	    //(T2_p->is_singleton() && T1_sibling_group != T1->get_component(0) && T2_p != T2->get_component(0)) {
	    singletons->push_front(T2_p);
	  }
	  if (T1_sibling_group->parent() != NULL) {
	    //Check if the contraction made a new sibling group
	    T1_sibling_group->parent()->recalculate_non_leaf_children();
	    if (T1_sibling_group->parent()->is_sibling_group()) {
	      sibling_groups->push_front(T1_sibling_group->parent());
	    }
	  }
	  if (!T1_sibling_group->is_sibling_group()) {
	    sibling_groups->remove(T1_sibling_group);
	  }
	}
      }

      /*
	4 branching case
	Step 6-8 in paper
	Cut above a1, b1, a2, b2
	Part 1: Get the LCA, this is the numbering to the top part
	Part 2: Get subset of sibling group that is descendant of this LCA
	Part 3: Sort them based on depth
	Part 4: Since this is the approximation, we cut deepest 2
      */

      // Case 3
      else {
	#ifdef DEBUG_APPROX
	cout << "Sibling group to be cutting: " << T1_sibling_group->str_subtree() << endl;
	#endif
	//vector of ints describing how many of the siblings are in its descendants, indexed by preorder number
	vector<int> descendant_count = T1_sibling_group->find_pseudo_lca_descendant_count(T2->max_preorder + 1);
	Node* arbitrary_lca = T1_sibling_group->find_arbitrary_lca(T2->components, descendant_count);
	vector<Node *> deepest_siblings;

	//All siblings are in different components, ie no path between them
	//Get depth of siblings from root of each component
	if (arbitrary_lca == NULL) {
	  vector<vector<Node *>> siblings_by_depth = vector<vector<Node *>>(10);
	  for (int i = 0; i != T2->components.size(); i++) {
	    T2->components[i]->get_deepest_siblings(descendant_count, siblings_by_depth);
	  }
	  deepest_siblings = contract_deepest_siblings(siblings_by_depth);
	}
	//Otherwise they share an LCA
	else {
	  vector<vector<Node *>> siblings_by_depth = arbitrary_lca->get_deepest_siblings(descendant_count);
	  deepest_siblings = contract_deepest_siblings(siblings_by_depth);
	}

	// Should assert here
	if (deepest_siblings.size() < 2) { cout << "improper length" << endl; }

	//Get the deepest two of the siblings
	Node *T2_a1 = deepest_siblings[0];
	Node *T2_a2 = deepest_siblings[1];
	#ifdef DEBUG_APPROX
	cout << "a1: " << T2_a1->str() << " a2: " << T2_a2->str() << endl;
	#endif

	bool cut_a1   = false;
	bool cut_a1_p = false;
	bool cut_a2   = false;
	bool cut_a2_p = false;

	if (T1_sibling_group->get_children().size() == 2) {
	  /*
	    7.1 case
	    Cut a1, pa1, a2 in F2, add 3 to num_cut
	    Consider adding 3 regardless if we actually cut 3,
	  */
	  cut_a1   = true;
	  cut_a1_p = true;
	  cut_a2   = true;
	  //num_cut += 3;
	  /*
	    7.2 case
	    Cut a1, a2, pa1, pa2, add 4 to num_cuts
	  */
	  if (previous_group == T1_sibling_group) {
	    cut_a2_p = true;
	    #ifdef DEBUG_APPROX
	    cout <<"Case 7.2" << endl;
	    #endif
	    //num_cut += 1;
	  }
	  else {
	    #if DEBUG_APPROX
	    cout << "Case 7.1" << endl;
	    #endif
	  }
	} // size == 2

	else if (T1_sibling_group->get_children().size() > 2) {
	  /*
	    7.3 case
	    If a2's parent's only sibling is part of the sibling group,
	    and a1's parent is a root or has a sibling that is not part of the sibling group
	    then cut a2 and a2_p otherwise a1 and a1_p
	  */
	  if (previous_group != T1_sibling_group) {
	    Node* T2_a2_p = T2_a2->parent();
	    bool x_2 = false;
	    bool a2_p_one_sibling = (T2_a2_p != NULL) &&
	      (T2_a2_p->parent() != NULL) &&
	      (T2_a2_p->parent()->get_children().size() == 2);
	    if (a2_p_one_sibling) {
	      list<Node *> group = T1_sibling_group->get_children();
	      //get the other one
	      Node *a2_p_sibling = T2_a2_p->parent()->get_children().front() == T2_a2_p ?
		T2_a2_p->parent()->get_children().back() :
		T2_a2_p->parent()->get_children().front();
	      //check if it is part of sibling group
	      bool a2_p_sibling_in_group = descendant_count[a2_p_sibling->get_preorder_number()] == -1;
	      if (a2_p_sibling_in_group) {
		Node* T2_a1_p = T2_a1->parent();
		bool a1_p_is_root = T2_a1_p->parent() == NULL;
		bool a1_p_sibling_not_in_group = false;
		if (!a1_p_is_root) {
		  list<Node*> a1_p_siblings = T2_a1_p->parent()->get_children();
		  for (list<Node*>::iterator i = a1_p_siblings.begin(); i != a1_p_siblings.end(); i++) {
		    if (descendant_count[(*i)->get_preorder_number()] != -1) {
		      a1_p_sibling_not_in_group = true;
		      break;
		    }
		  }
		}
		if (a2_p_one_sibling && a2_p_sibling_in_group && (a1_p_is_root || a1_p_sibling_not_in_group)) {
		  x_2 = true;

		}
	      }
	    }
	    #ifdef DEBUG_APPROX
	    cout << "Case 7.3" << endl;
	    #endif
	    //num_cut += 2;
	    if (x_2){

	      cut_a2   = true;
	      cut_a2_p = true;
	    }
	    else {
	      cut_a1   = true;
	      cut_a1_p = true;
	    }
	  }
	  /*7.4 case
	    cut a1 and a1_p
	  */
	  else if (previous_group == T1_sibling_group) {
	    #ifdef DEBUG_APPROX
	    cout << "Case 7.4" << endl;
	    #endif
	    cut_a1   = true;
	    cut_a1_p = true;
	    //num_cut += 2;
	  }
	}
	/*
	  Cutting section
	*/
	Node* T2_a2_p = NULL;
	if (cut_a1) {
	  Node *T2_a1_p = T2_a1->parent();
	  //singletons? components?
	  if (T2_a1_p != NULL) {
	    //Cut connections
	    T2_a1->cut_parent();
	    num_cut++;
	    //add as components
	    T2->add_component(T2_a1);
	    //Just cut 1 of two children of a1 parent
	    if (T2_a1_p->get_children().size() == 1) {
	      if (T2_a1_p->parent() == NULL) {
		T2_a1_p->contract(true);
		if (T2_a1_p->is_singleton() && T2_a1_p != T2->get_component(0)) {
		  singletons->push_front(T2_a1_p);
		}
	      }
	      else {
		Node* T2_b1 = T2_a1_p->get_children().front();
		T2_a1_p->contract(true);
		T2_a1_p = T2_b1; // <------------------------------------
	      }
	    }
	    //Check for singletons
	    if (T2_a1->is_singleton()) // wont ever be C0?
	      singletons->push_front(T2_a1);


	    if (cut_a1_p) {
	      if (T2_a1_p->parent() != NULL) {
		bool aborted_a2 = false;
		if (T2_a1_p->parent() == T2_a2->parent() &&
		    T2_a2->parent()->get_children().size() == 2)
		  {
		    cut_a2 = false; //cutting a1_p will cause a2 to get contracted up, so there is no more a2 to cut
		    cut_a2_p = true; //Instead we cut a2_p
		    aborted_a2 = true;
		  }
		Node *T2_a1_gp = T2_a1_p->parent();
		//Cut connections
		T2_a1_p->cut_parent();
		num_cut++;
		//add as components
		T2->add_component(T2_a1_p);
		//Just cut 1 of two children of a1 parent
		if (T2_a1_gp->get_children().size() == 1) {
		  if (T2_a1_gp->parent() == NULL) {
		    T2_a1_gp->contract(true);
		  }
		  else {
		    Node* T2_b1_p = T2_a1_gp->get_children().front();
		    T2_a1_gp->contract(true);
		    T2_a1_gp = T2_b1_p;
		  }
		  if (aborted_a2) { T2_a2_p = T2_a1_gp; }
		  if (T2_a1_gp != NULL && T2_a1_gp->is_singleton() && (T2_a1_gp != T2->get_component(0) || aborted_a2)) {
		    singletons->push_front(T2_a1_gp);
		  }
		}
		//Check for singletons
		if (T2_a1_p->is_singleton() && T2_a1_p != T2->get_component(0)) {
		  singletons->push_front(T2_a1_p);
		}
	      }
	    }
	  }//T2_a1_p() != NULL
	}
	if (cut_a2){
	  T2_a2_p = T2_a2->parent();
	  //could have cut a2's parent in previous steps, so could be singleton now
	  if (T2_a2->is_leaf() && T2_a2_p == NULL) {
	    singletons->push_front(T2_a2);
	  }
	  else if (T2_a2_p != NULL) {
	    //Cut connections
	    T2_a2->cut_parent();
	    num_cut++;
	    //add as components
	    T2->add_component(T2_a2);
	    //Just cut 1 of two children of a2 parent
	    if (T2_a2_p->get_children().size() == 1) {
	      if (T2_a2_p->parent() == NULL) {
		T2_a2_p->contract(true);
		if (T2_a2_p->is_singleton() && T2_a2_p != T2->get_component(0)) {
		  singletons->push_front(T2_a2_p);
		}
	      }
	      else {
		Node* T2_b2 = T2_a2_p->get_children().front();
		T2_a2_p->contract(true);
		T2_a2_p = T2_b2;
	      }
	    }
	    //Check for singletons
	    if (T2_a2->is_singleton())
	      singletons->push_front(T2_a2);

	  }
	}
	if (cut_a2_p) {
	  if (T2_a2_p != NULL && T2_a2_p->parent() != NULL) {
	    Node *T2_a2_gp = T2_a2_p->parent();
	    //Cut connections
	    T2_a2_p->cut_parent();
	    	    num_cut++;
	    //add as components
	    T2->add_component(T2_a2_p);
	    //Just cut 1 of two children of a2 parent
	    if (T2_a2_gp->get_children().size() == 1) {
	      if (T2_a2_gp->parent() == NULL) {
		T2_a2_gp->contract(true);
	      }
	      else {
		Node* T2_b2_p = T2_a2_gp->get_children().front();
		T2_a2_gp->contract(true);
		T2_a2_gp = T2_b2_p;
	      }

	      if (T2_a2_gp->is_singleton() && T2_a2_gp != T2->get_component(0)) {
		singletons->push_front(T2_a2_gp);
	      }
	    }
	    //Check for singletons
	    if (T2_a2_p->is_singleton()) {
	      singletons->push_front(T2_a2_p);
	    }

	  }
	}


      }//else
      //delete identical_sibling_groups;
      previous_group = T1_sibling_group;
    }//!sibling_groups->empty()

  } //while(!sibling_groups->empty() && !singletons->empty()
  // if the first component of the forests differ then we have cut p
  /*
  if (T1->get_component(0)->get_twin() != T2->get_component(0)) {
    if (!T1->contains_rho()) {
      T1->add_rho();
      T2->add_rho();
    }
    else
      // hack to ignore rho when it shouldn't be in a cluster
      num_cut -=3;
  }
  */
  if (save_forests) {
    *F1 = new Forest(T1);
    *F2 = new Forest(T2);
  }
  //if (num_cut) num_cut--;
  return num_cut;
}

#ifdef DEBUG_CASE_COUNTER
struct mult_case_counter {
  int case_71, case_72, case_73, case_74;
  int case_81, case_82, case_83, case_84;
  int case_85, case_86, case_87;
};
static mult_case_counter case_counter = {};
void print_mult_case_count() {
  cout << "Stats:" << endl;
  cout << "\tCase 7.1:" << case_counter.case_71 << endl;
  cout << "\tCase 7.2:" << case_counter.case_72 << endl;
  cout << "\tCase 7.3:" << case_counter.case_73 << endl;
  cout << "\tCase 7.4:" << case_counter.case_74 << endl;
  cout << "\tCase 8.1:" << case_counter.case_81 << endl;
  cout << "\tCase 8.2:" << case_counter.case_82 << endl;
  cout << "\tCase 8.3:" << case_counter.case_83 << endl;
  cout << "\tCase 8.4:" << case_counter.case_84 << endl;
  cout << "\tCase 8.5:" << case_counter.case_85 << endl;
  cout << "\tCase 8.6:" << case_counter.case_86 << endl;
  cout << "\tCase 8.7:" << case_counter.case_87 << endl;
}
  
#endif

int rSPR_branch_and_bound_mult_range(Forest *T1, Forest *T2, int start_k){
  return rSPR_branch_and_bound_mult_range(T1, T2, start_k, MAX_SPR);
}
int rSPR_branch_and_bound_mult_range(Forest *T1, Forest *T2, int start_k, int end_k){
  int exact_spr = -1;
  int k;
  for (k = start_k; k <= end_k; k++) {
    Forest *F1 = new Forest(T1);
    Forest *F2 = new Forest(T2);
    if (!sync_twins(F1,F2)) {
      exact_spr = 0;
      continue;
    }
    if (LEAF_REDUCTION) {
      reduction_leaf_mult(F1, F2);
    }
    #ifdef DEBUG
    cout << "Trying K = " << k << endl << "------------------" << endl;
    #else
    cout << k << " " << endl;
    #endif
    exact_spr = rSPR_branch_and_bound_mult(F1, F2, k);
    #ifdef DEBUG
    cout << "------------------" << endl;
    cout << "Finished K = " << k << " return value : " << exact_spr << endl;   
    #endif
    if (exact_spr >= 0) {
      F1->swap(T1);
      F2->swap(T2);
    }
    delete F1;
    delete F2;

    if (exact_spr >= 0) {    
      break;
    }
  }
  #ifdef DEBUG_CASE_COUNTER
  print_mult_case_count();
  #endif
  if (k > end_k) {
    k = -1;
  }
  return k;
}
int rSPR_branch_and_bound_mult(Forest *T1, Forest *T2, int k){

  if (!sync_twins(T1,T2)) {
    return 0;      
  }
  T2->max_preorder = T2->components[0]->get_max_preorder_number(0);//preorder_number(0);          
  list<Node *> *sibling_groups = T1->find_sibling_groups();
  //sibling_groups->push_front(sibling_groups->back());
  //sibling_groups->pop_back();
  list<Node *> singletons     = T1->find_singletons();
  
  list<pair<Forest,Forest>> AFs = list<pair<Forest,Forest>>();
  //list<Node *> protected_stack = list<Node*>();
  int num_ties = 2;
  int final_k = rSPR_branch_and_bound_mult_hlpr(T1, T2, k, sibling_groups, &singletons, NULL, &AFs, &num_ties);  

  //print AFs
  if (!AFs.empty() && final_k > -1) {
    if (ALL_MAFS) {
      cout << endl << endl << "FOUND ANSWER" << endl;
      // TODO: this is a cheap hack
      for (list<pair<Forest,Forest> >::iterator x = AFs.begin(); x != AFs.end(); x++) {
	cout << "\tT1: ";
	x->first.print_components();
	cout << "\tT2: ";
	x->second.print_components();
      }
    }
      AFs.front().first.swap(T1);
  AFs.front().second.swap(T2);
  sync_twins(T1,T2);

  }

  
  delete sibling_groups;
  if (final_k >= 0)
  return k - final_k;
  else
    return final_k;
    
}


class bb_mult_recurse_data {
public:
  Forest *T1;
  Forest *T2;
  list<Node*> *sibling_groups;
  list<Node*> *singletons;
  map<Node *, Node*> node_map;
};

/* Generates a copy of the data used to recurse on rSPR_branch_and_bound_mult_hlpr so 
   that original forests aren't clobbered,
   updates pointers accordingly */
//TODO: 1 new not 5
bb_mult_recurse_data *generate_mult_recurse_data(Forest *T1, Forest *T2, list<Node*> *sibling_groups, list<Node*> *singletons) {

  bb_mult_recurse_data *data = new bb_mult_recurse_data();

  data->node_map = map<Node*, Node*>();
  data->T1 = new Forest(T1, &data->node_map);
  data->T2 = new Forest(T2, &data->node_map);
  //TODO: smarter way of syncing
  //for (auto n = node_map.begin(); n != node_map.end(); n++) {
    //(*n).first->set_twin(node_map[(*n).first->get_twin()]);
    /*
      if ((*n).second->is_leaf()) {
      (*n).second->set_twin(node_map[(*n).first->get_twin()]);
      node_map[(*n).first->get_twin()]->set_twin((*n).second);
      cout << "Synced twin : " << (*n).second->str() << " -> " << (*n).second->get_twin()->str() << endl;}*/
  //}
  sync_twins(data->T1, data->T2);
  data->sibling_groups = new list<Node*>();
  for (auto n = sibling_groups->begin(); n != sibling_groups->end(); n++) {
    data->sibling_groups->push_back(data->node_map[*n]);
  }
  data->singletons = new list<Node*>();
  for (auto n = singletons->begin(); n != singletons->end(); n++) {
    data->singletons->push_back(data->node_map[*n]);
  }
  return data;
}
//Cuts to_cut, adds to components, conditionally adds to singletons
//Assumes parent is not null and no Null parameters
void mult_cut_and_cleanup(Node* to_cut, Forest *T2, list<Node*> *singletons) {
  //Node* T2_a1_next = next_data->node_map[T2_a1];
  Node* to_cut_p = to_cut->parent();
  //Cut connections  
  to_cut->cut_parent();
  //add as components
  T2->add_component(to_cut);
  //Just cut 1 of two children of a1 parent
  if (to_cut_p->get_children().size() == 1) {
    if (to_cut_p->parent() == NULL) {
      to_cut_p->contract(true);
      if (to_cut_p->is_singleton() && to_cut_p != T2->get_component(0)) {
      singletons->push_front(to_cut_p);
      }
    }
    else {
      Node* to_cut_b = to_cut_p->get_children().front();
      to_cut_p->contract(true);
      if (to_cut_b->is_singleton() && to_cut_b != T2->get_component(0)) {
	singletons->push_front(to_cut_b);
      }
    }
  }
  //Check for singletons
  if (to_cut->is_singleton())
    singletons->push_front(to_cut);  
}
//cuts everything except node, possibly expanding, adds to components, conditionally adds to singletons
//Assumes parent is not null and no Null parameters
//NOTE: potentially unsafe for preorder numbers
void mult_cut_all_except_and_cleanup(Node* T2_a1, Forest *T2, list<Node*> *singletons) {
  Node* T2_b1;
  Node* parent = T2_a1->parent();
  if (parent->get_children().size() == 2) {
    T2_b1 = parent->get_children().front() == T2_a1 ?
      parent->get_children().back() :
      parent->get_children().front();		
  }
  else {
    list<Node*> all_but_a1 = list<Node*>(parent->get_children());
    all_but_a1.remove(T2_a1);	    
    T2_b1 = parent->expand_children_out(all_but_a1);
    T2_b1->set_preorder_number(parent->get_preorder_number());
  }
  //hack for now
  //We immediately contract a1 up so we know the parent's preorder number is available
  //Cut connections

  T2_b1->cut_parent();

  //add as components
  T2->add_component(T2_b1);
	    
  //Just cut 1 of two children of a2 parent (will always be in this case?)
  if (parent->get_children().size() == 1) {
    if (parent->parent() == NULL) {
      parent->contract(true);
      if (parent->is_singleton() && parent != T2->get_component(0)) {
	singletons->push_front(parent);
      }
    }
    else {
      parent = parent->contract(true);
      if (T2_a1->is_singleton() && T2_a1 != T2->get_component(0))
	singletons->push_front(T2_a1);
    }
  }
  if (T2_b1->is_singleton())
    singletons->push_front(T2_b1);
}

//Adds rho and a certain singleton to cut. Basically its for realising when we have cut rho but it is already a singleton so we explicitly continue on 
#define MULT_RHO_CUT_AND_RESOLVE(rho_to_singleton, node_to_protect) {			\
  bb_mult_recurse_data *next_data = generate_mult_recurse_data(T1, T2, sibling_groups, singletons); \
  Node* protect = NULL;							\
  next_data->T1->add_rho();						\
  next_data->T2->add_rho();						\
  next_data->singletons->push_front(next_data->node_map[rho_to_singleton]);	\
  int result_k = rSPR_branch_and_bound_mult_hlpr(next_data->T1, next_data->T2, k - 1, next_data->sibling_groups, next_data->singletons, protect, AFs, num_ties); \
  delete next_data->T1;							\
  delete next_data->T2;						\
  delete next_data->sibling_groups;				\
  delete next_data->singletons;					\
  delete next_data;						\
  if (result_k > best_k) {					\
    best_k = result_k;						\
  }								\
  }								\

//TODO (Ben): try using a routine instead of a macro, I suspect it will slow down because of
//            stack and parameter passing though
#define MULT_BB_CUT_AND_RESOLVE(nodes_to_cut, nodes_to_exclude_cutting, node_to_protect) { \
  bb_mult_recurse_data *next_data = generate_mult_recurse_data(T1, T2, sibling_groups, singletons);\
  int num_cuts = 0;							\
  for (int i = 0; i < nodes_to_cut.size(); i++) {			\
    Node* T2_ax_next = next_data->node_map[nodes_to_cut[i]];		\
    mult_cut_and_cleanup(T2_ax_next, next_data->T2, next_data->singletons);\
    num_cuts++;								\
  }									\
  for (int i = 0; i < nodes_to_exclude_cutting.size(); i++) {		\
    Node* T2_ax_next = next_data->node_map[nodes_to_exclude_cutting[i]]; \
    mult_cut_all_except_and_cleanup(T2_ax_next, next_data->T2, next_data->singletons); \
    num_cuts++;								\
  }									\
  Node* protect = next_data->node_map[node_to_protect];			\
  int result_k = rSPR_branch_and_bound_mult_hlpr(next_data->T1, next_data->T2, k - num_cuts, next_data->sibling_groups, next_data->singletons, protect, AFs, num_ties); \
  delete next_data->T1;							\
  delete next_data->T2;							\
  delete next_data->sibling_groups;					\
  delete next_data->singletons;						\
  delete next_data;							\
  if (result_k > best_k) {						\
    best_k = result_k;							\
    if (!ALL_MAFS && best_k > -1) {					\
    }									\
  }									\
}

//TODO: UndoMachine, then cleanup all constructors, bb_mult_recurse_data relying on copies of trees
int rSPR_branch_and_bound_mult_hlpr(Forest *T1, Forest *T2,
				    int k,
				    list<Node*> *sibling_groups, list<Node*> *singletons,
				    Node *protected_node, list<pair<Forest,Forest>> *AFs,
				    int* num_ties) {
  //run out of cuts if k < 0
  if (k < 0) {
    return k;
  }
  Node* previous_group = sibling_groups->back();
  int best_k = -1;
  while(!singletons->empty() || !sibling_groups->empty()) {
	  
    // Case 1 - Remove singletons
    while(!singletons->empty()) {
      
      Node *T2_a = singletons->back();
      singletons->pop_back();
      #ifdef DEBUG
      cout << "Handling singleton: " << T2_a->str() << endl;
      #endif
      
      Node *T1_a = T2_a->get_twin();
      Node *T1_a_p = T1_a->parent();
      
      if (T1_a_p == NULL)
	continue;      

      // find twin in T1
      //If we have added component 0, then we have made a B cut and added rho
      if (T2_a == T2->get_component(0)){// && T1_a != T1->get_component(0)) {
	if (!T1->contains_rho()) {
	  T1->add_rho();
	  T2->add_rho();
	  k--;
	  //continue;
	}
      }

      bool is_sibling_group = T1_a_p->is_sibling_group();
      // cut the edge above T1_a
      T1_a->cut_parent();
      if (!T1_a->is_leaf()) {
	T1_a_p->decrement_non_leaf_children();
      }
      T1->add_component(T1_a);
      
      //only contract if one node
      if (T1_a_p->get_children().size() == 1) {
	//If we contract this node, and it used to be a sibling group
	//then it is not a sibling group anymore
	if (is_sibling_group) {
	  sibling_groups->remove(T1_a_p);
          #ifdef DEBUG
	  cout << "Removed " << T1_a_p->str() << " from sibling groups" << endl;
	  #endif
	}
	
	Node *possible_previous_sibling = T1_a_p->get_children().front();
	bool was_sibling_group = possible_previous_sibling->is_sibling_group();
	
	Node *T1_a_gp = T1_a_p->parent();
	Node *T1_new_a_p = T1_a_p;
	
	//If the child is a sibling group, there is a possibility contract()
	//will delete it, so we need to update it in the sibling groups
	if (T1_a_p->parent() == NULL) {
	  if (was_sibling_group){
	    list<Node*>::iterator i = find(sibling_groups->begin(), sibling_groups->end(), possible_previous_sibling);
	    *i = T1_new_a_p;	    
	  }
	}
	else {
	  T1_new_a_p = possible_previous_sibling;
	}
       
        T1_a_p->contract(true);
	T1_new_a_p->recalculate_non_leaf_children();
	//After contracting this, the grandparent may be a sibling group now
	if (T1_a_gp != NULL) {
	  T1_a_gp->recalculate_non_leaf_children(); //can we tell what this would be instead of recalculating?
	  if (T1_a_gp->is_sibling_group()) {	 
	    sibling_groups->push_front(T1_a_gp);	 
	  }
	}
      }

    }//!singletons->empty()

    //NOTE: we know there are no singletons left here
    if(!sibling_groups->empty()) {
      //Find identical sibling groups
      //Get the first group that has identical sibling groups, otherwise default to the group on the back
      list<Node*>::reverse_iterator i = sibling_groups->rbegin();
      Node *T1_sibling_group = sibling_groups->back();
      list<list<Node*>> identical_sibling_groups;
      T1_sibling_group->find_identical_sibling_groups(&identical_sibling_groups);
      bool found_identical = false;
      for (; i != sibling_groups->rend(); i++ ){
	(*i)->find_identical_sibling_groups(&identical_sibling_groups);
	if (identical_sibling_groups.size() > 0) {
	  T1_sibling_group = (*i);
	  found_identical = true;
	  break;
	}
      }
      #ifdef MULT_PICK_LARGEST_GROUP
      if (!found_identical) {
	int max_size = 0;
	Node* largest_group = NULL;
	for (auto i = sibling_groups->begin(); i != sibling_groups->end(); i++) {
	  if ((*i)->get_children().size() > max_size) {
	    max_size = (*i)->get_children().size();
	    largest_group = (*i);
	  }
	}
	T1_sibling_group = largest_group;
      }
      #endif
      #ifdef DEBUG
      cout << "K = " << k << endl;
      cout << "F2: ";
      T2->print_components();
      cout << endl;
      cout << "F1: ";
      T1->print_components();
      cout << endl;
      #endif

      /* 
	 Case where a subset of the group have the same parent both in T1 and T2
      */
      // Case 2 - Contract identical sibling pair
      if (identical_sibling_groups.size() > 0) {	  		
	list<list<Node *>>::iterator i;
	for (i = identical_sibling_groups.begin(); i != identical_sibling_groups.end(); i++) {
	  //Contract the groups
	  list<Node *> T2_group = (*i);
	  Node *T2_p = T2_group.front()->parent();
	  #ifdef DEBUG
	  //cout << "Contracting T1... " << endl;
	  #endif
	  Node *T1_group_new = T1_sibling_group->contract_twin_group(&T2_group);
	  #ifdef DEBUG
	  //cout << "Contracting T2... " << endl;
	  #endif
	  Node *T2_group_new = T2_p->contract_sibling_group(&T2_group);

	  //Maintain twins
	  T1_group_new->set_twin(T2_group_new);
	  T2_group_new->set_twin(T1_group_new);			

	  // check if T2_p is a singleton after the contraction
	  if (T2_p->is_singleton() && T2_p != T2->components[0]){
	    singletons->push_front(T2_p);
	  }
	  if (T1_sibling_group->parent() != NULL) {
	    //Check if the contraction made a new sibling group
	    T1_sibling_group->parent()->recalculate_non_leaf_children();
	    if (T1_sibling_group->parent()->is_sibling_group()) {
	      sibling_groups->push_front(T1_sibling_group->parent());
	      #ifdef DEBUG
	      //cout << "Added new sibling group after contraction: " << T1_sibling_group->parent()->str_subtree() << endl;
	      #endif
	    }
	  }
	  if (!T1_sibling_group->is_sibling_group()) {
	    sibling_groups->remove(T1_sibling_group);
	      #ifdef DEBUG
	    //cout << "Removed new sibling group after contraction: " << T1_sibling_group->str() << endl;
	      #endif

	  }	  
	}
      }

      /*
	4 branching case
	Step 6-8 in paper
	Cut above a1, b1, a2, b2
	Part 1: Get the LCA, this is the numbering to the top part
	Part 2: Get subset of sibling group that is descendant of this LCA
	Part 3: Sort them based on depth
	Part 4: Since this is the approximation, we cut deepest 2
      */

      // Case 3
      else {

	//Check branch and bound
	if (BB) {
	  //copies for the approx so we dont clobber this tree
	  map<Node*, Node*> approx_map = map<Node*, Node*>();
	  Forest T1_approx = Forest(T1, &approx_map);
	  Forest T2_approx = Forest(T2, &approx_map);
	  sync_twins(&T1_approx, &T2_approx);
	  list<Node*> sibling_group_approx = list<Node*>();
	  for (auto n = sibling_groups->begin(); n != sibling_groups->end(); n++) {
	    sibling_group_approx.push_back(approx_map[*n]);
	  }
	  list<Node*> singletons_approx = list<Node*>();
	  for (auto n = singletons->begin(); n != singletons->end(); n++) {
	    singletons_approx.push_back(approx_map[*n]);
	  }
	  int approx_spr = rSPR_worse_3_mult_approx_hlpr(&T1_approx, &T2_approx, &singletons_approx, &sibling_group_approx, NULL, NULL, false);
	  //TODO: Sometimes approx returns 1 over the right amount. For example 4 for a 1 cut tree or 16 for a 3 cut tree
	  //If approx_spr > 3k then we will not have enough cuts
	  if (approx_spr > 3*k + 1) {
#ifdef DEBUG
	    cout << "approx failed approx k = " << approx_spr  <<  endl;
#endif
	    return -1;
	  }
	}

	//Finding deepest siblings
	#ifdef DEBUG
	cout << "Sibling group to be cutting: " << T1_sibling_group->str_subtree() << endl;
	#endif
	//vector of ints describing how many of the siblings are in its descendants, indexed by preorder number
	vector<int> descendant_count = T1_sibling_group->find_pseudo_lca_descendant_count(T2->max_preorder + 1);
	Node* arbitrary_lca = T1_sibling_group->find_arbitrary_lca(T2->components, descendant_count);
	vector<Node *> deepest_siblings;
	vector<vector<Node*>> siblings_by_depth;
	//Get depth of siblings from root of each component
	//If the lca is null, all siblings are in different components, ie no path between them
	if (arbitrary_lca == NULL) {	
	  siblings_by_depth = vector<vector<Node *>>(10);
	  for (int i = 0; i != T2->components.size(); i++) {
	    //if preorder is -1 then this is rho
	    if (T2->components[i]->get_preorder_number() == -1) { continue; }
	    T2->components[i]->get_deepest_siblings(descendant_count, siblings_by_depth);
	  }
	}
	//Otherwise they share an LCA
	else {
	  siblings_by_depth = arbitrary_lca->get_deepest_siblings(descendant_count);
	}
	map<Node*, int> s_map = map<Node*, int>();
	//Get deepest siblings returns sparse vector, this compacts it
	deepest_siblings = contract_deepest_siblings(siblings_by_depth, &s_map);      
	// Should assert here
	if (deepest_siblings.size() < 2) { cout << "improper length" << endl; }

	//Get the deepest two of the siblings
	Node *T2_a1 = deepest_siblings[0];
	Node *T2_a2 = deepest_siblings[1];
	#ifdef DEBUG
	cout << "a1: " << T2_a1->str() << " a2: " << T2_a2->str() << endl;
	#endif
	best_k = -1;

	//MULT_4_BRANCH is naive approach of making 4 cuts every time
	if (!MULT_4_BRANCH) {
	//step 7
	/*
	  if a1 == prot, x = 2 else x = 1
	*/
	if (protected_node != NULL) {
	  Node* T2_ax;
	  if (T2_a1 != protected_node) {
	    T2_ax = T2_a1;
	  }
	  else {
	    T2_ax = T2_a2;
	  }
	  
	  bool a0_descendant_of_lca = false;
	  Node* stepper = protected_node;
	  while (stepper->parent() != NULL && stepper->parent() != arbitrary_lca) {
	    stepper = stepper->parent();
	  }
	  if (stepper->parent() == arbitrary_lca) {
	    a0_descendant_of_lca = true;
	  }

	  //TODO (Ben) do one iteration of pl parent for has_aj_child and has_non_aj_child,
	  //           as well as any other flags, before step 7 check
	  bool pl_has_aj_child = false;
	  if (arbitrary_lca != NULL) {
	    Node* pl = arbitrary_lca->parent();	  
	    if (pl != NULL) {
	      for (auto i = pl->get_children().begin(); i != pl->get_children().end(); i++) {
		if (*i != arbitrary_lca && descendant_count[(*i)->get_preorder_number()] == -1) {
		  pl_has_aj_child = true;
		  break;
		}
	      }
	    }
	  }
	  //step 7.1
	  if (arbitrary_lca == NULL) {
	    #ifdef DEBUG
	    cout << "Case 7.1 T2_ax: " << T2_ax->str() << endl;
	    #endif
	    #ifdef DEBUG_CASE_COUNTER
	    case_counter.case_71++;
	    #endif
	    /* 
	       recurse on cut ax prot protected node
	    */
	    //7.1 : cut ax
	    if (T2_ax->parent() != NULL) {
	      vector<Node*> to_cut = {T2_ax};
	      vector<Node*> to_cut_except = {};
	      MULT_BB_CUT_AND_RESOLVE(to_cut, to_cut_except, protected_node);
	    }
	  }

	  //step 7.2
	  else if (!a0_descendant_of_lca) {
	    /*
	      recurse on cut a1's B's up to LCA prot protected node

	    */
	    #ifdef DEBUG
	    cout << "Case 7.2a cut all B1's Protected Node: " << protected_node->str() <<  endl;
	    #endif
	    #ifdef DEBUG_CASE_COUNTER
	    case_counter.case_72++;
	    #endif

	    {
	      vector<Node*> to_cut = {};
	      vector<Node*> to_cut_except = {};
	      Node* stepper = T2_a1;
		
	      while(stepper->parent() != arbitrary_lca) { 
		to_cut_except.push_back(stepper);
		stepper = stepper->parent();
	      }
	      //TODO: use list to push_front or figure out how to add from top to bottom 
	      reverse(to_cut_except.begin(), to_cut_except.end());
	      MULT_BB_CUT_AND_RESOLVE(to_cut, to_cut_except, protected_node);
	    }
	    
	    if (s_map[T2_a1] > 1 ||
		s_map[deepest_siblings[deepest_siblings.size()-1]] > 0){
	      /*
		recurse on cut a1 prot protected node
	      */
#ifdef DEBUG
	    cout << "Case 7.2b Cut a1 Protected Node: " << protected_node->str() <<  endl;
#endif
	      //7.2b cut a1
	      if (T2_a1->parent() != NULL) {
		vector<Node*> to_cut = {T2_a1};
		vector<Node*> to_cut_except = {};
		MULT_BB_CUT_AND_RESOLVE(to_cut, to_cut_except, protected_node);
	      }     
	    }
	     /*

	      Not part of the outlined special cases

	     */
	    //if (T2_a1->parent()->get_children().size() == 2)
	    {
	      #ifdef DEBUG
	      cout << "Case 7.2c Cut a2 T2_ax: " << T2_ax->str() << " Protected Node: " << protected_node->str() << endl;
	      #endif
	      vector<Node*> to_cut = {T2_a2};
	      vector<Node*> to_cut_except = {};
	      MULT_BB_CUT_AND_RESOLVE(to_cut, to_cut_except, protected_node);	      
	    }
	  }	
	    

	  //TODO: Is component 0 considered a root?
	  //If we do get component 0 can we consider this if rho has been added?
	  //step 7.3
	  else if (deepest_siblings.size() == 2 &&
		   T2_ax->parent() == arbitrary_lca &&
		   a0_descendant_of_lca &&
		   (
		    (arbitrary_lca->parent() == NULL && arbitrary_lca != T2->get_component(0)) ||
		    pl_has_aj_child
		    )) {
	    /*
	      recurse on cut Bx prot protected node
	    */
	    //7.3 cut Bx

	    #ifdef DEBUG
	    cout << "Case 7.3 Cut T2_bx T2_ax: " << T2_ax->str() << " Protected Node: " << protected_node->str() << endl;
	    #endif
	    #ifdef DEBUG_CASE_COUNTER
	    case_counter.case_73++;
	    #endif

	    {
	      vector<Node*> to_cut = {};
	      vector<Node*> to_cut_except = {T2_ax};
	      MULT_BB_CUT_AND_RESOLVE(to_cut, to_cut_except, protected_node);
	    }
	    /*

	      This is NOT in the cases. Figure out whats happening with 8.2 with r = 2

	     */
	    /*
	    {
	      cout << "Case 7.3 Cut T2_ax: " << T2_ax->str() << " Protected Node: " << protected_node->str() << endl;
	      vector<Node*> to_cut = {T2_ax};
	      vector<Node*> to_cut_except = {};
	      MULT_BB_CUT_AND_RESOLVE(to_cut, to_cut_except, protected_node);
	      }
	    */
	  }


	  //step 7.4
	  else if (a0_descendant_of_lca) {
	    /*
	      recurse on cut ax prot protected_node
	                 if m > 2 and r > 2 recurse on 
			    cut all Bx's then B`x
	     */
	    //7.4 cut ax
	    #ifdef DEBUG
	    cout << "Case 7.4 T2_ax: " << T2_ax->str() << " Protected Node: " << protected_node->str() << endl;
	    #endif
	    #ifdef DEBUG_CASE_COUNTER
	    case_counter.case_74++;
	    #endif

	    {
	      vector<Node*> to_cut = {T2_ax};
	      vector<Node*> to_cut_except = {};
	      MULT_BB_CUT_AND_RESOLVE(to_cut, to_cut_except, protected_node);
	    }
	    if (T1_sibling_group->get_children().size() > 2 && deepest_siblings.size() > 2) {
	      #ifdef DEBUG
	      cout << "Case 7.4b T2_ax: " << T2_ax->str() << " Protected Node: " << protected_node->str() << endl;
	      #endif
	      //7.4 cut all Bx B`x
	      vector<Node*> to_cut = {};
	      vector<Node*> to_cut_except = {};
	      Node* stepper = T2_ax;
		
	      while(stepper->parent() != arbitrary_lca) { //Do we know ax is descendant of lca?
		to_cut_except.push_back(stepper);
		stepper = stepper->parent();
	      }
	      //TODO: use list to push_front or figure out how to add from top to bottom 
	      reverse(to_cut_except.begin(), to_cut_except.end());
	      to_cut_except.push_back(T2_ax);
	      MULT_BB_CUT_AND_RESOLVE(to_cut, to_cut_except, protected_node);	      
	    }
	    /*

	      Not part of the outlined special cases

	     */
	    {
	      vector<Node*> to_cut = {T2_a2};
	      vector<Node*> to_cut_except = {};
	      MULT_BB_CUT_AND_RESOLVE(to_cut, to_cut_except, protected_node);	     
	    }
	    {
	      vector<Node*> to_cut = {};
	      vector<Node*> to_cut_except = {T2_a1};
	      MULT_BB_CUT_AND_RESOLVE(to_cut, to_cut_except, protected_node);	     
	    }

	    /*
	    if (T1_sibling_group->get_children().size() == 2 &&
		deepest_siblings.size() == 2 &&
		T2_a1 == T2_ax &&
		T2_a2 != protected_node) {
	      vector<Node*> to_cut = {T2_a2};
	      vector<Node*> to_cut_except = {};
	      MULT_BB_CUT_AND_RESOLVE(to_cut, to_cut_except, protected_node);
	      }*/
	  }
	}

	//step 8
	else {
	  bool all_but_ar_s1 = true;
	  for (int i = 0; i < deepest_siblings.size() - 1; i++) {
	    if (s_map[deepest_siblings[i]] != 1) {
	      all_but_ar_s1 = false;
	      break;
	    }
	  }
	  bool all_s1 = s_map[deepest_siblings[deepest_siblings.size()-1]] == 1 &&
	    all_but_ar_s1;

	  bool lca_p_contains_sibling = false;	    
	  if (arbitrary_lca != NULL && arbitrary_lca->parent() != NULL) {
	    for (list<Node*>::iterator i = arbitrary_lca->parent()->get_children().begin();
		 i != arbitrary_lca->parent()->get_children().end();
		 i++) {
	      if (descendant_count[(*i)->get_preorder_number()] == -1) {
		lca_p_contains_sibling = true;
		break;
	      }
	    }
	  }

	  //step 8.1
	  if (arbitrary_lca == NULL) {
	    /*
	      recurse on cut a1 no prot,
	      cut a2 no prot
	    */
#ifdef DEBUG
	    cout << "Case 8.1a cut a1" << endl;
#endif
#ifdef DEBUG_CASE_COUNTER
	    case_counter.case_81++;
#endif

	    //8.1 : Cut a1
	    if (T2_a1->parent() != NULL) {
	      vector<Node*> to_cut = {T2_a1};
	      vector<Node*> to_cut_except = {};
	      MULT_BB_CUT_AND_RESOLVE(to_cut, to_cut_except, NULL);
	    }
	    else if (T2_a1 == T2->get_component(0)) {
	      if (!T1->contains_rho()) {
		MULT_RHO_CUT_AND_RESOLVE(T2_a1, NULL);
	      }
	    }
#ifdef DEBUG
	    cout << "Case 8.1b cut a2" << endl;
#endif
	    //8.1 : Cut a2
	    if (T2_a2->parent() != NULL) {
	      vector<Node*> to_cut = {T2_a2};
	      vector<Node*> to_cut_except = {};
	      MULT_BB_CUT_AND_RESOLVE(to_cut, to_cut_except, NULL);
	    }
	    else if (T2_a2 == T2->get_component(0)) {
	      if (!T1->contains_rho()) {
		MULT_RHO_CUT_AND_RESOLVE(T2_a2, NULL);
	      }
	    }
	  
	  }
	  //step 8.2
	  else if (all_but_ar_s1 &&
		   deepest_siblings[deepest_siblings.size()-1]->parent() == arbitrary_lca) {
	    /*
	      recurse on cut all B's except shallowest
	      for each sibling except for shallowest,
	      cut all other B's protect ai
		           
	    */
	    //8.2 B's
	    {
#ifdef DEBUG
	    cout << "Case 8.2a cut B1 - B(r-1)" << endl;
#endif
#ifdef DEBUG_CASE_COUNTER
	    case_counter.case_82++;
#endif

	      vector<Node*> to_cut = {};
	      vector<Node*> to_cut_except = {};
	      for (int i = 0; i < deepest_siblings.size() - 1; i++) {
		to_cut_except.push_back(deepest_siblings[i]);
	      }
	      MULT_BB_CUT_AND_RESOLVE(to_cut, to_cut_except, NULL);
	    }
	    //all other B's except ai
	    //Doesnt explicitly say this in the paper, but
	    //this would only happen if r > 2


	    if (deepest_siblings.size() > 2){
	      for (int i = 0; i < deepest_siblings.size() - 1; i++) {
#ifdef DEBUG
	      cout << "Case 8.2b for all ai cut all other B" << endl;
#endif

		vector<Node*> to_cut = {};
		vector<Node*> to_cut_except = {};
		for (int j = 0; j < deepest_siblings.size() - 1; j++) {
		  if (i != j) {
		    to_cut_except.push_back(deepest_siblings[j]);
		  }
		}
		MULT_BB_CUT_AND_RESOLVE(to_cut, to_cut_except, deepest_siblings[i]);
	      }
	      }
	    /*
	    else if (T1_sibling_group->parent() != NULL && T1_sibling_group->parent()->get_children().size() == 2) {
	      Node* gp = T1_sibling_group->parent();
	      Node* aunt = gp->get_children().front() == T1_sibling_group ?
		gp->get_children().back() :
		gp->get_children().front();
	      if (aunt->is_leaf()) {
		Node* a1_p = T2_a1->parent();
		vector<Node*> a1_p_leaves = a1_p->find_leaves();
		Node* T2_aunt = aunt->get_twin();
		for (int i = 0; i < a1_p_leaves.size(); i++) {
		  if (a1_p_leaves[i] == T2_aunt) {
#ifdef DEBUG
		    cout << "Case 8.2c cut a2" << endl;
#endif
		    vector<Node*> to_cut = {T2_a2};
		    vector<Node*> to_cut_except = {};
		    MULT_BB_CUT_AND_RESOLVE(to_cut, to_cut_except, NULL);
		    break;
		  }
		}

	      }
	      }*/


	      else if (deepest_siblings.size() == 2) {
		  vector<Node*> to_cut = {T2_a2};
		  vector<Node*> to_cut_except = {};
		  MULT_BB_CUT_AND_RESOLVE(to_cut, to_cut_except, NULL);
	      }

	  }
	  

	  //step 8.3
	  else if (T1_sibling_group->get_children().size() == 2 &&
		   s_map[T2_a1] + s_map[T2_a2] >= 2) {
	    /*recurse on cut a1 no prot,
	      cut a2 no prot,
	      cut all along a1 to a2 no prot
	    */
	    //8.3 : Cut a1
	    if (T2_a1->parent() != NULL) {
#ifdef DEBUG
	      cout << "Case 8.3a cut a1" << endl;
#endif
#ifdef DEBUG_CASE_COUNTER
	    case_counter.case_83++;
#endif

	      vector<Node*> to_cut = {T2_a1};
	      vector<Node*> to_cut_except = {};
	      MULT_BB_CUT_AND_RESOLVE(to_cut, to_cut_except, NULL);
	    }
	    //8.3 : Cut a2
	    if (T2_a2->parent() != NULL) {
#ifdef DEBUG
	    cout << "Case 8.3b cut a2" << endl;
#endif
	      vector<Node*> to_cut = {T2_a2};
	      vector<Node*> to_cut_except = {};
	      MULT_BB_CUT_AND_RESOLVE(to_cut, to_cut_except, NULL);
	    }
	    //8.3 : All along path from a1, a2
	    {
#ifdef DEBUG
	    cout << "Case 8.3c cut all B1's and B2's" << endl;
#endif
	      vector<Node*> to_cut = {};
	      vector<Node*> to_cut_except = {};
	      Node* stepper = T2_a1;
		
	      while(stepper->parent() != arbitrary_lca) { //no need to check for null! we know theres a path
		to_cut_except.push_back(stepper);
		stepper = stepper->parent();
	      }
	      stepper = T2_a2;
	      while(stepper->parent() != arbitrary_lca) { //no need to check for null! we know theres a path
		to_cut_except.push_back(stepper);
		stepper = stepper->parent();
	      }
	      /*
		Cut from top to bottom,
		Cutting from bottom to top causes the contraction to invalidate
		the node, since contraction is implemented by cutting the parent, then giving
		the child to the parent
	      */
	      //TODO: use list to push_front or figure out how to add from top to bottom 
	      reverse(to_cut_except.begin(), to_cut_except.end());
	      MULT_BB_CUT_AND_RESOLVE(to_cut, to_cut_except, NULL);
	    }	      
	  }

	  //step 8.4
	  else if (T1_sibling_group->get_children().size() > 2 &&
		   deepest_siblings.size() == 2 &&
		   s_map[T2_a1] == 1 &&
		   s_map[T2_a2] == 1) {
	    /* 
	       recurse on cut a1 and a2 no prot
	       cut b1 and b2 no prot
	       cut all except b1 off of lca, b1
	       cut all except b2 off of lca, b2
	    */
	    //8.4 : Cut a1 and a2
	    {
#ifdef DEBUG
	    cout << "Case 8.4a cut a1 & a2" << endl;
#endif
#ifdef DEBUG_CASE_COUNTER
	    case_counter.case_84++;
#endif

	      vector<Node*> to_cut = {T2_a1, T2_a2};
	      vector<Node*> to_cut_except = {};
	      MULT_BB_CUT_AND_RESOLVE(to_cut, to_cut_except, NULL);
	    }
	    //8.4 : Cut a1 and a2's B's
	    {
#ifdef DEBUG
	    cout << "Case 8.4b cut B1 & B2" << endl;
#endif
	      vector<Node*> to_cut = {};
	      vector<Node*> to_cut_except = {T2_a1, T2_a2};
	      MULT_BB_CUT_AND_RESOLVE(to_cut, to_cut_except, NULL);
	    }

	    //8.4 : All except a1->p, then b1
	    {
#ifdef DEBUG
	    cout << "Case 8.4c cut B1 & B`1" << endl;
#endif
	      vector<Node*> to_cut = {};
	      vector<Node*> to_cut_except = {T2_a1->parent(), T2_a1};
	      MULT_BB_CUT_AND_RESOLVE(to_cut, to_cut_except, T2_a2);
	    }
	    //8.4 : All except a2->p, then b2
	    {
#ifdef DEBUG
	    cout << "Case 8.4d cut B2 & B`2" << endl;
#endif
	      vector<Node*> to_cut = {};
	      vector<Node*> to_cut_except = {T2_a2->parent(), T2_a2};
	      MULT_BB_CUT_AND_RESOLVE(to_cut, to_cut_except, T2_a1);
	    }
	    Node* pl = arbitrary_lca->parent();
	    bool pl_has_non_aj_child = false;
	    if (pl != NULL) {
	      for (auto i = pl->get_children().begin(); i != pl->get_children().end(); i++) {
		if (*i != arbitrary_lca && descendant_count[(*i)->get_preorder_number()] != -1) {
		  pl_has_non_aj_child = true;
		  break;
		}
	      }
	    }
	    bool gpl_has_non_aj_child = false;
	    //place this conditional so we do not unnecessarily iterate
	    //If pl_has_non_aj_child then the if will evaluate to true later on anyways
	    if (!pl_has_non_aj_child) {
	      if (pl != NULL) {
		Node* gpl = pl->parent();		    		    
		if (gpl != NULL) {
		  for (auto i = gpl->get_children().begin(); i != gpl->get_children().end(); i++) {
		    if (*i != pl && descendant_count[(*i)->get_preorder_number()] != -1) {
		      gpl_has_non_aj_child = true;
		      break;
		    }
		  }
		}
	      }
	    }
	    //8.4 special case
	    if (arbitrary_lca->parent() == NULL ||
		pl_has_non_aj_child ||
		gpl_has_non_aj_child) {
	      //8.4 Cut a1 prot a2
	      {
#ifdef DEBUG
	      cout << "Case 8.4e cut a1" << endl;
#endif
		vector<Node*> to_cut = {T2_a1};
		vector<Node*> to_cut_except = {};
		MULT_BB_CUT_AND_RESOLVE(to_cut, to_cut_except, T2_a2);	
	      }
	      //8.4 Cut a2 prot a1
	      {
#ifdef DEBUG
	      cout << "Case 8.4f a2" << endl;
#endif
		vector<Node*> to_cut = {T2_a2};
		vector<Node*> to_cut_except = {};
		MULT_BB_CUT_AND_RESOLVE(to_cut, to_cut_except, T2_a1);	
	      }

	    }
	  }

	  //step 8.5
	  else if (T1_sibling_group->get_children().size() > 2 &&
		   deepest_siblings.size() > 2 &&
		   all_s1) {
	    /*
	      recurse on cut all a1 through ar no prot,
	      cut all b1 through br no prot,
	      for each ai, cut all a's except ai prot ai
	      for each ai, cut all B's except for ai's B prot ai
	    */
	    //8.5 all ai
	    {
#ifdef DEBUG
	    cout << "Case 8.5a cut all a1-ar" << endl;
#endif
#ifdef DEBUG_CASE_COUNTER
	    case_counter.case_85++;
#endif	    
	      vector<Node*> to_cut = {};
	      vector<Node*> to_cut_except = {};
	      for (int i = 0; i < deepest_siblings.size(); i++) {
		to_cut.push_back(deepest_siblings[i]);
	      }
	      MULT_BB_CUT_AND_RESOLVE(to_cut, to_cut_except, NULL);
	    }
	    //8.5 cut all bi
	    {
#ifdef DEBUG
	      cout << "Case 8.5b cut all B1-Br" << endl;
#endif
	      vector<Node*> to_cut = {};
	      vector<Node*> to_cut_except = {};
	      for (int i = 0; i < deepest_siblings.size(); i++) {
		to_cut_except.push_back(deepest_siblings[i]);
	      }
	      MULT_BB_CUT_AND_RESOLVE(to_cut, to_cut_except, NULL);
	    }
	    //8.5 for each ai cut all a's except ai
	    {	      
	      for (int i = 0; i < deepest_siblings.size(); i++) {
#ifdef DEBUG
		cout << "Case 8.5c cut all a1-ar except ai" << endl;
#endif
		vector<Node*> to_cut = {};
		vector<Node*> to_cut_except = {};
		for (int j = 0; j < deepest_siblings.size(); j++) {
		  if (i != j) {
		    to_cut.push_back(deepest_siblings[j]);
		  }
		}
		MULT_BB_CUT_AND_RESOLVE(to_cut, to_cut_except, deepest_siblings[i]);//TODO protect ai
	      }
	    }
	    //8.5 for each ai cut all b's except ai's
	    {	      
	      for (int i = 0; i < deepest_siblings.size(); i++) {
#ifdef DEBUG
		cout << "Case 8.5d cut all B1-Br except Bi" << endl;
#endif
		vector<Node*> to_cut = {};
		vector<Node*> to_cut_except = {};
		for (int j = 0; j < deepest_siblings.size(); j++) {
		  if (i != j) {
		    to_cut_except.push_back(deepest_siblings[j]);
		  }
		}
		MULT_BB_CUT_AND_RESOLVE(to_cut, to_cut_except, deepest_siblings[i]);//TODO protect ai
	      }
	    }
	  }

	  //step 8.6
	  else if (T1_sibling_group->get_children().size() > 2 &&
		   deepest_siblings.size() == 2 &&
		   s_map[T2_a1] >= 2 &&
		   T2_a2->parent() == arbitrary_lca && //equivalent to s_map[T2_a2] == 0
		   (
		    arbitrary_lca->parent() == NULL ||
		    lca_p_contains_sibling
		    )) {
	    /*
	      recurse on cut a1 no prot,
	      cut all B's leading up to LCA from a1, no prot,
	      cut B2 (essentially a1's whole branch) no prot
	    */
	    //if (T2_a1->parent() != NULL) {
	    //8.6 cut a1
	    {
#ifdef DEBUG	     
	      cout << "Case 8.6a cut a1" << endl;
#endif
#ifdef DEBUG_CASE_COUNTER
	    case_counter.case_86++;
#endif
	      vector<Node*> to_cut = {T2_a1};
	      vector<Node*> to_cut_except = {};
	      MULT_BB_CUT_AND_RESOLVE(to_cut, to_cut_except, NULL);
	    }
	    //8.6 cut all B1s
	    {
#ifdef DEBUG	     
	      cout << "Case 8.6b cut all B1's" << endl;
#endif
	      vector<Node*> to_cut = {};
	      vector<Node*> to_cut_except = {};
	      Node* stepper = T2_a1;
		
	      while(stepper->parent() != arbitrary_lca) { 
		to_cut_except.push_back(stepper);
		stepper = stepper->parent();
	      }
	      //TODO: use list to push_front or figure out how to add from top to bottom 
	      reverse(to_cut_except.begin(), to_cut_except.end());

	      MULT_BB_CUT_AND_RESOLVE(to_cut, to_cut_except, NULL);
	    }
	    //8.6 cut B2
	    {
#ifdef DEBUG	     
	      cout << "Case 8.6c cut B2" << endl;
#endif
	      vector<Node*> to_cut = {};
	      vector<Node*> to_cut_except = {T2_a2};
	      MULT_BB_CUT_AND_RESOLVE(to_cut, to_cut_except, NULL);
	    }	      
	  }

	  //step 8.7
	  else if (T1_sibling_group->get_children().size() > 2 &&
		   s_map[T2_a1] >= 2) {
	    /*
	      recurse on cut a1 no prot,
	      cut a2 prot a1,
	      cut all B1's leading up to LCA no prot,
			   
	    */
	    //8.7 cut a1
	    {
#ifdef DEBUG	     
	      cout << "Case 8.7a cut a1" << endl;
#endif
#ifdef DEBUG_CASE_COUNTER
	    case_counter.case_87++;
#endif
	      vector<Node*> to_cut = {T2_a1};
	      vector<Node*> to_cut_except = {};
	      MULT_BB_CUT_AND_RESOLVE(to_cut, to_cut_except, NULL);
	    }	      
	    //8.7 cut a2
	    {
#ifdef DEBUG	     
	      cout << "Case 8.7b cut a2" << endl;
#endif
	      vector<Node*> to_cut = {T2_a2};
	      vector<Node*> to_cut_except = {};
	      MULT_BB_CUT_AND_RESOLVE(to_cut, to_cut_except, T2_a1);
	    }
	    //8.7 cut all B1s
	    {
#ifdef DEBUG	     
	      cout << "Case 8.7c cut all B1's" << endl;
#endif
	      vector<Node*> to_cut = {};
	      vector<Node*> to_cut_except = {};
	      Node* stepper = T2_a1;
		
	      while(stepper->parent() != arbitrary_lca) { 
		to_cut_except.push_back(stepper);
		stepper = stepper->parent();
	      }
	      //TODO: use list to push_front or figure out how to add from top to bottom 
	      reverse(to_cut_except.begin(), to_cut_except.end());

	      MULT_BB_CUT_AND_RESOLVE(to_cut, to_cut_except, NULL);		
	    }
	    //8.7 cut all B2s
	    {
#ifdef DEBUG
	    cout << "Case 8.7d";
	    if (deepest_siblings.size() == 2) {
	      cout << "2" << endl;
	    }
	    else {
	      cout << "n" << endl;
	    }
#endif

	      vector<Node*> to_cut = {};
	      vector<Node*> to_cut_except = {};
	      Node* stepper = T2_a2;
		
	      while(stepper->parent() != arbitrary_lca) { 
		to_cut_except.push_back(stepper);
		stepper = stepper->parent();
	      }
	      //TODO: use list to push_front or figure out how to add from top to bottom 
	      reverse(to_cut_except.begin(), to_cut_except.end());
	      //8.7 cut all B2's, cut B`2, otherwise r > 2
	      if(deepest_siblings.size() == 2) {
		to_cut_except.push_back(T2_a2); // cut B`2 at the end
	      }
	      //This exception should be caught by 8.6 case.
	      /*
		if (to_cut_except.size() == 0) {
		cout << "\n\n\nto_cut_except empty in 8.7n. Parent is lca ERROR \n\n\n"; 
		}*/

	      MULT_BB_CUT_AND_RESOLVE(to_cut, to_cut_except, T2_a1);	
	    }
	  }

	}

	}
	else { //4_MULT_BRANCH

	  //To debug: set all to true. Make all cuts all the time
	  //Otherwise uses same logic as approximation algorithm
	  bool cut_a1 = true;
	  bool cut_b1 = true;
	  bool cut_a2 = true;
	  bool cut_b2 = true;

	  if (T1_sibling_group->get_children().size() == 2) {
	    /*
	      7.1 case
	      Cut a1, pa1, a2 in F2, add 3 to num_cut
	    */	
	    cut_a1   = true;
	    cut_b1 = true;
	    cut_a2   = true;
	    /*
	      7.2 case
	      Cut a1, a2, pa1, pa2, add 4 to num_cuts
	    */
	    if (previous_group == T1_sibling_group) {
	      cut_b2 = true;
	    }
	  } // size == 2
      
	  else if (T1_sibling_group->get_children().size() > 2) {
	    /*
	      7.3 case
	      If a2's parent's only sibling is part of the sibling group, 
	      and a1's parent is a root or has a sibling that is not part of the sibling group
	      then cut a2 and a2_p otherwise a1 and a1_p
	    */
	    if (previous_group != T1_sibling_group) {
	      Node* T2_a2_p = T2_a2->parent();
	      bool x_2 = false;
	      bool a2_p_one_sibling = (T2_a2_p != NULL) &&
		(T2_a2_p->parent() != NULL) &&
		(T2_a2_p->parent()->get_children().size() == 2);	  
	      if (a2_p_one_sibling) {
		list<Node *> group = T1_sibling_group->get_children();
		//get the other one
		Node *a2_p_sibling = T2_a2_p->parent()->get_children().front() == T2_a2_p ?
		  T2_a2_p->parent()->get_children().back() :
		  T2_a2_p->parent()->get_children().front();	
		//check if it is part of sibling group
		bool a2_p_sibling_in_group = descendant_count[a2_p_sibling->get_preorder_number()] == -1;
		if (a2_p_sibling_in_group) {
		  Node* T2_a1_p = T2_a1->parent();
		  bool a1_p_is_root = T2_a1_p->parent() == NULL;
		  bool a1_p_sibling_not_in_group = false;
		  if (!a1_p_is_root) {
		    list<Node*> a1_p_siblings = T2_a1_p->parent()->get_children();
		    for (list<Node*>::iterator i = a1_p_siblings.begin(); i != a1_p_siblings.end(); i++) {
		      if (descendant_count[(*i)->get_preorder_number()] != -1) {
			a1_p_sibling_not_in_group = true;
			break;
		      }
		    }
		  }
		  if (a2_p_one_sibling && a2_p_sibling_in_group && (a1_p_is_root || a1_p_sibling_not_in_group)) {
		    x_2 = true;
		  }
		}
	      }
	      if (x_2){
		cut_a2 = true;
		cut_b2 = true;
		//cut_a1 = true;
		//cut_b1 = true;
	      }
	      else {
		cut_a1 = true;
		cut_b1 = true;
		cut_a2 = true;
		//cut_b2 = true;
	      }
	    }
	    /*7.4 case
	      cut a1 and a1_p
	    */
	    else if (previous_group == T1_sibling_group) {
	      cut_a1   = true;
	      cut_b1 = true;
	    }
	  }
	  /*
	    Cutting section
	  */
	  Node *T2_a1_p = T2_a1->parent();

	  if (cut_a1) {	  
	    if (T2_a1_p != NULL) {
#ifdef DEBUG
	      cout << "Case Cut a1" << endl;
#endif
	      vector<Node*> to_cut = {T2_a1};
	      vector<Node*> to_cut_except = {};
	      MULT_BB_CUT_AND_RESOLVE(to_cut, to_cut_except, NULL);
	    }
	    //mimics case 8.1, should add rho
	    else if (T2_a1 == T2->get_component(0) && arbitrary_lca == NULL) {
	      if (!T1->contains_rho()) {
		MULT_RHO_CUT_AND_RESOLVE(T2_a1, NULL);
	      }
	    }
	  }
	  if (cut_b1) {
	    if (T2_a1_p != NULL) {
#ifdef DEBUG
	      cout << "Case Cut b1" << endl;
#endif
	      vector<Node*> to_cut = {};
	      vector<Node*> to_cut_except = {T2_a1};
	      MULT_BB_CUT_AND_RESOLVE(to_cut, to_cut_except, NULL);
	    }
	  }
	  Node *T2_a2_p = T2_a2->parent();
	  if (cut_a2){
#ifdef DEBUG
	      cout << "Case Cut a2" << endl;
#endif
	    if (T2_a2_p != NULL) {
	      vector<Node*> to_cut = {T2_a2};
	      vector<Node*> to_cut_except = {};
	      MULT_BB_CUT_AND_RESOLVE(to_cut, to_cut_except, NULL);
	    }
	    else if (T2_a2 == T2->get_component(0) && arbitrary_lca == NULL) {
	      
	      if (!T1->contains_rho()) {
		MULT_RHO_CUT_AND_RESOLVE(T2_a2, NULL);
	      }
	    }

	  }
	  if (cut_b2) {
	    if (T2_a2_p != NULL) {
#ifdef DEBUG
	      cout << "Case Cut b2" << endl;
#endif
	      vector<Node*> to_cut = {};
	      vector<Node*> to_cut_except = {T2_a2};
	      MULT_BB_CUT_AND_RESOLVE(to_cut, to_cut_except, NULL);
	    }
	  }
	}
	sibling_groups->pop_back();
	return best_k;
      }//else (cutting)
      previous_group = T1_sibling_group;
    }//!sibling_groups->empty()

  } //while(!sibling_groups->empty() && !singletons->empty()

  //Made it to end, so add to results
  if (k >= 0) {
    //if (PREFER_RHO && !AFs->empty() && !AFs->front().first.contains_rho() && T1->contains_rho()) {
    if (true && !AFs->empty() && !AFs->front().first.contains_rho() && T1->contains_rho()) {
      if (!ALL_MAFS)
	AFs->clear();
      AFs->push_front(make_pair(Forest(T1),Forest(T2)));
      *num_ties = 2;
    }
    else if (ALL_MAFS || AFs->empty()) {
      AFs->push_back(make_pair(Forest(T1),Forest(T2)));
    }
    //else if (!PREFER_RHO || AFs->front().first.contains_rho() == T1->contains_rho()) {
    else if (false || AFs->front().first.contains_rho() == T1->contains_rho()) {
      if (rand() < RAND_MAX/ *num_ties) {
	AFs->clear();
	AFs->push_back(make_pair(Forest(T1),Forest(T2)));
      }
      (*num_ties)++;
    }
  }
  return k;
}


/* rSPR_3_approx
 * Calculate an approximate maximum agreement forest and SPR distance
 * RETURN At most 3 times the rSPR distance
 * NOTE: destructive. The computed forests replace T1 and T2.
 */
int rSPR_3_approx(Forest *T1, Forest *T2) {
	// find sibling pairs of T1
	// match up nodes of T1 and T2
	if (!sync_twins(T1, T2))
return 0;
	// find singletons of T2
	list<Node *> *sibling_pairs = T1->find_sibling_pairs();
	list<Node *> singletons = T2->find_singletons();
	int ans = rSPR_3_approx_hlpr(T1, T2, &singletons, sibling_pairs);
	delete sibling_pairs;
	return ans;
}

// rSPR_3_approx recursive helper function
int rSPR_3_approx_hlpr(Forest *T1, Forest *T2, list<Node *> *singletons,
list<Node *> *sibling_pairs) {
	//return rSprAlgorithm3Approx::rSPR_3_approx_hlpr_Inline(APPROX_CHECK_COMPONENT,T1, T2, singletons,sibling_pairs);


	int num_cut = 0;
	while(!singletons->empty() || !sibling_pairs->empty()) {
// Case 1 - Remove singletons
while(!singletons->empty()) {


	Node *T2_a = singletons->back();
	singletons->pop_back();
  // find twin in T1
	Node *T1_a = T2_a->get_twin();
	// if this is in the first component of T_2 then
	// it is not really a singleton.
	if (T2_a == T2->get_component(0))
		continue;

	Node *T1_a_parent = T1_a->parent();
	if (T1_a_parent == NULL)
		continue;
	bool potential_new_sibling_pair = T1_a_parent->is_sibling_pair();
	// cut the edge above T1_a
	T1_a->cut_parent();
	T1->add_component(T1_a);
	if (T1_a->get_sibling_pair_status() > 0)
		T1_a->clear_sibling_pair(sibling_pairs);
	//delete(T1_a);

	Node *node = T1_a_parent->contract();
	if (node != NULL && potential_new_sibling_pair && node->is_sibling_pair()){
		node->rchild()->add_to_front_sibling_pairs(sibling_pairs, 2);
		node->lchild()->add_to_front_sibling_pairs(sibling_pairs, 1);
	}

}
if(!sibling_pairs->empty()) {
	Node *T1_a = sibling_pairs->back();
	sibling_pairs->pop_back();
	Node *T1_c = sibling_pairs->back();
	sibling_pairs->pop_back();
	T1_a->clear_sibling_pair_status();
	T1_c->clear_sibling_pair_status();
	if (T1_a->parent() == NULL || T1_a->parent() != T1_c->parent()) {
		continue;
	}
	Node *T1_ac = T1_a->parent();
	// lookup in T2 and determine the case
	Node *T2_a = T1_a->get_twin();
	Node *T2_c = T1_c->get_twin();

	// Case 2 - Contract identical sibling pair
	if (T2_a->parent() != NULL && T2_a->parent() == T2_c->parent()) {
		Node *T2_ac = T2_a->parent();
		T1_ac->contract_sibling_pair();
		T2_ac->contract_sibling_pair();
		T1_ac->set_twin(T2_ac);
		T2_ac->set_twin(T1_ac);
		T1->add_deleted_node(T1_a);
		T1->add_deleted_node(T1_c);
		T2->add_deleted_node(T2_a);
		T2->add_deleted_node(T2_c);

		// check if T2_ac is a singleton
		if (T2_ac->is_singleton() && !T1_ac->is_singleton() && T2_ac != T2->get_component(0))
			singletons->push_back(T2_ac);
		// check if T1_ac is part of a sibling pair
		if (T1_ac->parent() != NULL && T1_ac->parent()->is_sibling_pair()) {
			T1_ac->parent()->lchild()->add_to_sibling_pairs(sibling_pairs, 1);
			T1_ac->parent()->rchild()->add_to_sibling_pairs(sibling_pairs, 2);
		}
	}
	// Case 3
	else {

		//  ensure T2_a is below T2_c
		if (T2_a->get_depth() < T2_c->get_depth()) {
			swap(&T1_a, &T1_c);
			swap(&T2_a, &T2_c);
		}
		else if (T2_a->get_depth() == T2_c->get_depth()) {
			if (T2_a->parent() && T2_c->parent() &&
					(T2_a->parent()->get_depth() <
					T2_c->parent()->get_depth())) {
			swap(&T1_a, &T1_c);
			swap(&T2_a, &T2_c);
			}
		}

		// get T2_b
		Node *T2_ab = T2_a->parent();
		Node *T2_b = T2_ab->rchild();
		if (T2_b == T2_a)
			T2_b = T2_ab->lchild();
		// cut T1_a, T1_c, T2_a, T2_b, T2_c

		bool cut_b_only = false;
		if (T2_a->parent() != NULL && T2_a->parent()->parent() != NULL && T2_a->parent()->parent() == T2_c->parent()) {
			cut_b_only = true;
			T1_a->add_to_sibling_pairs(sibling_pairs,1);
			T1_c->add_to_sibling_pairs(sibling_pairs,2);
		}

		if (!cut_b_only) {
			T1_a->cut_parent();
			T1_c->cut_parent();
			// contract parents
			Node *node = T1_ac->contract();
			// check for T1_ac sibling pair
			if (node != NULL && node && node->is_sibling_pair()){
				node->lchild()->add_to_sibling_pairs(sibling_pairs,1);
				node->rchild()->add_to_sibling_pairs(sibling_pairs,2);
			}
		}

		bool same_component = true;
		if (APPROX_CHECK_COMPONENT)
			same_component = (T2_a->find_root() == T2_c->find_root());

		if (!cut_b_only) {
			T2_a->cut_parent();
			num_cut++;
		}
		bool cut_b = false;
		if (same_component && T2_ab->parent() != NULL) {
			T2_b->cut_parent();
			num_cut++;
			cut_b = true;
		}
		// T2_b will move up after contraction
		else {
			T2_b = T2_b->parent();
		}
		// check for T2 parents as singletons
		Node *node = T2_ab->contract();
		if (node != NULL && node->is_singleton()
				&& node != T2->get_component(0))
			singletons->push_back(node);

		// if T2_c is gone then its replacement is in singleton list
		// contract might delete old T2_c, see where it is
		bool add_T2_c = true;
		T2_c = T1_c->get_twin();
		// ignore T2_c if it is a singleton
		if (T2_c != node && T2_c->parent() != NULL && !cut_b_only) {

			Node *T2_c_parent = T2_c->parent();
			T2_c->cut_parent();
			num_cut++;
			node = T2_c_parent->contract();
			if (node != NULL && node->is_singleton()
					&& node != T2->get_component(0))
				singletons->push_back(node);
		}
		else {
			add_T2_c = false;
		}


		if (!cut_b_only)
			T1->add_component(T1_a);
		if (!cut_b_only)
			T1->add_component(T1_c);
		// put T2 cut parts into T2
		if (!cut_b_only) {
			T2->add_component(T2_a);
		}
		// may have already been added
		if (cut_b) {
			T2->add_component(T2_b);
		}
		// problem if c is deleted
		if (add_T2_c) {
			T2->add_component(T2_c);
		}

		// may have already been added
		if (T2_b->is_leaf())
			singletons->push_back(T2_b);

	}
}
	}
// if the first component of the forests differ then we have to cut p
if (T1->get_component(0)->get_twin() != T2->get_component(0)) {
	num_cut++;
	T1->add_rho();
	T2->add_rho();
}
return num_cut;
}
/*******************************************************************************
	RSPR WORSE_3_APPROX
*******************************************************************************/

/* rSPR_worse_3_approx
 * Calculate an approximate maximum agreement forest and SPR distance
 * RETURN At most 3 times the rSPR distance
 * NOTE: destructive. The computed forests replace T1 and T2.
 * T1 must be a binary tree. T2 can be a multifurcating forest.
 */
int rSPR_worse_3_approx(Forest *T1, Forest *T2) {
	return rSPR_worse_3_approx(T1, T2, true);
}

int rSPR_worse_3_approx(Forest *T1, Forest *T2, bool sync) {
	// match up nodes of T1 and T2
	if (sync) {
if (!sync_twins(T1, T2))
	return 0;
	}
//	cout << "T1: "; T1->print_components();
//	cout << "T2: "; T2->print_components();
	// find sibling pairs of T1
	list<Node *> *sibling_pairs = T1->find_sibling_pairs();
	// find singletons of T2
	list<Node *> singletons = T2->find_singletons();
	list<pair<Forest,Forest> > AFs = list<pair<Forest,Forest> >();

	Forest *F1;
	Forest *F2;

	int ans = rSPR_worse_3_approx_hlpr(T1, T2, &singletons, sibling_pairs, &F1, &F2, true);

	F1->swap(T1);
	F2->swap(T2);
	sync_twins(T1,T2);


	delete sibling_pairs;
	delete F1;
	delete F2;
	return ans;
}

int rSPR_worse_3_approx_distance_only(Forest *T1, Forest *T2) {
if (!sync_twins(T1, T2))
	return 0;
	list<Node *> *sibling_pairs = T1->find_sibling_pairs();
	list<Node *> singletons = T2->find_singletons();
	list<pair<Forest,Forest> > AFs = list<pair<Forest,Forest> >();

	int ans = rSPR_worse_3_approx_hlpr(T1, T2, &singletons, sibling_pairs, NULL, NULL, false);

	delete sibling_pairs;
	return ans;
}

int rSPR_worse_3_approx(Node *subtree, Forest *T1, Forest *T2) {
	return rSPR_worse_3_approx(subtree, T1, T2, true);
}

int rSPR_worse_3_approx(Node *subtree, Forest *T1, Forest *T2, bool sync) {
	// match up nodes of T1 and T2
	if (sync) {
if (!sync_twins(T1, T2))
	return 0;
	}
//	cout << "T1: "; T1->print_components();
//	cout << "T2: "; T2->print_components();
	// find sibling pairs of T1
	list<Node *> *sibling_pairs = subtree->find_sibling_pairs();
	// find singletons of T2
	list<Node *> singletons = T2->find_singletons();
	list<pair<Forest,Forest> > AFs = list<pair<Forest,Forest> >();

	Forest *F1;
	Forest *F2;

	int ans = rSPR_worse_3_approx_hlpr(T1, T2, &singletons, sibling_pairs, &F1, &F2, true);

	F1->swap(T1);
	F2->swap(T2);
	sync_twins(T1,T2);


	delete sibling_pairs;
	delete F1;
	delete F2;
	return ans;
}

// rSPR_worse_3_approx recursive helper function
int rSPR_worse_3_approx_hlpr(Forest *T1, Forest *T2, list<Node *> *singletons, list<Node *> *sibling_pairs, Forest **F1, Forest **F2, bool save_forests) {
	#ifdef DEBUG_APPROX
cout << "rSPR_worse_3_approx_hlpr" << endl;
			cout << "\tT1: ";
			T1->print_components_with_twins();
			cout << "\tT2: ";
			T2->print_components_with_twins();
			cout << "sibling pairs:";
			for (list<Node *>::iterator i = sibling_pairs->begin(); i != sibling_pairs->end(); i++) {
				cout << "  ";
				(*i)->print_subtree_hlpr();
			}
			cout << endl;
	#endif
	int num_cut = 0;
	UndoMachine um = UndoMachine();
	while(!singletons->empty() || !sibling_pairs->empty()) {
// Case 1 - Remove singletons
while(!singletons->empty()) {
	#ifdef DEBUG_APPROX
		cout << "Case 1" << endl;
	#endif

	Node *T2_a = singletons->back();
	singletons->pop_back();
	// find twin in T1
	Node *T1_a = T2_a->get_twin();
	// if this is in the first component of T_2 then
	// it is not really a singleton.
	// TODO: problem when we cluster and have a singleton as the
	//		first comp of T2
	//    NEED TO MODIFY CUTTING?
	// 		HERE AND IN BB?
	if (T2_a == T2->get_component(0))
		continue;

	Node *T1_a_parent = T1_a->parent();
	if (T1_a_parent == NULL)
		continue;
	bool potential_new_sibling_pair = T1_a_parent->is_sibling_pair();
	// cut the edge above T1_a
	um.add_event(new CutParent(T1_a));
	T1_a->cut_parent();
	um.add_event(new AddComponent(T1));
	T1->add_component(T1_a);
	//if (T1_a->get_sibling_pair_status() > 0)
	//	T1_a->clear_sibling_pair(sibling_pairs);
	//delete(T1_a);

	ContractEvent(&um, T1_a_parent);
	Node *node = T1_a_parent->contract();
	if (node != NULL && potential_new_sibling_pair &&
			node->is_sibling_pair()){
		um.add_event(new AddToFrontSiblingPairs(sibling_pairs));
		sibling_pairs->push_front(node->rchild());
		sibling_pairs->push_front(node->lchild());
	}

	#ifdef DEBUG_APPROX
			cout << "\tT1: ";
			T1->print_components();
			cout << "\tT2: ";
			T2->print_components();
	#endif
}
if(!sibling_pairs->empty()) {
	/*
	if (PREORDER_SIBLING_PAIRS) {
		T1->get_component(0)->preorder_number();
		list<Node *>::iterator c;
		list<Node *>::iterator best_sib = sibling_pairs->end();
		int best_prenum = INT_MAX;
		list<Node *>::iterator T1_a_i;
		list<Node *>::iterator T1_c_i;
		for(c = sibling_pairs->begin(); c != sibling_pairs->end(); ) {
				T1_c_i = c;
				T1_c = *c;
				c++;
				T1_a_i = c;
				T1_a = *c;
				c++;
				cout << T1_a->str_subtree() << endl;
				cout << T1_c->str_subtree() << endl;
//					if (T1_a->parent() == NULL || T1_a->parent() != T1_c->parent()) {
//						cout << "invalid" << endl;
//						sibling_pairs->erase(T1_c_i);
//						sibling_pairs->erase(T1_a_i);
//						um.add_event(new PopSiblingPair(T1_a, T1_c, sibling_pairs));
//						continue;
//					}
//					else {
				//int prenum = T1_a->parent()->get_preorder_number();
				int prenum = T1_a->get_preorder_number();
				cout << "prenum=" << prenum << endl;
				cout << "old_prenum=" << best_prenum << endl;
				if (prenum < best_prenum) {
					best_sib = T1_c_i; 
					best_prenum = prenum;
				}
				cout << "new_prenum=" << best_prenum << endl;
//					}
		}
		cout << endl;
		if (best_prenum == INT_MAX)
			continue;
		else {
			T1_c_i = best_sib;
			T1_c = *T1_c_i;
			best_sib++;
			T1_a_i = best_sib;
			T1_a = *T1_a_i;
			sibling_pairs->erase(T1_a_i);
			sibling_pairs->erase(T1_c_i);
		}
	}
	else {
	*/
	Node *T1_a = sibling_pairs->back();
	sibling_pairs->pop_back();
	Node *T1_c = sibling_pairs->back();
	sibling_pairs->pop_back();
	um.add_event(new PopSiblingPair(T1_a, T1_c, sibling_pairs));

	//if (T1_a->get_sibling_pair_status() == 0 ||
	//		T1_c->get_sibling_pair_status() == 0) {
	//	continue;
//			}

	//T1_a->clear_sibling_pair_status();
	//T1_c->clear_sibling_pair_status();
	if (T1_a->parent() == NULL || T1_c->parent() == NULL || T1_a->parent() != T1_c->parent()) {
		continue;
	}
	if (!T1_a->can_be_sibling() || !T1_c->can_be_sibling()
	|| num_cut >= INT_MAX - 3) {
		continue;
	}
	Node *T1_ac = T1_a->parent();
	// lookup in T2 and determine the case
	Node *T2_a = T1_a->get_twin();
	Node *T2_c = T1_c->get_twin();

	#ifdef DEBUG_APPROX
		cout << "Fetching sibling pair" << endl;
		T1_ac->print_subtree();
		cout << "T2_a" << ": ";
		cout << " d=" << T2_a->get_depth() << " ";
		T2_a->print_subtree();
		cout << "T1_c" << ": ";
		T1_c->print_subtree();
		cout << "T2_c" << ": ";
		cout << " d=" << T2_c->get_depth() << " ";
		T2_c->print_subtree();
	#endif

	// Case 2 - Contract identical sibling pair
	if (T2_a->parent() != NULL && T2_a->parent() == T2_c->parent()) {
		#ifdef DEBUG_APPROX
			cout << "Case 2" << endl;
			T1->print_components();
			T2->print_components();
		#endif
		Node *T2_ac = T2_a->parent();
		um.add_event(new ContractSiblingPair(T1_ac));
		T1_ac->contract_sibling_pair_undoable();
		um.add_event(new ContractSiblingPair(T2_ac, T2_a, T2_c, &um));
		Node *T2_ac_new = T2_ac->contract_sibling_pair_undoable(T2_a, T2_c);
		if (T2_ac_new != NULL && T2_ac_new != T2_ac) {
			T2_ac = T2_ac_new;
			um.add_event(new CreateNode(T2_ac));
			um.add_event(new ContractSiblingPair(T2_ac));
			T2_ac->contract_sibling_pair_undoable();
		}
		um.add_event(new SetTwin(T1_ac));
		um.add_event(new SetTwin(T2_ac));
		T1_ac->set_twin(T2_ac);
		T2_ac->set_twin(T1_ac);
		//T2_ac->fix_contracted_order();
		//T1->add_deleted_node(T1_a);
		//T1->add_deleted_node(T1_c);
		//T2->add_deleted_node(T2_a);
		//T2->add_deleted_node(T2_c);

		// check if T2_ac is a singleton
		//if (T2_ac->is_singleton() && !T1_ac->is_singleton() && T2_ac != T2->get_component(0))
		if (T2_ac->is_singleton() && T1_ac != T1->get_component(0) && T2_ac != T2->get_component(0))
			singletons->push_back(T2_ac);
		// check if T1_ac is part of a sibling pair
		if (T1_ac->parent() != NULL && T1_ac->parent()->is_sibling_pair()) {
			um.add_event(new AddToSiblingPairs(sibling_pairs));
			sibling_pairs->push_back(T1_ac->parent()->lchild());
			sibling_pairs->push_back(T1_ac->parent()->rchild());
		}
	}
	// Case 3
	else {
		#ifdef DEBUG_APPROX
			cout << "Case 3" << endl;
		#endif
		
		//  ensure T2_a is below T2_c
		if ((T2_a->get_depth() < T2_c->get_depth()
				&& T2_c->parent() != NULL)
				|| T2_a->parent() == NULL) {
			#ifdef DEBUG_APPROX
				cout << "swapping" << endl;
			#endif
		
			swap(&T1_a, &T1_c);
			swap(&T2_a, &T2_c);

		}
		else if (T2_a->get_depth() == T2_c->get_depth()) {
			if (T2_a->parent() && T2_c->parent() &&
					(T2_a->parent()->get_depth() <
					T2_c->parent()->get_depth())) {
			swap(&T1_a, &T1_c);
			swap(&T2_a, &T2_c);
			}
		}

		// get T2_b
		bool multi_node = false;
		Node *T2_ab = T2_a->parent();
		Node *T2_b = T2_ab;
		if (T2_ab->get_children().size() > 2) {
			multi_node = true;
		}
		else {
			T2_b = T2_ab->rchild();
			if (T2_b == T2_a)
				T2_b = T2_ab->lchild();
		}

		#ifdef DEBUG_APPROX
		cout << "T2_b" << ": ";
		cout.flush();
		T2_b->print_subtree();
	#endif
		// cut T1_a, T1_c, T2_a, T2_b, T2_c

		bool cut_a_only = false;
		bool cut_b_only = false;
		bool cut_c_only = false;
		bool cut_b_only_if_not_a_or_c = false;
		if (APPROX_CUT_ONE_B && T2_a->parent() != NULL && T2_a->parent()->parent() != NULL && T2_a->parent()->parent() == T2_c->parent() && !multi_node
						&& (!APPROX_EDGE_PROTECTION || !T2_b->is_protected())) {
			cut_b_only = true;
			um.add_event(new AddToSiblingPairs(sibling_pairs));
			sibling_pairs->push_back(T1_c);
			sibling_pairs->push_back(T1_a);
		}
	if (APPROX_CUT_TWO_B && !cut_b_only && T1_ac->parent() != NULL
						&& (!APPROX_EDGE_PROTECTION || !T2_b->is_protected())) {
		Node *T1_s = T1_ac->get_sibling();
		if (T1_s->is_leaf()) {
			Node *T2_l = T2_a->parent()->parent();
			if (T2_l != NULL) {
				if (T2_c->parent() != NULL && T2_c->parent()->parent() == T2_l
						&& T2_a->parent()->get_children().size() > 2
						&& T2_c->parent()->get_children().size() > 2) {
					if (T2_l->get_sibling() == T1_s->get_twin()) {
						cut_b_only=true;
					}
					else if (T2_l->parent() == NULL &&
							(T2->contains_rho() ||
							 T2->get_component(0) != T2_l)) {
						cut_b_only_if_not_a_or_c=true;
					}
				}
				else if ((T2_l = T2_l->parent()) != NULL
						&& T2_c->parent() == T2_l
						&& T2_a->parent()->get_children().size() > 2
						&& T2_a->parent()->parent()->get_children().size() > 2) {
					if (T2_l->get_sibling() == T1_s->get_twin()) {
						cut_b_only=true;
					}
					else if (T2_l->parent() == NULL &&
							(T2->contains_rho() ||
							 T2->get_component(0) != T2_l)) {
						cut_b_only_if_not_a_or_c=true;
					}
				}
			}
		}
	}
	if (APPROX_REVERSE_CUT_ONE_B && !cut_b_only && T1_ac->parent() != NULL) {
		Node *T1_s = T1_ac->get_sibling();
		if (T1_s->is_leaf()) {
			if (T1_s->get_twin()->parent() == T2_a->parent()//) {
						&& (!APPROX_EDGE_PROTECTION || !T2_c->is_protected())) {
				cut_c_only=true;
			}
			else if (T1_s->get_twin()->parent() == T2_c->parent()//) {
						&& (!APPROX_EDGE_PROTECTION || !T2_a->is_protected())
							&& T2_c->parent()->get_children().size() <= 2) {
				cut_a_only=true;
			}
		}
		else if (APPROX_REVERSE_CUT_ONE_B_2) {
			if (T2_c->parent() != NULL
				&& rSprUtility::chain_match(T1_s, T2_c->get_sibling(), T2_a) //)
						&& (!APPROX_EDGE_PROTECTION || !T2_a->is_protected()))
			cut_a_only = true;
		}
	}
	if (APPROX_CUT_TWO_B_ROOT && cut_a_only == false && cut_c_only == false
			&& cut_b_only_if_not_a_or_c == true) {
		cut_b_only = true;
	}
	/*
	if (CUT_LOST) {
		if (T1_a->num_lost_children() > 0
				|| T2_a->num_lost_children() > 0) {
			cut_a_only = true;
			cut_b_only = false;
			cut_c_only = false;
			num_cut-=3;
		}
		else if (T1_c->num_lost_children() > 0
				|| T2_c->num_lost_children() > 0) {
			cut_a_only = false;
			cut_b_only = false;
			cut_c_only = true;
			num_cut-=3;
		}
		else if (T2_b->is_leaf()) {
			if (T2_b->num_lost_children() > 0
				|| T2_b->get_twin()->num_lost_children() > 0) {
			cut_a_only = false;
			cut_b_only = true;
			cut_c_only = false;
			num_cut-=3;
			}
		}
	}
	*/


		Node *node;

		bool cut_a = false;
		bool cut_c = false;
		if (!cut_b_only || T2_a->parent()->get_children().size() > 2) {
			if (!cut_c_only &&
					(!APPROX_EDGE_PROTECTION
					 	|| (!T2_a->is_protected()
							&& (T2_a->parent()->parent() != NULL
								|| !T2_b->is_protected()
								|| T2_a->parent()->get_children().size() > 2)))) {
//					|| cut_a_only)) {
				um.add_event(new CutParent(T1_a));
				T1_a->cut_parent();
				cut_a = true;

				ContractEvent(&um, T1_ac);
				node = T1_ac->contract();
			}
			else
				node = T1_ac;
			if (!cut_a_only &&
					(!APPROX_EDGE_PROTECTION
					 	|| (!T2_c->is_protected()
						&& (T2_c->parent() == NULL
								|| T2_c->parent()->parent() != NULL
								|| !T2_c->get_sibling()->is_protected()
								|| T2_c->parent()->get_children().size() > 2)))) {// &&
//					|| cut_c_only)) {
				um.add_event(new CutParent(T1_c));
				T1_c->cut_parent();
				cut_c = true;

				if (node) {
					ContractEvent(&um, node);
					node = node->contract();
				}
			}

			// contract parents
			// check for T1_ac sibling pair
			if (node && node->is_sibling_pair()){
				um.add_event(new AddToSiblingPairs(sibling_pairs));
				sibling_pairs->push_back(node->lchild());
				sibling_pairs->push_back(node->rchild());
			}
		}

		bool same_component = true;
		if (APPROX_CHECK_COMPONENT && !cut_a_only && !cut_c_only)
			same_component = (T2_a->find_root() == T2_c->find_root());

		Node *T2_ab_parent = T2_ab->parent();
		node = T2_ab;
		if (cut_a) {
			um.add_event(new CutParent(T2_a));
			T2_a->cut_parent();

			//ContractEvent(&um, T2_ab);
			//node = T2_ab->contract();
		}
		bool cut_b = false;
		if (same_component && T2_ab_parent != NULL
				&& !cut_a_only && !cut_c_only
				&& (!APPROX_EDGE_PROTECTION
					|| (!T2_b->is_protected() ))) {
//							&& (T2_b->parentT2_a->parent()->parent() != NULL
//								|| !T2_a->is_protected())))) {
//					|| cut_b_only)) {
			if (multi_node) {
				T2_b = T2_ab;
				um.add_event(new CutParent(T2_ab));
				T2_ab->cut_parent();
				if (T2_a->parent() != NULL) {
					um.add_event(new CutParent(T2_a));
					T2_a->cut_parent();
					um.add_event(new AddChild(T2_a));
					T2_ab_parent->add_child(T2_a);
				}
				else
					node = T2_ab_parent;
			}
			else {
				um.add_event(new CutParent(T2_b));
				T2_b->cut_parent();
				//ContractEvent(&um, node);
				//node = node->contract();
			}
			cut_b = true;
		}
		// T2_b will move up after contraction
		else if (!multi_node) {
			T2_b = T2_b->parent();
		}
		if (node != NULL) {
			ContractEvent(&um, node);
			node = node->contract();
		// check for T2 parents as singletons
		if (node != NULL && node->is_singleton()
				&& node != T2->get_component(0))
			singletons->push_back(node);
		}

		// if T2_c is gone then its replacement is in singleton list
		// contract might delete old T2_c, see where it is
		bool add_T2_c = true;
		T2_c = T1_c->get_twin();
		// ignore T2_c if it is a singleton
		if (cut_c && T2_c != node && T2_c->parent() != NULL) {
			Node *T2_c_parent = T2_c->parent();
			um.add_event(new CutParent(T2_c));
			T2_c->cut_parent();
			ContractEvent(&um, T2_c_parent);
			node = T2_c_parent->contract();
			if (node != NULL && node->is_singleton()
					&& node != T2->get_component(0))
				singletons->push_back(node);
		}
		else {
			add_T2_c = false;
		}

		
		if (cut_a) {
			um.add_event(new AddComponent(T1));
			T1->add_component(T1_a);
			um.add_event(new AddComponent(T2));
			T2->add_component(T2_a);
		}
		if (cut_c) {
			um.add_event(new AddComponent(T1));
			T1->add_component(T1_c);
		}
		if (cut_b) {
			um.add_event(new AddComponent(T2));
			T2->add_component(T2_b);
		}
		// problem if c is deleted
		if (add_T2_c) {
			um.add_event(new AddComponent(T2));
			T2->add_component(T2_c);
		}

		// may have already been added
		if (T2_b->is_leaf() && cut_b)
			singletons->push_back(T2_b);

		num_cut+=3;

		if (cut_a == false && cut_b == false && cut_c == false) {
			num_cut = INT_MAX-3;
		}

	}
}
	}
// if the first component of the forests differ then we have cut p
if (T1->get_component(0)->get_twin() != T2->get_component(0)) {
	if (!T1->contains_rho()) {
		um.add_event(new AddRho(T1));
		um.add_event(new AddRho(T2));
		T1->add_rho();
		T2->add_rho();
	}
	else
		// hack to ignore rho when it shouldn't be in a cluster
		num_cut -=3;
}
if (save_forests) {
	*F1 = new Forest(T1);
	*F2 = new Forest(T2);
}
 
#ifdef DEBUG_APPROX
#ifdef DEBUG_UNDO
 while(um.num_events() > 0) {
		cout << "Undo step " << um.num_events() << endl;
		cout << "T1: ";
		T1->print_components();
		cout << "T2: ";
		T2->print_components();
			cout << "sibling pairs:";
			for (list<Node *>::iterator i = sibling_pairs->begin(); i != sibling_pairs->end(); i++) {
				cout << "  ";
				(*i)->print_subtree_hlpr();
			}
			cout << endl;
	 um.undo();
	cout << endl;
 }
#else
 um.undo_all();
#endif
#else
 um.undo_all();
#endif

 
//		 for(int i = 0; i < T1->num_components(); i++)
//		 	T1->get_component(i)->fix_parents();
//		 for(int i = 0; i < T2->num_components(); i++)
//		 	T2->get_component(i)->fix_parents();
return num_cut;
}

/*******************************************************************************
	RSPR WORSE_3_APPROX_BINARY
*******************************************************************************/

int rSPR_worse_3_approx_binary(Forest *T1, Forest *T2, bool sync) {
	// match up nodes of T1 and T2
	if (sync) {
if (!sync_twins(T1, T2))
	return 0;
	}
	// find sibling pairs of T1
	list<Node *> *sibling_pairs = T1->find_sibling_pairs();
	// find singletons of T2
	list<Node *> singletons = T2->find_singletons();
	list<pair<Forest,Forest> > AFs = list<pair<Forest,Forest> >();

	Forest *F1;
	Forest *F2;

	int ans = rSPR_worse_3_approx_binary_hlpr(T1, T2, &singletons, sibling_pairs, &F1, &F2, true);

	F1->swap(T1);
	F2->swap(T2);
	sync_twins(T1,T2);


	delete sibling_pairs;
	delete F1;
	delete F2;
	return ans;
}

// rSPR_worse_3_approx_binary recursive helper function
int rSPR_worse_3_approx_binary_hlpr(Forest *T1, Forest *T2, list<Node *> *singletons, list<Node *> *sibling_pairs, Forest **F1, Forest **F2, bool save_forests) {
	#ifdef DEBUG_APPROX
cout << "rSPR_worse_3_approx_binary_hlpr" << endl;
			cout << "\tT1: ";
			T1->print_components_with_twins();
			cout << "\tT2: ";
			T2->print_components_with_twins();
			cout << "sibling pairs:";
			for (list<Node *>::iterator i = sibling_pairs->begin(); i != sibling_pairs->end(); i++) {
				cout << "  ";
				(*i)->print_subtree_hlpr();
			}
			cout << endl;
	#endif
	int num_cut = 0;
	UndoMachine um = UndoMachine();
	while(!singletons->empty() || !sibling_pairs->empty()) {
// Case 1 - Remove singletons
while(!singletons->empty()) {
	#ifdef DEBUG_APPROX
		cout << "Case 1" << endl;
	#endif

	Node *T2_a = singletons->back();
	singletons->pop_back();
	// find twin in T1
	Node *T1_a = T2_a->get_twin();
	// if this is in the first component of T_2 then
	// it is not really a singleton.
	// TODO: problem when we cluster and have a singleton as the
	//		first comp of T2
	//    NEED TO MODIFY CUTTING?
	// 		HERE AND IN BB?
	if (T2_a == T2->get_component(0))
		continue;

	Node *T1_a_parent = T1_a->parent();
	if (T1_a_parent == NULL)
		continue;
	bool potential_new_sibling_pair = T1_a_parent->is_sibling_pair();
	// cut the edge above T1_a
	um.add_event(new CutParent(T1_a));
	T1_a->cut_parent();
	um.add_event(new AddComponent(T1));
	T1->add_component(T1_a);
	//if (T1_a->get_sibling_pair_status() > 0)
	//	T1_a->clear_sibling_pair(sibling_pairs);
	//delete(T1_a);

	ContractEvent(&um, T1_a_parent);
	Node *node = T1_a_parent->contract();
	if (node != NULL && potential_new_sibling_pair && node->is_sibling_pair()){
		um.add_event(new AddToFrontSiblingPairs(sibling_pairs));
		sibling_pairs->push_front(node->rchild());
		sibling_pairs->push_front(node->lchild());
	}

	#ifdef DEBUG_APPROX
			cout << "\tT1: ";
			T1->print_components();
			cout << "\tT2: ";
			T2->print_components();
	#endif
}
if(!sibling_pairs->empty()) {
	Node *T1_a = sibling_pairs->back();
	sibling_pairs->pop_back();
	Node *T1_c = sibling_pairs->back();
	sibling_pairs->pop_back();
	um.add_event(new PopSiblingPair(T1_a, T1_c, sibling_pairs));

	//if (T1_a->get_sibling_pair_status() == 0 ||
	//		T1_c->get_sibling_pair_status() == 0) {
	//	continue;
	//}

	//T1_a->clear_sibling_pair_status();
	//T1_c->clear_sibling_pair_status();
	if (T1_a->parent() == NULL || T1_c->parent() == NULL || T1_a->parent() != T1_c->parent()) {
		continue;
	}
	Node *T1_ac = T1_a->parent();
	// lookup in T2 and determine the case
	Node *T2_a = T1_a->get_twin();
	Node *T2_c = T1_c->get_twin();

	#ifdef DEBUG_APPROX
		cout << "Fetching sibling pair" << endl;
		T1_ac->print_subtree();
		cout << "T2_a" << ": ";
		cout << " d=" << T2_a->get_depth() << " ";
		T2_a->print_subtree();
		cout << "T1_c" << ": ";
		T1_c->print_subtree();
		cout << "T2_c" << ": ";
		cout << " d=" << T2_c->get_depth() << " ";
		T2_c->print_subtree();
	#endif

	// Case 2 - Contract identical sibling pair
	if (T2_a->parent() != NULL && T2_a->parent() == T2_c->parent()) {
		#ifdef DEBUG_APPROX
			cout << "Case 2" << endl;
			T1->print_components();
			T2->print_components();
		#endif
		Node *T2_ac = T2_a->parent();
		um.add_event(new ContractSiblingPair(T1_ac));
		um.add_event(new ContractSiblingPair(T2_ac));
		T1_ac->contract_sibling_pair_undoable();
		T2_ac->contract_sibling_pair_undoable();
		um.add_event(new SetTwin(T1_ac));
		um.add_event(new SetTwin(T2_ac));
		T1_ac->set_twin(T2_ac);
		T2_ac->set_twin(T1_ac);
		//T1->add_deleted_node(T1_a);
		//T1->add_deleted_node(T1_c);
		//T2->add_deleted_node(T2_a);
		//T2->add_deleted_node(T2_c);

		// check if T2_ac is a singleton
		if (T2_ac->is_singleton() && !T1_ac->is_singleton() && T2_ac != T2->get_component(0))
			singletons->push_back(T2_ac);
		// check if T1_ac is part of a sibling pair
		if (T1_ac->parent() != NULL && T1_ac->parent()->is_sibling_pair()) {
			um.add_event(new AddToSiblingPairs(sibling_pairs));
			sibling_pairs->push_back(T1_ac->parent()->lchild());
			sibling_pairs->push_back(T1_ac->parent()->rchild());
		}
	}
	// Case 3
	else {
		#ifdef DEBUG_APPROX
			cout << "Case 3" << endl;
		#endif
		
		//  ensure T2_a is below T2_c
		if ((T2_a->get_depth() < T2_c->get_depth()
				&& T2_c->parent() != NULL)
				|| T2_a->parent() == NULL) {
			#ifdef DEBUG_APPROX
				cout << "swapping" << endl;
			#endif
		
			swap(&T1_a, &T1_c);
			swap(&T2_a, &T2_c);

		}
		else if (T2_a->get_depth() == T2_c->get_depth()) {
			if (T2_a->parent() && T2_c->parent() &&
					(T2_a->parent()->get_depth() <
					T2_c->parent()->get_depth())) {
			swap(&T1_a, &T1_c);
			swap(&T2_a, &T2_c);
			}
		}

		// get T2_b
		Node *T2_ab = T2_a->parent();
		Node *T2_b = T2_ab->rchild();
		if (T2_b == T2_a)
			T2_b = T2_ab->lchild();

		#ifdef DEBUG_APPROX
		cout << "T2_b" << ": ";
		cout.flush();
		T2_b->print_subtree();
	#endif
		// cut T1_a, T1_c, T2_a, T2_b, T2_c

		bool cut_b_only = false;
		if (T2_a->parent() != NULL && T2_a->parent()->parent() != NULL && T2_a->parent()->parent() == T2_c->parent()) {
			cut_b_only = true;
			um.add_event(new AddToSiblingPairs(sibling_pairs));
			sibling_pairs->push_back(T1_c);
			sibling_pairs->push_back(T1_a);
		}

		Node *node;

		if (!cut_b_only) {
			um.add_event(new CutParent(T1_a));
			T1_a->cut_parent();

			ContractEvent(&um, T1_ac);
			node = T1_ac->contract();

			um.add_event(new CutParent(T1_c));
			T1_c->cut_parent();


			ContractEvent(&um, node);
			node = node->contract();

			// contract parents
			// check for T1_ac sibling pair
			if (node && node->is_sibling_pair()){
				um.add_event(new AddToSiblingPairs(sibling_pairs));
				sibling_pairs->push_back(node->lchild());
				sibling_pairs->push_back(node->rchild());
			}
		}

		bool same_component = true;
		if (APPROX_CHECK_COMPONENT)
			same_component = (T2_a->find_root() == T2_c->find_root());

		Node *T2_ab_parent = T2_ab->parent();
		node = T2_ab;
		if (!cut_b_only) {
			um.add_event(new CutParent(T2_a));
			T2_a->cut_parent();

			//ContractEvent(&um, T2_ab);
			//node = T2_ab->contract();
		}
		bool cut_b = false;
		if (same_component && T2_ab_parent != NULL) {
			um.add_event(new CutParent(T2_b));
			T2_b->cut_parent();
			//ContractEvent(&um, node);
			//node = node->contract();
			cut_b = true;
		}
		// T2_b will move up after contraction
		else {
			T2_b = T2_b->parent();
		}
			ContractEvent(&um, node);
			node = node->contract();
		// check for T2 parents as singletons
		if (node != NULL && node->is_singleton()
				&& node != T2->get_component(0))
			singletons->push_back(node);

		// if T2_c is gone then its replacement is in singleton list
		// contract might delete old T2_c, see where it is
		bool add_T2_c = true;
		T2_c = T1_c->get_twin();
		// ignore T2_c if it is a singleton
		if (T2_c != node && T2_c->parent() != NULL && !cut_b_only) {

			Node *T2_c_parent = T2_c->parent();
			um.add_event(new CutParent(T2_c));
			T2_c->cut_parent();
			ContractEvent(&um, T2_c_parent);
			node = T2_c_parent->contract();
			if (node != NULL && node->is_singleton()
					&& node != T2->get_component(0))
				singletons->push_back(node);
		}
		else {
			add_T2_c = false;
		}

		
		if (!cut_b_only) {
			um.add_event(new AddComponent(T1));
			T1->add_component(T1_a);
			um.add_event(new AddComponent(T1));
			T1->add_component(T1_c);
			// put T2 cut parts into T2
			um.add_event(new AddComponent(T2));
			T2->add_component(T2_a);
			// may have already been added
		}
		if (cut_b) {
			um.add_event(new AddComponent(T2));
			T2->add_component(T2_b);
		}
		// problem if c is deleted
		if (add_T2_c) {
			um.add_event(new AddComponent(T2));
			T2->add_component(T2_c);
		}

		// may have already been added
		if (T2_b->is_leaf())
			singletons->push_back(T2_b);

		num_cut+=3;

	}
}
	}
// if the first component of the forests differ then we have cut p
if (T1->get_component(0)->get_twin() != T2->get_component(0)) {
	if (!T1->contains_rho()) {
		um.add_event(new AddRho(T1));
		um.add_event(new AddRho(T2));
		T1->add_rho();
		T2->add_rho();
	}
	else
		// hack to ignore rho when it shouldn't be in a cluster
		num_cut -=3;
}
if (save_forests) {
	*F1 = new Forest(T1);
	*F2 = new Forest(T2);
}
 

#ifdef DEBUG_UNDO
 while(um.num_events() > 0) {
		cout << "Undo step " << um.num_events() << endl;
		cout << "T1: ";
		T1->print_components();
		cout << "T2: ";
		T2->print_components();
			cout << "sibling pairs:";
			for (list<Node *>::iterator i = sibling_pairs->begin(); i != sibling_pairs->end(); i++) {
				cout << "  ";
				(*i)->print_subtree_hlpr();
			}
			cout << endl;
	 um.undo();
	cout << endl;
 }
#else
 um.undo_all();
#endif

 
//		 for(int i = 0; i < T1->num_components(); i++)
//		 	T1->get_component(i)->fix_parents();
//		 for(int i = 0; i < T2->num_components(); i++)
//		 	T2->get_component(i)->fix_parents();
return num_cut;
}


int rSPR_branch_and_bound(Forest *T1, Forest *T2) {
	return rSPR_branch_and_bound_range(T1, T2, MAX_SPR);
}


int rSPR_branch_and_bound_range(Forest *T1, Forest *T2, int end_k) {
	string problem_key;
	map<string,ProblemSolution>::iterator i;

	if (MEMOIZE) {
problem_key = T1->str() + ":" + T2->str();
i = memoized_clusters.find(problem_key);
if (i != memoized_clusters.end()) {
	//cout << "already solved: " << endl;
	//cout << problem_key << endl;
	//cout << i->second.T2 << endl;
	//cout << "start" << endl;
	Forest *new_T1 = build_finished_forest(i->second.T1);
	//cout << "middle" << endl;
	Forest *new_T2 = build_finished_forest(i->second.T2);
	//cout << "end" << endl;
	T1->swap(new_T1);
	T2->swap(new_T2);
	sync_twins(T1, T2);
	delete new_T1;
	delete new_T2;
	return i->second.k;
}
	}
	Forest F1 = Forest(T1);
	Forest F2 = Forest(T2);
	int approx_spr = rSPR_worse_3_approx(&F1, &F2);
	int min_spr = approx_spr / 3;
	int exact_spr = rSPR_branch_and_bound_range(T1, T2, min_spr, end_k);
	if (MEMOIZE && exact_spr >= 0 && i == memoized_clusters.end()) {
//string solution_key = T1->str() + ":" + T2->str();
memoized_clusters.insert(make_pair(problem_key,
		ProblemSolution(T1,T2,exact_spr)));
	}

	return exact_spr;
}
	
int rSPR_branch_and_bound_range(Forest *T1, Forest *T2, int start_k,
int end_k) {
	int exact_spr = -1;
	bool in_main = MAIN_CALL;
	MAIN_CALL = false;
	int k;
	for(k = start_k; k <= end_k; k++) {
if (in_main) {
	cout << " " << k;
	cout.flush();
}
//Forest F1 = Forest(T1);
//Forest F2 = Forest(T2);
//exact_spr = rSPR_branch_and_bound(&F1, &F2, k);
exact_spr = rSPR_branch_and_bound(T1,T2, k);
//if (exact_spr >= 0 || k == end_k) {
if (exact_spr >= 0) {
//			F1.swap(T1);
//			F2.swap(T2);
	break;
}
	}
	if (in_main)
cout << endl;
	if (k > end_k)
k = -1;
	return k;
}

int rSPR_branch_and_bound(Forest *T1, Forest *T2, int k) {
	return rSPR_branch_and_bound(T1, T2, k, NULL, NULL);
}

/* rSPR_branch_and_bound
 * Calculate a maximum agreement forest and SPR distance
 * Uses a branch and bound optimization to not explore paths
 * guaranteed to be incorrect based on rspr_3_approx
 * RETURN The rSPR distance
 * NOTE: destructive. The computed forests replace T1 and T2.
 */
int rSPR_branch_and_bound(Forest *T1, Forest *T2, int k,
		map<string, int> *label_map,
		map<int, string> *reverse_label_map) {
	// find sibling pairs of T1
//	cout << "foo1" << endl;
	if (!sync_twins(T1, T2))
return 0;
	if (PREORDER_SIBLING_PAIRS &&
			T1->get_component(0)->get_preorder_number() == -1) {
		T1->get_component(0)->preorder_number();
		T2->get_component(0)->preorder_number();
}
	if (DEEPEST_PROTECTED_ORDER
			&& T1->get_component(0)->get_edge_pre_start() == -1) {
		T1->get_component(0)->edge_preorder_interval();
		T2->get_component(0)->edge_preorder_interval();
	}

	set<SiblingPair> *sibling_pairs;
	list<Node *> singletons;
	list<pair<Forest,Forest> > AFs = list<pair<Forest,Forest> >();
	sibling_pairs = find_sibling_pairs_set(T1);
	singletons = T2->find_singletons();
	list<Node *> protected_stack = list<Node *>();
	int num_ties = 2;


	int final_k = 
rSPR_branch_and_bound_hlpr(T1, T2, k, sibling_pairs, &singletons, false, &AFs, &protected_stack, &num_ties);

//		cout << "foo" << endl;
	// TODO: this is a cheap hack
	if (!AFs.empty()) {
if (ALL_MAFS
#ifdef DEBUG
		|| true
#endif
		) {
	cout << endl << endl << "FOUND ANSWERS" << endl;
	// TODO: this is a cheap hack
	for (list<pair<Forest,Forest> >::iterator x = AFs.begin(); x != AFs.end(); x++) {
		if (label_map != NULL && reverse_label_map != NULL) {
			x->first.numbers_to_labels(reverse_label_map);
			x->second.numbers_to_labels(reverse_label_map);
		}
		cout << "\tT1: ";
		x->first.print_components();
		cout << "\tT2: ";
		x->second.print_components();
		if (label_map != NULL && reverse_label_map != NULL) {
			x->first.labels_to_numbers(label_map, reverse_label_map);
			x->second.labels_to_numbers(label_map, reverse_label_map);
		}
	}
}
AFs.front().first.swap(T1);
AFs.front().second.swap(T2);
sync_twins(T1,T2);
	}
	if (final_k >= 0)
final_k = k - final_k;
	delete sibling_pairs;
	return final_k;
}

void add_sibling_pair(set<SiblingPair> *sibling_pairs, Node *a, Node *c, UndoMachine *um) {
	SiblingPair sp = SiblingPair(a,c);
	pair< set<SiblingPair>::iterator, bool> ins = 
	sibling_pairs->insert(sp);
	if (ins.second == false) {
um->add_event(new RemoveSetSiblingPairs(sibling_pairs, *(ins.first)));
sibling_pairs->erase(ins.first);
ins = sibling_pairs->insert(sp);
	}
	um->add_event(new AddToSetSiblingPairs(sibling_pairs, *(ins.first)));
}

SiblingPair pop_sibling_pair(set<SiblingPair> *sibling_pairs, UndoMachine *um) {
	set<SiblingPair>::iterator s = sibling_pairs->begin();
	SiblingPair spair = SiblingPair(*s); 
	um->add_event(new RemoveSetSiblingPairs(sibling_pairs, spair));
	sibling_pairs->erase(s);
	return spair;
}

SiblingPair pop_sibling_pair(set<SiblingPair>::iterator s, set<SiblingPair> *sibling_pairs, UndoMachine *um) {
	SiblingPair spair = SiblingPair(*s); 
	um->add_event(new RemoveSetSiblingPairs(sibling_pairs, spair));
	sibling_pairs->erase(s);
	return spair;
}

inline int rSPR_branch_and_bound_hlpr(Forest *T1, Forest *T2, int k,
set<SiblingPair> *sibling_pairs, list<Node *> *singletons,
bool cut_b_only, list<pair<Forest,Forest> > *AFs,
list<Node *> *protected_stack, int *num_ties) {
	return rSPR_branch_and_bound_hlpr(T1, T2, k, sibling_pairs,
			singletons, cut_b_only, AFs, protected_stack, num_ties, NULL, NULL);
}

// rSPR_branch_and_bound recursive helper function
int rSPR_branch_and_bound_hlpr(Forest *T1, Forest *T2, int k,
set<SiblingPair> *sibling_pairs, list<Node *> *singletons,
bool cut_b_only, list<pair<Forest,Forest> > *AFs,
list<Node *> *protected_stack, int *num_ties, Node *prev_T1_a, Node *prev_T1_c){
	return rSprBB::rSPR_branch_and_bound_hlpr_Inline(
	PREFER_RHO,
	LEAF_REDUCTION,
	PREFER_NONBRANCHING,
	DEEPEST_ORDER,
	DEEPEST_PROTECTED_ORDER,
	CHECK_MERGE_DEPTH,
	EDGE_PROTECTION,
	CUT_LOST,
	CUT_ALL_B,
	CUT_ONE_B,
	CUT_ONE_AB,
	CUT_TWO_B,
	REVERSE_CUT_ONE_B,
	REVERSE_CUT_ONE_B_3,
	REVERSE_CUT_ONE_B_2,
	CUT_TWO_B_ROOT,
	BB,
	CLUSTER_REDUCTION,
	NUM_CLUSTERS,
	MAX_CLUSTERS,
	EDGE_PROTECTION_TWO_B,
	CUT_AC_SEPARATE_COMPONENTS,
	ABORT_AT_FIRST_SOLUTION,
	ALL_MAFS,

	&add_sibling_pair,
	&is_nonbranching,
	&pop_sibling_pair,
	&pop_sibling_pair,
	&rSPR_worse_3_approx_hlpr,
	&rSPR_branch_and_bound_range,
	&rSPR_branch_and_bound_hlpr,
	&rSPR_branch_and_bound_hlpr,
	&rSPR_branch_and_bound_hlpr,



	T1, T2, k, sibling_pairs, singletons, cut_b_only, AFs, protected_stack, num_ties, prev_T1_a, prev_T1_c);
}


// REFERS TO INLINED AWAY!
int rSPR_branch_and_bound_simple_clustering(Node *T1, Node *T2, bool verbose, map<string, int> *label_map, map<int, string> *reverse_label_map) {
	return rSPR_branch_and_bound_simple_clustering(T1,T2, verbose, label_map, reverse_label_map, -1, -1, NULL, NULL);
}
// REFERS TO INLINED AWAY - KEY A!
int rSPR_branch_and_bound_simple_clustering(Node *T1, Node *T2, bool verbose, map<string, int> *label_map, map<int, string> *reverse_label_map, int min_k, int max_k) {
	return rSPR_branch_and_bound_simple_clustering(T1,T2, verbose, label_map, reverse_label_map, min_k, max_k, NULL, NULL);
}
// INLINED AWAY TO rSprAlgorithmBB.h - REFER TO KEY A!
int rSPR_branch_and_bound_simple_clustering(Node *T1, Node *T2, bool verbose, map<string, int> *label_map, map<int, string> *reverse_label_map, int min_k, int max_k, Forest **out_F1, Forest **out_F2)
{
	return rSprBB::rSPR_branch_and_bound_simple_clustering_Inline(
	MAX_SPR, PREFER_RHO, MULTIFURCATING, CLUSTER_TUNE, COUNT_LOSSES, SHOW_CLUSTERS, MIN_SPR, CLUSTER_MAX_SPR, INITIAL_TREE_FRACTION, SPLIT_APPROX, SPLIT_APPROX_THRESHOLD, CLAMP,

	&rSPR_worse_3_mult_approx,
	&rSPR_worse_3_approx,
	&rSPR_worse_3_approx,
	&reduction_leaf_mult,
	&reduction_leaf,
	&rSPR_branch_and_bound_mult,
	&rSPR_branch_and_bound,
	&find_subtree_of_approx_distance,
	&rSPR_branch_and_bound_hlpr,

	T1, T2, verbose, label_map, reverse_label_map, min_k, max_k, out_F1, out_F2);
}

// INLINED AWAY TO rSprAlgorithmBB.h - REFER TO KEY B!
int rSPR_branch_and_bound_simple_clustering(Forest *T1, Forest *T2, bool verbose, map<string, int> *label_map, map<int, string> *reverse_label_map) {
	return rSprBB::rSPR_branch_and_bound_simple_clustering_Inline(
	MAX_SPR, PREFER_RHO, MULTIFURCATING, CLUSTER_TUNE, CLAMP,
	&rSPR_worse_3_mult_approx,
	&rSPR_worse_3_approx,
	&rSPR_branch_and_bound_mult_range,
	&rSPR_branch_and_bound_range,

	T1, T2, verbose, label_map, reverse_label_map

	);
}

/*Joel's part*/
// Refers to Inlined Away - REFER TO KEY B!
int rSPR_branch_and_bound_simple_clustering(Forest *T1, Forest *T2, bool verbose){
	return rSPR_branch_and_bound_simple_clustering(T1, T2, false, NULL, NULL);
}
// Refers to Inlined Away - REFER TO KEY B!
int rSPR_branch_and_bound_simple_clustering(Forest *T1, Forest *T2){
	return rSPR_branch_and_bound_simple_clustering(T1, T2, false, NULL, NULL);
}
// Refers to Inlined Away - REFER TO KEY A!
int rSPR_branch_and_bound_simple_clustering(Node *T1, Node *T2, bool verbose) {
	return rSPR_branch_and_bound_simple_clustering(T1, T2, false, NULL, NULL, -1, -1);
}
// Refers to Inlined Away - REFER TO KEY A!
int rSPR_branch_and_bound_simple_clustering(Node *T1, Node *T2, bool verbose, int min_k, int max_k) {
	return rSPR_branch_and_bound_simple_clustering(T1, T2, false, NULL, NULL, min_k, max_k);
}
// Refers to Inlined Away - REFER TO KEY A!

int rSPR_branch_and_bound_simple_clustering(Node *T1, Node *T2) {
	return rSPR_branch_and_bound_simple_clustering(T1, T2, false);
}

// Refers to Inlined Away - REFER TO KEY A!
int rSPR_branch_and_bound_simple_clustering(Node *T1, Node *T2, Forest **out_F1, Forest **out_F2) {
	return rSPR_branch_and_bound_simple_clustering(T1,T2, false, NULL, NULL, -1, -1, out_F1, out_F2);
}

// T1 and T2 are assumed to already be synced
void reduction_leaf_mult(Forest *T1, Forest* T2) {
  
  list<Node *> *sibling_groups = T1->find_sibling_groups();
  while (!sibling_groups->empty()) {
    //Get a sibling group with identical siblings
    list<Node*>::reverse_iterator i = sibling_groups->rbegin();
    Node *T1_sibling_group = sibling_groups->back();
    list<list<Node*>> identical_sibling_groups;
    T1_sibling_group->find_identical_sibling_groups(&identical_sibling_groups);
    for (; i != sibling_groups->rend(); i++ ){
      (*i)->find_identical_sibling_groups(&identical_sibling_groups);
      if (identical_sibling_groups.size() > 0) {
	T1_sibling_group = (*i);
	break;
      }
    }

    // Checked all sibling groups, found none with identical groups in T2
    // Therefore there are no more contractions to be made

    if (identical_sibling_groups.size() == 0) {
      delete sibling_groups;
      return;
      }
    // Contract them
    else {
      list<list<Node *>>::iterator i;
      for (i = identical_sibling_groups.begin(); i != identical_sibling_groups.end(); i++) {
	list<Node *> T2_group = (*i);
	Node *T2_p = T2_group.front()->parent();
	Node *T1_group_new = T1_sibling_group->contract_twin_group(&T2_group);
	Node *T2_group_new = T2_p->contract_sibling_group(&T2_group);
	  
	T1_group_new->set_twin(T2_group_new);
	T2_group_new->set_twin(T1_group_new);			

	if (T1_sibling_group->parent() != NULL) {
	  //Check if the contraction made a new sibling group
	  T1_sibling_group->parent()->recalculate_non_leaf_children();
	  if (T1_sibling_group->parent()->is_sibling_group()) {
	    sibling_groups->push_front(T1_sibling_group->parent());
	  }
	}
	//Check if this contraction removed a sibling group
	if (!T1_sibling_group->is_sibling_group()) {
	  sibling_groups->remove(T1_sibling_group);
	}	  
      }
    }
  }
  delete sibling_groups;
}
  
// T1 and T2 are assumed to already be synced
void reduction_leaf(Forest *T1, Forest *T2) {
	reduction_leaf(T1, T2, NULL);
}

void reduction_leaf(Forest *T1, Forest *T2, UndoMachine *um) {
	list<Node *> *sibling_pairs = T1->find_sibling_pairs();
	Node *T1_a;
	Node *T1_c;
	while (!sibling_pairs->empty()) {
		T1_a = sibling_pairs->front();
		sibling_pairs->pop_front();
		T1_c = sibling_pairs->front();
		sibling_pairs->pop_front();
		// shouldn't happen here
		if (T1_a->parent() == NULL || T1_a->parent() != T1_c->parent()) {
				continue;
		}
		Node *T2_a = T1_a->get_twin();
		Node *T2_c = T1_c->get_twin();
		if (T2_a->parent() != NULL && T2_a->parent() == T2_c->parent()) {
			Node *T1_ac = T1_a->parent();
			Node *T2_ac = T2_a->parent();
			T1_ac->contract_sibling_pair_undoable();
			Node *T2_ac_new = T2_ac->contract_sibling_pair_undoable(T2_a, T2_c);
			if (T2_ac_new != NULL && T2_ac_new != T2_ac) {
				T2_ac = T2_ac_new;
				T2_ac->contract_sibling_pair_undoable();
			}

			T1_ac->set_twin(T2_ac);
			T2_ac->set_twin(T1_ac);

			// check if T2_ac is a singleton
			// also shouldn't happen
//				if (T2_ac->is_singleton() && !T1_ac->is_singleton() && T2_ac != T2->get_component(0))

			// check if T1_ac is part of a sibling pair
			if (T1_ac->parent() != NULL &&
					T1_ac->parent()->is_sibling_pair()) {
				sibling_pairs->push_back(T1_ac->parent()->lchild());
				sibling_pairs->push_back(T1_ac->parent()->rchild());
			}
			#ifdef DEBUG
				cout << "\tT1: ";
				T1->print_components();
				cout << "\tT2: ";
				T2->print_components();
			#endif
		}
	}
	delete sibling_pairs;
}

/* return true if T1_node matches the chain between T2_node and
	 T2_node_end
*/


int rSPR_total_distance(Node *T1, vector<Node *> &gene_trees) {
	return rSPR_total_distance(T1, gene_trees, NULL);
}

int rSPR_total_distance(Node *T1, vector<Node *> &gene_trees,
		vector<int> *original_scores) {
	int total = 0;
	MAIN_CALL = false;
	int end = gene_trees.size();
//	T1->preorder_number();
	#pragma omp parallel for reduction(+ : total) firstprivate(PREFER_RHO)  // firstprivate(IN_SPLIT_APPROX)
//	for(int j = 0; j < 10; j++)
//	cout << "T1: " << T1->str_subtree() << endl;
	for(int i = 0; i < end; i++) {
			//		cout << i << endl;
	  cout << "Trying tree #" << i << " : " << gene_trees[i]->str_subtree() << endl;
	  int k = rSPR_branch_and_bound_simple_clustering(T1, gene_trees[i], VERBOSE);

	  MULTIFURCATING = true;
		//MULT_4_BRANCH = true;
		int mult_k = rSPR_branch_and_bound_simple_clustering(T1, gene_trees[i], VERBOSE);
		MULTIFURCATING = false;
		//MULT_4_BRANCH = false;
		
		if (k != mult_k) {
		  cout << "BINARY DOES NOT MATCH MULT" << endl;
		  cout << "T1: " << T1->str_subtree() << endl;;
		  cout << "BINARY k = " << k << " mult_k = " << mult_k << endl;
		  break;
		  }
		else {
		  cout << "\tMATCHES: k = " << k << endl;
		}
		
		//cout << "\t: k = " << k << endl;
//		k *= mylog2(gene_trees[i]->size());

		if (original_scores != NULL)
			(*original_scores)[i] = k;
		total += k;
//		cout << "T2: " << gene_trees[i]->str_subtree() << endl;
//		cout << " k: " << k << endl;
		if (FIND_RATE) {
			if (k > 0) {
				int size = gene_trees[i]->find_leaves().size();
//				cout << k << endl;
//				cout << size << endl;
				cout << "rate=" << (float)k / size << endl;
//				cout << T1->str_subtree() << endl;
//				cout << gene_trees[i]->str_subtree() << endl;
			}
		}
//		Forest F1 = Forest(T1);
//		Forest F2 = Forest(gene_trees[i]);
//		total += rSPR_branch_and_bound(&F1, &F2);
	}
	return total;
}

void rSPR_pairwise_distance(Node *T1, vector<Node *> &gene_trees) {
	rSPR_pairwise_distance(T1, gene_trees, 0, gene_trees.size());
}

void rSPR_pairwise_distance(Node *T1, vector<Node *> &gene_trees, bool APPROX) {
	rSPR_pairwise_distance(T1, gene_trees, 0, gene_trees.size(), APPROX);
}

void rSPR_pairwise_distance(Node *T1, vector<Node *> &gene_trees, int start, int end) {
	rSPR_pairwise_distance(T1, gene_trees, start, end, false);
}

void rSPR_pairwise_distance(Node *T1, vector<Node *> &gene_trees, int start, int end, bool approx) {
	MAIN_CALL = false;
//	T1->preorder_number();
	vector<int> distances = vector<int>(end-start);
	#pragma omp parallel for shared(distances) firstprivate(PREFER_RHO)
	for(int i = start; i < end; i++) {
		int k;
		if (approx) {
			Forest F1 = Forest(T1);
			Forest F2 = Forest(gene_trees[i]);
			k = rSPR_worse_3_approx_distance_only(&F1, &F2)/3;
		}
		else {
			k = rSPR_branch_and_bound_simple_clustering(T1, gene_trees[i]);
		}
		distances[i-start] = k;
	}

	cout << distances[0];
	for(int i = 1; i < end-start; i++) {
		cout << "," << distances[i];
	}
	cout << "\n";
}


void rSPR_pairwise_distance(Node *T1, vector<Node *> &gene_trees, int max_spr) {
	rSPR_pairwise_distance(T1, gene_trees, max_spr, 0, (int)gene_trees.size());
}

void rSPR_pairwise_distance(Node *T1, vector<Node *> &gene_trees, int max_spr, int start, int end) {
	MAIN_CALL = false;
//	T1->preorder_number();
	vector<int> distances = vector<int>(end-start);
	#pragma omp parallel for shared(distances) firstprivate(PREFER_RHO)
	for(int i = start; i < end; i++) {
		Forest F1 = Forest(T1);
		Forest F2 = Forest(gene_trees[i]);
		int k = rSPR_branch_and_bound_range(&F1, &F2, 0, max_spr);
		distances[i-start] = k;
	}

	cout << distances[0];
	for(int i = 1; i < end-start; i++) {
		cout << "," << distances[i];
	}
	cout << "\n";
}

void rSPR_pairwise_distance_unrooted(Node *T1, vector<Node *> &gene_trees) {
	rSPR_pairwise_distance_unrooted(T1, gene_trees, 0, gene_trees.size());
}

void rSPR_pairwise_distance_unrooted(Node *T1, vector<Node *> &gene_trees, bool approx) {
	rSPR_pairwise_distance_unrooted(T1, gene_trees, 0, gene_trees.size(), approx);
}

void rSPR_pairwise_distance_unrooted(Node *T1, vector<Node *> &gene_trees, int start, int end) {
	rSPR_pairwise_distance_unrooted(T1, gene_trees, 0, gene_trees.size(), false);
}

void rSPR_pairwise_distance_unrooted(Node *T1, vector<Node *> &gene_trees, int start, int end, bool approx) {
	MAIN_CALL = false;
	T1->preorder_number();
	vector<int> distances = vector<int>(end-start);
	#pragma omp parallel for shared(distances) firstprivate(PREFER_RHO)
	for(int i = start; i < end; i++) {
		int best_k = INT_MAX;
		Node *T2_copy = new Node(*(gene_trees[i]));
		vector<Node *> descendants = 
				T2_copy->find_descendants();
		for(int j = 0; j < descendants.size(); j++) {
			T2_copy->reroot(descendants[j]);
			T2_copy->set_depth(0);
			T2_copy->fix_depths();
			T2_copy->preorder_number();
	//				cout << i << "," << j << endl;
	//				cout << T1->str_subtree() << endl;
	//				cout << gene_trees[i]->str_subtree() << endl;
			int k;
			if (approx) {
				Forest F1 = Forest(T1);
				Forest F2 = Forest(T2_copy);
				k = rSPR_worse_3_approx(&F1, &F2) / 3;
			}
			else {
				k = rSPR_branch_and_bound_simple_clustering(T1, T2_copy, false);
			}
			if (k < best_k) {
				best_k = k;
			}
		}
		distances[i-start] = best_k;
		T2_copy->delete_tree();
	}

	cout << distances[0];
	for(int i = 1; i < end-start; i++) {
		cout << "," << distances[i];
	}
	cout << "\n";
}

void rSPR_pairwise_distance_unrooted(Node *T1, vector<Node *> &gene_trees, int max_spr) {
	rSPR_pairwise_distance_unrooted(T1, gene_trees, max_spr, 0, (int)gene_trees.size());
}

void rSPR_pairwise_distance_unrooted(Node *T1, vector<Node *> &gene_trees, int max_spr, int start, int end) {
	MAIN_CALL = false;
	T1->preorder_number();
	vector<int> distances = vector<int>(end-start);
	#pragma omp parallel for shared(distances) firstprivate(PREFER_RHO)
	for(int i = start; i < end; i++) {
		int best_k = -1;
		Node *T2_copy = new Node(*(gene_trees[i]));
		vector<Node *> descendants = 
				T2_copy->find_descendants();
		for(int j = 0; j < descendants.size(); j++) {
			T2_copy->reroot(descendants[j]);
			T2_copy->set_depth(0);
			T2_copy->fix_depths();
			T2_copy->preorder_number();
	//				cout << i << "," << j << endl;
	//				cout << T1->str_subtree() << endl;
	//				cout << gene_trees[i]->str_subtree() << endl;
			Forest F1 = Forest(T1);
			Forest F2 = Forest(T2_copy);
			int k = rSPR_branch_and_bound_range(&F1, &F2, 0, max_spr);
			if ((best_k == -1) || (k < best_k && k >= 0)) {
				best_k = k;
			}
		}
		distances[i-start] = best_k;
		T2_copy->delete_tree();
	}

	cout << distances[0];
	for(int i = 1; i < end-start; i++) {
		cout << "," << distances[i];
	}
	cout << "\n";
}

int rSPR_total_distance_precomputed(Node *T1, vector<Node *> &gene_trees,
		vector<int> *original_scores, vector<int> *new_original_scores, Node *old_T1) {
	int total = 0;
	MAIN_CALL = false;
	int end = gene_trees.size();
//	T1->preorder_number();
	#pragma omp parallel for reduction(+ : total) firstprivate(PREFER_RHO)  // firstprivate(IN_SPLIT_APPROX)
	for(int i = 0; i < end; i++) {
		// check that the SPR move affects the projection of T1
		Forest F1 = Forest(T1);
		Forest F2 = Forest(gene_trees[i]);
		Forest F1_old = Forest(old_T1);
		sync_twins(&F1, &F2);
		sync_twins(&F1, &F1_old);
		int k = 0;
		if (original_scores == NULL
				|| rSPR_worse_3_approx(&F1, &F1_old) > 0) {
			k = rSPR_branch_and_bound_simple_clustering(T1, gene_trees[i], VERBOSE);
		}
		else {
			k = (*original_scores)[i];
		}
		if (new_original_scores != NULL) {
			(*new_original_scores)[i] = k;
		}

		total += k;
	}
	return total;
}


int rf_total_distance(Node *T1, vector<Node *> &gene_trees) {
	int total = 0;
	int end = gene_trees.size();
	#pragma omp parallel for reduction(+ : total) firstprivate(PREFER_RHO)  // firstprivate(IN_SPLIT_APPROX)
	for(int i = 0; i < end; i++) {
			//		cout << i << endl;
		int k = rf_distance(T1, gene_trees[i]);
		total += k;
	}
	return total;
}

int rf_total_distance_unrooted(Node *T1, vector<Node *> &gene_trees) {
	int total = 0;
	int end = gene_trees.size();
	#pragma omp parallel for reduction(+ : total) firstprivate(PREFER_RHO)  // firstprivate(IN_SPLIT_APPROX)
	for(int i = 0; i < end; i++) {
		int best_k = INT_MAX;
		Node T2_copy = Node(*(gene_trees[i]));
		vector<Node *> descendants = 
				T2_copy.find_descendants();
		for(int j = 0; j < descendants.size(); j++) {
			T2_copy.reroot(descendants[j]);
			T2_copy.set_depth(0);
			T2_copy.fix_depths();
			T2_copy.preorder_number();
			int k = rf_distance(T1, &T2_copy);
			if (k < best_k) {
				best_k = k;
			}
		}
		total += best_k;
	}
	return total;
}



void rf_pairwise_distance(Node *T1, vector<Node *> &gene_trees) {
	rf_pairwise_distance(T1, gene_trees, 0, gene_trees.size());
}

void rf_pairwise_distance(Node *T1, vector<Node *> &gene_trees, int start, int end) {
	MAIN_CALL = false;
//	T1->preorder_number();
	vector<int> distances = vector<int>(end-start);
	#pragma omp parallel for shared(distances) firstprivate(PREFER_RHO)
	for(int i = start; i < end; i++) {
		int k = rf_distance(T1, gene_trees[i]);
		distances[i-start] = k;
	}

	cout << distances[0];
	for(int i = 1; i < end-start; i++) {
		cout << "," << distances[i];
	}
	cout << "\n";
}

void rf_pairwise_distance_unrooted(Node *T1, vector<Node *> &gene_trees) {
	rf_pairwise_distance_unrooted(T1, gene_trees, 0, gene_trees.size());
}

void rf_pairwise_distance_unrooted(Node *T1, vector<Node *> &gene_trees, int start, int end) {
	MAIN_CALL = false;
	T1->preorder_number();
	vector<int> distances = vector<int>(end-start);
	#pragma omp parallel for shared(distances) firstprivate(PREFER_RHO)
	for(int i = start; i < end; i++) {
		int best_k = INT_MAX;
		Node T2_copy = Node(*(gene_trees[i]));
		vector<Node *> descendants = 
				T2_copy.find_descendants();
		for(int j = 0; j < descendants.size(); j++) {
			T2_copy.reroot(descendants[j]);
			T2_copy.set_depth(0);
			T2_copy.fix_depths();
			T2_copy.preorder_number();
	//				cout << i << "," << j << endl;
	//				cout << T1->str_subtree() << endl;
	//				cout << gene_trees[i]->str_subtree() << endl;
			int k = rf_distance(T1, &T2_copy);
			if (k < best_k) {
				best_k = k;
			}
		}
		distances[i-start] = best_k;
	}

	cout << distances[0];
	for(int i = 1; i < end-start; i++) {
		cout << "," << distances[i];
	}
	cout << "\n";
}

int rSPR_total_distance(Node *T1, vector<Node *> &gene_trees, int threshold) {
	int total = 0;
	MAIN_CALL = false;
	int end = gene_trees.size();
	T1->preorder_number();
	#pragma omp parallel for reduction(+ : total) firstprivate(PREFER_RHO)  // firstprivate(IN_SPLIT_APPROX)
	for(int i = 0; i < end; i++) {
		int k = rSPR_branch_and_bound_simple_clustering(T1, gene_trees[i], VERBOSE);
//		k *= mylog2(gene_trees[i]->size());
		total += k;
//		if (total > threshold) {
//			break;
//		}
	}
	return total;
}

/*Joel's part*/
int rSPR_total_distance(Forest *T1, vector<Node *> &gene_trees){
	int total = 0;
	#pragma omp parallel for reduction(+ : total) firstprivate(PREFER_RHO)
	for(int i = 0; i < gene_trees.size(); i++) {
		Forest T2 = Forest(gene_trees[i]);
		total += rSPR_branch_and_bound_simple_clustering(&T2, T1, VERBOSE);
	}
	return total;
}

int rSPR_total_approx_distance(Forest *T1, vector<Node *> &gene_trees) {
	int total = 0;
	#pragma omp parallel for reduction(+ : total)
	for(int i = 0; i < gene_trees.size(); i++) {
		Forest F1 = Forest(T1);
		Forest F2 = Forest(gene_trees[i]);
//		cout << i << endl;
//		cout << T1->str_subtree() << endl;
//		cout << gene_trees[i]->str_subtree() << endl;
		//total += rSPR_worse_3_approx(&F2, &F1)/3;
		total += rSPR_worse_3_approx(&F2, &F1)/3;
	}
	return total;
}

int rSPR_total_distance_unrooted(Node *T1, vector<Node *> &gene_trees) {
	return rSPR_total_distance_unrooted(T1, gene_trees, INT_MAX, NULL);
}

int rSPR_total_distance_unrooted(Node *T1, vector<Node *> &gene_trees,
		int threshold) {
	return rSPR_total_distance_unrooted(T1, gene_trees, threshold, NULL);
}

int rSPR_total_distance_unrooted(Node *T1, vector<Node *> &gene_trees, int threshold, vector<int> *original_scores) {
	//cout << "rSPR_total_distance_unrooted" << endl;
	int total = 0;
	MAIN_CALL = false;
	T1->preorder_number();
	#pragma omp parallel for reduction(+ : total) firstprivate(PREFER_RHO) firstprivate(MAX_SPR) firstprivate(MIN_SPR)
	for(int i = 0; i < gene_trees.size(); i++) {
//		cout << "T1: " << T1->str_subtree() << endl;
//		cout << "T2: " << gene_trees[i]->str_subtree() << endl;
		Forest f1 = Forest(T1);
		//f1.print_components();
		Forest f2 = Forest(gene_trees[i]);
		//f2.print_components();
		if (!sync_twins(&f1, &f2))
			continue;
		if (f2.get_component(0)->get_children().size() > 2) {
			f2.get_component(0)->fixroot();
			f2.get_component(0)->set_depth(0);
			f2.get_component(0)->fix_depths();
			f2.get_component(0)->preorder_number();
		}
		//f1.print_components();
		//f2.print_components();
		int size = f2.get_component(0)->size();
		int best_distance = INT_MAX;
		int old_max = MAX_SPR;
		bool done = false;
		int NO_CLUSTER_ROUNDS=15;
//		cout << "boo" << endl;
		if (!UNROOTED_MIN_APPROX) {
//			for(int k = 0; !done; k++) {
			int best_min_spr = INT_MAX;
			vector<Node *> descendants = 
				f2.get_component(0)->find_descendants();
/*
			for(int j = 0; j < descendants.size(); j++) {
				f2.get_component(0)->reroot(descendants[j]);
				f2.get_component(0)->set_depth(0);
				f2.get_component(0)->fix_depths();
				f2.get_component(0)->preorder_number();
				Forest F1 = Forest(f1);
				Forest F2 = Forest(f2);
				int distance = rSPR_worse_3_approx(&F1, &F2)/3;
				if (distance < best_min_spr)
					best_min_spr = distance;
			}
			*/
		
			int min_spr = 0;
			if (best_min_spr < INT_MAX)
				min_spr = best_min_spr;
			for(int k = min_spr; k <= NO_CLUSTER_ROUNDS; k++) {
//			for(int k = min_spr; !done; k++) {
////				cout << k << endl;
				MIN_SPR=k;
				MAX_SPR=k;
//				Node *original_lc = f2.get_component(0)->lchild();
////					f2.print_components();
////					cout << endl;
//				vector<Node *> descendants = 
//				f2.get_component(0)->find_descendants();
				for(int j = 0; j < descendants.size(); j++) {
//					cout << "J=" << j << endl;
//					cout << i << "," << k << "," << j << endl;
//					cout << "rooting at: " << descendants[j]->str_subtree() << endl;
					//f2.get_component(0)->reroot(original_lc);
					f2.get_component(0)->reroot(descendants[j]);
					f2.get_component(0)->set_depth(0);
					f2.get_component(0)->fix_depths();
					f2.get_component(0)->preorder_number();
////					f2.print_components();
////					cout << endl;
	//			cout << T1->str_subtree() << endl;
	//			cout << gene_trees[i]->str_subtree() << endl;
					int distance;
					Forest *F1 = new Forest(f1);
					Forest *F2 = new Forest(f2);
					if (k <= NO_CLUSTER_ROUNDS)
						distance = rSPR_branch_and_bound_range(F1, F2, MIN_SPR, MAX_SPR);
//					else
//						break;
//						distance = rSPR_branch_and_bound_simple_clustering(F1->get_component(0), F2->get_component(0), VERBOSE, k, k);
					if (distance < 0)
						distance = k+1;
					delete F1;
					delete F2;
					if (distance <= k) {
						best_distance = distance;
						k=NO_CLUSTER_ROUNDS+1;
						done = true;
						break;
					}
				}
////				cout << endl;
				//f2.get_component(0)->reroot(original_lc);
//				cout << endl;
			}
			MAX_SPR=old_max;
			MIN_SPR=0;
			if (!done) {
				vector<Node *> descendants = 
					f2.get_component(0)->find_descendants();
				for(int j = 0; j < descendants.size(); j++) {
					f2.get_component(0)->reroot(descendants[j]);
					f2.get_component(0)->set_depth(0);
					f2.get_component(0)->fix_depths();
					f2.get_component(0)->preorder_number();
	//				cout << i << "," << j << endl;
	//				cout << T1->str_subtree() << endl;
	//				cout << gene_trees[i]->str_subtree() << endl;
					int distance = rSPR_branch_and_bound_simple_clustering(f1.get_component(0), f2.get_component(0), VERBOSE);
					if (distance <= best_distance) {
							best_distance = distance;
					}
				}
			}
	//		cout << "best_distance: " << best_distance << endl;
			if (best_distance == INT_MAX)
				best_distance = 0;
			total += best_distance;
			if (original_scores != NULL)
				(*original_scores)[i] = best_distance;
	//		cout << "total: " << total << endl;
		}
		else {
			int best_approx = INT_MAX;
			Node *best_rooting = f2.get_component(0)->lchild();
			int num_ties = 2;
			vector<Node *> descendants = 
				f2.get_component(0)->find_descendants();
			int NUM_ROOTINGS = 0;
//			int NUM_ROOTINGS = 6;
//			if (descendants.size() > 70)
//				NUM_ROOTINGS = descendants.size() / 10;
//			NUM_ROOTINGS = sqrt(descendants.size());
			if (NUM_ROOTINGS > 0 && descendants.size() < NUM_ROOTINGS + 1) {
				vector<Node *> rand_descendants =
					random_select(descendants, NUM_ROOTINGS);
				descendants = rand_descendants;
				descendants.push_back(f2.get_component(0)->lchild());
			}
			for(int j = 0; j < descendants.size(); j++) {
				f2.get_component(0)->reroot(descendants[j]);
					f2.get_component(0)->set_depth(0);
					f2.get_component(0)->fix_depths();
					f2.get_component(0)->preorder_number();
				//Forest F1 = Forest(f1);
				//Forest F2 = Forest(f2);
				int distance = rSPR_worse_3_approx_distance_only(&f1, &f2)/3;
				if (distance < best_approx) {
					best_approx = distance;
					best_rooting = descendants[j];
					num_ties = 2;
				}
				else if (distance == best_approx) {
					int r = rand();
					if (r < RAND_MAX/num_ties) {
						best_approx = distance;
						best_rooting = descendants[j];
					}
					num_ties++;
				}
			}
			f2.get_component(0)->reroot(best_rooting);
					f2.get_component(0)->set_depth(0);
					f2.get_component(0)->fix_depths();
					f2.get_component(0)->preorder_number();
			int k;
			if (best_approx > 20)
				k = rSPR_branch_and_bound_simple_clustering(f1.get_component(0), f2.get_component(0), VERBOSE);
			else
					k = rSPR_branch_and_bound_range(&f1, &f2, best_approx/3, best_approx);
			total += k;
		if (original_scores != NULL)
			(*original_scores)[i] = k;
		}
//		if (total > threshold)
//			break;
	}
	return total;
}

int rSPR_total_approx_distance_unrooted(Node *T1, vector<Node *> &gene_trees) {
	int total = 0;
	MAIN_CALL = false;
	#pragma omp parallel for reduction(+: total)
	for(int i = 0; i < gene_trees.size(); i++) {
		Forest f1 = Forest(T1);
		Forest f2 = Forest(gene_trees[i]);
		if (!sync_twins(&f1, &f2))
			continue;
		if (f2.get_component(0)->get_children().size() > 2) {
			f2.get_component(0)->fixroot();
			f2.get_component(0)->set_depth(0);
			f2.get_component(0)->fix_depths();
			f2.get_component(0)->preorder_number();
		}
		int size = f2.get_component(0)->size();
		int best_distance = INT_MAX;
		vector<Node *> descendants = 
			f2.get_component(0)->find_descendants();
		for(int j = 0; j < descendants.size(); j++) {
			f2.get_component(0)->reroot(descendants[j]);
					f2.get_component(0)->set_depth(0);
					f2.get_component(0)->fix_depths();
					f2.get_component(0)->preorder_number();
			Forest F1 = Forest(f1);
			Forest F2 = Forest(f2);

			int distance = rSPR_worse_3_approx(&F1, &F2)/3;
			if (distance < best_distance)
				best_distance = distance;
		}
		if (best_distance == INT_MAX)
			best_distance = 0;
		total += best_distance;
	}
	return total;
}

int rSPR_total_approx_distance(Node *T1, vector<Node *> &gene_trees) {
	return rSPR_total_approx_distance(T1, gene_trees, INT_MAX);
}

int rSPR_total_approx_distance(Node *T1, vector<Node *> &gene_trees,
		int threshold) {
	int total = 0;
	MAIN_CALL = false;
	#pragma omp parallel for reduction(+ : total)
	for(int i = 0; i < gene_trees.size(); i++) {
		Forest F1 = Forest(T1);
		Forest F2 = Forest(gene_trees[i]);
//		cout << i << endl;
//		cout << T1->str_subtree() << endl;
//		cout << gene_trees[i]->str_subtree() << endl;
		total += rSPR_worse_3_approx(&F1, &F2)/3;
//		if (total > threshold)
//			break;
	}
	return total;
}




// --------------------------------------------------------------------
// FROM HERE ON OUT; MOST IS MOVED TO rSprUtility.h!
// --------------------------------------------------------------------

string itos(int i) {
	return rSprUtility::itos_Inline(i);
}

Node *find_subtree_of_approx_distance_hlpr(Node *n, Forest *F1, Forest *F2, int target_size) {
	return rSprUtility::find_subtree_of_approx_distance_hlpr_Inline(&find_subtree_of_approx_distance_hlpr,&rSPR_worse_3_approx, n, F1, F2, target_size);
}

// Moved to rSprUtility.
Node *find_subtree_of_approx_distance(Node *n, Forest *F1, Forest *F2, int target_size) {
	return rSprUtility::find_subtree_of_approx_distance_Inline(
		&find_subtree_of_approx_distance_hlpr, &rSPR_worse_3_approx, n, F1, F2, target_size);
}


// Moved to rSprUtility!
Node *find_best_root(Node *T1, Node *T2, double *best_root_b_acc) {
	return rSprUtility::find_best_root_Inline(&find_best_root_hlpr, T1, T2, best_root_b_acc);
}



Node *find_best_root(Node *T1, Node *T2) {
	return rSprUtility::find_best_root_Inline(&find_best_root, T1,T2);
}


// Moved to rSprUtlity
double find_best_root_acc(Node *T1, Node *T2) {
	return rSprUtility::find_best_root_acc_Inline(&find_best_root, T1, T2);
}


// Moved to rSprUtility!
void find_best_root_hlpr(Node *T2, int pre_separator, int group_1_total,
		int group_2_total, Node **best_root, double *best_root_b_acc) {
	rSprUtility::find_best_root_hlpr_Inline(
		&find_best_root_hlpr,
		T2,pre_separator,group_1_total, group_2_total, best_root, best_root_b_acc);
}



/*	class child_ba_comp {
		private:
			int g_1_total;
			int g_2_total;
			bool d;
		public:
			child_ba_comp(int group_1_total, int group_2_total, bool direction) {
				g_1_total = group_1_total;
				g_2_total = group_2_total;
				d = direction;
			}
			bool operator()(const pair<int,int> x,const pair<int,int> y) {
				if (d) {
					return ((x.first / x.second * (x.first + x.second)) - 
							(y.first / y.second * (y.first + y.second)));
				}
				else {
					return ((x.second / x.first * (x.first + x.second)) - 
							(y.second / y.first * (y.first + y.second)));
				}
			}
	};
	*/


//P Moved to rSprUtility
void find_best_root_hlpr(Node *n, int pre_separator, int group_1_total,
		int group_2_total, Node **best_root, double *best_root_b_acc,
		int *p_group_1_descendants, int *p_group_2_descendants, int *num_ties) {
	rSprUtility::find_best_root_hlpr_Inline(
	&find_best_root_hlpr,
	n, pre_separator, group_1_total,
	group_2_total,best_root, best_root_b_acc,
	p_group_1_descendants, p_group_2_descendants, num_ties);
}



//P Moved to rSprUtility
Node *find_random_root(Node *T1, Node *T2) {
	return rSprUtility::find_random_root_Inline(T1,T2);
}

//P Moved to rSprUtility!
Node *find_best_root_rspr(Node *T1, Node *T2) {
	return rSprUtility::find_best_root_rspr_Inline(&rSPR_branch_and_bound_simple_clustering, T1, T2);
}




// assume already sync_twins and preorder numbered
// Moved to rSprUtility!
bool contains_bipartition(Node *n, int pre_start, int pre_end,
		int group_1_total, int group_2_total, int *p_group_1_descendants,
		int *p_group_2_descendants) {
	return rSprUtility::contains_bipartition_Inline(&contains_bipartition,
		n, pre_start, pre_end, group_1_total, group_2_total, p_group_1_descendants, p_group_2_descendants);
}


// ToDo; Apply Inline Calculus
// root the tree based on an outgroup
// returns false if the outgroup is not found or not a clade
bool outgroup_root(Node *T, set<string, StringCompare> outgroup) {
	vector<int> num_in = vector<int>();
	vector<int> num_out = vector<int>();
	count_in_out(T, num_in, num_out, outgroup);
	list<Node *>::iterator c;
	int pre = T->get_preorder_number();
	bool clean_split = true;
	if (num_out[pre] == 0)
		return false;
	for(c = T->get_children().begin(); c != T->get_children().end(); c++) {
		int c_pre = (*c)->get_preorder_number();
		if (num_out[pre] == num_out[c_pre]) {
			// split the out_group
			return outgroup_root(*c, num_in, num_out);
		}
		else if (num_in[pre] == num_in[c_pre]) {
			// split the in_group
			return outgroup_root(*c, num_out, num_in);
		}
		else if (num_out[c_pre] > 0 && num_in[c_pre] > 0)
			clean_split = false;
	}
	if (clean_split)
		return outgroup_reroot(T, num_in, num_out);
	else
		return false;
}


// IS NOT TRIVIALLY MOVED
bool outgroup_root(Node *n, vector<int> &num_in, vector<int> &num_out) {
	list<Node *>::iterator c;
	int pre = n->get_preorder_number();
	bool clean_split = true;
	for(c = n->get_children().begin(); c != n->get_children().end(); c++) {
		int c_pre = (*c)->get_preorder_number();
		if (num_out[pre] == num_out[c_pre]) {
			// split the out_group
			return outgroup_root(*c, num_in, num_out);
		}
		else if (num_out[c_pre] > 0 && num_in[c_pre] > 0)
			clean_split = false;
	}
	if (clean_split)
		return outgroup_reroot(n, num_in, num_out);
	else
		return false;
}







// IS NOT TRIVIALLY MOVED
bool outgroup_reroot(Node *n, vector<int> &num_in, vector<int> &num_out) {
	Node *T = n->find_root();
	if (num_in[n->get_preorder_number()] == 0) {
		if (!n->is_leaf()) {
			T->reroot(n);
			T->set_depth(0);
			T->fix_depths();
			T->preorder_number();
		}
		return true;
	}
	Node *new_split = new Node("");
	n->add_child(new_split);
	list<Node *>::iterator c = n->get_children().begin();
	while(c != n->get_children().end()) {
		Node *child = *c;
		c++;
		int c_pre = child->get_preorder_number();
		if (c_pre >= 0 && num_in[c_pre] == 0) {
			//child->cut_parent();
			new_split->add_child(child);
		}
	}
	T->reroot(new_split);
	T->set_depth(0);
	T->fix_depths();
	T->preorder_number();
	return true;
}

// Moved to rSpr Utility!
void count_in_out(Node *n, vector<int> &num_in, vector<int> &num_out,
		set<string, StringCompare> &outgroup) {
	rSprUtility::count_in_out_Inline(&count_in_out,n, num_in, num_out, outgroup);
}


// STUCK DUE TO ENUMERATION!
void modify_bipartition_support(Node *T1, Node *T2, enum RELAXATION relaxed) {
	Forest F1 = Forest(T1);
	Forest F2 = Forest(T2);
	if (!sync_twins(&F1, &F2)) {
		return;
	}
	Node *n = F1.get_component(0);
	vector<int> *F1_descendant_counts = n->find_leaf_counts();
	modify_bipartition_support(n, &F1, &F2, T1, T2, F1_descendant_counts,
			relaxed);
	delete F1_descendant_counts;
}

// STUCK DUE TO ENUMERATION!
void modify_bipartition_support(Node *n, Forest *F1, Forest *F2,
		Node *T1, Node *T2, vector<int> *F1_descendant_counts, enum RELAXATION relaxed) {
	if (n->is_leaf())
		return;
	list<Node *>::iterator c;
	for(c = n->get_children().begin(); c != n->get_children().end(); c++) {
		modify_bipartition_support(*c, F1, F2, T1, T2, F1_descendant_counts,
				relaxed);
	}
	Node *t = T1->find_by_prenum(n->get_preorder_number());
	if (t->parent() == NULL)
		return;
	int pre_start = n->get_edge_pre_start();
	int pre_end = n->get_edge_pre_end();
	int group_1_total = (*F1_descendant_counts)[n->get_preorder_number()];
	int group_2_total = (*F1_descendant_counts)[F1->get_component(0)->get_preorder_number()] - group_1_total;
	if (group_2_total >= 2) {
		if (contains_bipartition(F2->get_component(0), pre_start, pre_end,
				group_1_total, group_2_total, NULL, NULL)) {
			t->a_inc_support();
			t->a_inc_support_normalization();
			// relaxed
			if (false && relaxed == ALL_RELAXED) {
				int stop_pre = 0;
				if (n->parent() != NULL) {
					stop_pre = n->parent()->get_preorder_number();
				while ((t = t->parent()) != NULL
						&& t->get_preorder_number() != stop_pre)
					t->a_inc_support();
					t->a_inc_support_normalization();
				}
			}
		}
		else {
//			t->a_dec_support();
			t->a_inc_support_normalization();
			// relaxed
			if (relaxed == ALL_RELAXED || relaxed == NEGATIVE_RELAXED) {
				int stop_pre = 0;
				if (n->parent() != NULL) {
					stop_pre = n->parent()->get_preorder_number();
				while ((t = t->parent()) != NULL
						&& t->get_preorder_number() != stop_pre)
//					t->a_dec_support();
					t->a_inc_support_normalization();
				}
			}
		}
	}
}

// Moved to rSprUtility!
int rf_distance(Node *T1, Node *T2) {
	return rSprUtility::rf_distance(&count_differing_bipartitions, T1,T2);
}


// Moved TO rSprUtility!
int count_differing_bipartitions(Node *n) {
	return rSprUtility::count_differing_bipartitions_Inline(&count_differing_bipartitions, n);
}




// MOVED TO rSprUtility!
bool is_nonbranching(Forest *T1, Forest *T2, Node *T1_a, Node *T1_c, Node *T2_a, Node *T2_c) {
    return rSprUtility::is_nonbranching_Inline(
    CUT_ONE_B, CUT_TWO_B, CUT_TWO_B_ROOT, REVERSE_CUT_ONE_B, REVERSE_CUT_ONE_B_2,
    T1, T2, T1_a, T1_c, T2_a, T2_c
    );

}



// MOVED TO rSprUtility!
void strip_whitespace(string &str) {
	rSprUtility::strip_whitespace_Inline(str);
}


// MOVED TO rSprUtility!
void strip_trailing_whitespace(string &str) {
	rSprUtility::strip_trailing_whitespace_Inline(whitespaces,str);
}

// MOVED TO rSprUtility!
//randomizes T2 with count number of sprs. If T1 equals T2 at the beginning,
//then the spr distance between T1 and T2 is equal to (or possibly less than) count
//Assumes they are already synced, assumes there are valid sprs to be made
void randomize_tree_with_spr(Forest* T1, Forest* T2, int count) {

	rSprUtility::randomize_tree_with_spr_Inline(T1, T2, count);
}
