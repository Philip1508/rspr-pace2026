
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
#include <functional>



#include "../../DataStructures/Forest.h"
#include "../../DataStructures/ClusterForest.h"
#include "../..//UndecipheredComponents/LCA.h"
#include "../..//DataStructures/ClusterInstance.h"
#include "../..//DataStructures/SiblingPair.h"
#include "../..//DataStructures/UndoMachine.h"

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


namespace rSprAlgorithm3Approx {



/*	__attribute__((always_inline)) inline int rSPR_3_approx_hlpr_Inline(
		bool &APPROX_CHECK_COMPONENT,
		Forest *T1, Forest *T2, list<Node *> *singletons,
list<Node *> *sibling_pairs) {
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
*/


__attribute__((always_inline)) inline int rSPR_worse_3_mult_approx_hlpr_Inline(Forest *T1, Forest *T2, list<Node *> *singletons, list<Node *> *sibling_groups, Forest **F1, Forest **F2, bool &save_forests) {

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




}