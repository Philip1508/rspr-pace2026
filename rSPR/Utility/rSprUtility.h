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

#include "../rspr.h"


#include "../../DataStructures/Forest.h"
#include "../../DataStructures/ClusterForest.h"
#include "../..//UndecipheredComponents/LCA.h"
#include "../..//DataStructures/ClusterInstance.h"
#include "../..//DataStructures/SiblingPair.h"
#include "../..//DataStructures/UndoMachine.h"



namespace rSprUtility
{




	bool chain_match(Node *T1_node, Node *T2_node, Node *T2_node_end) {
		Node *T1_pendant;
		Node *T2_pendant;
		bool pendant_found = false;
		if (T2_node->is_leaf())
			return false;
		// T1_node is a leaf
		if (T1_node->is_leaf()) {
			T1_pendant = T1_node;
			if (T1_pendant->get_twin() == T2_node->lchild()) {
				if (T2_node->rchild() == T2_node_end)
					return true;
			}
			else if (T1_pendant->get_twin() == T2_node->rchild()) {
				if (T2_node->lchild() == T2_node_end)
					return true;
			}
			return false;
		}
		// T1_pendant is T1_node->lchild()
		T1_pendant = T1_node->lchild();
		if (T1_pendant->is_leaf()) {
			T2_pendant = T2_node->lchild();
			if (T2_pendant->is_leaf() && T1_pendant->get_twin() == T2_pendant) {
				return chain_match(T1_pendant->get_sibling(),
						T2_pendant->get_sibling(), T2_node_end);
			}
			T2_pendant = T2_node->rchild();
			if (T2_pendant->is_leaf() && T1_pendant->get_twin() == T2_pendant) {
				return chain_match(T1_pendant->get_sibling(),
						T2_pendant->get_sibling(), T2_node_end);
			}
		}
		// T1_pendant is T1_node->rchild()
		if (T1_pendant->is_leaf()) {
			T2_pendant = T2_node->lchild();
			if (T2_pendant->is_leaf() && T1_pendant->get_twin() == T2_pendant) {
				return chain_match(T1_pendant->get_sibling(),
						T2_pendant->get_sibling(), T2_node_end);
			}
			T2_pendant = T2_node->rchild();
			if (T2_pendant->is_leaf() && T1_pendant->get_twin() == T2_pendant) {
				return chain_match(T1_pendant->get_sibling(),
						T2_pendant->get_sibling(), T2_node_end);
			}
		}
		return false;
	}

 __attribute__((always_inline)) inline bool is_nonbranching_Inline(
 	bool &CUT_ONE_B, bool &CUT_TWO_B, bool &CUT_TWO_B_ROOT, bool &REVERSE_CUT_ONE_B, bool &REVERSE_CUT_ONE_B_2,
 	Forest *T1, Forest *T2, Node *T1_a, Node *T1_c, Node *T2_a, Node *T2_c) {
	if ((T2_a->get_depth() < T2_c->get_depth()
			&& T2_c->parent() != NULL)
			|| T2_a->parent() == NULL) {
		swap(&T1_a, &T1_c);
		swap(&T2_a, &T2_c);
	}
	else if (T2_a->get_depth() == T2_c->get_depth()) {
		if (T2_a->parent() && T2_c->parent() &&
				(T2_a->parent()->get_depth() <
				T2_c->parent()->get_depth()
				//|| (T2_c->parent()->parent()
				//&& T2_c->parent()->parent() == T2_a->parent())
				)) {
		swap(&T1_a, &T1_c);
		swap(&T2_a, &T2_c);
		}
	}
	int num_protected = T2_a->is_protected() + T2_c->is_protected();
	if (T2_a->parent()->get_children().size() == 2)
		num_protected += T2_a->get_sibling()->is_protected();
	if (num_protected >= 2)
		return true;
	if (CUT_ONE_B) {
		if (T2_a->parent()->parent() == T2_c->parent()
			&& T2_c->parent() != NULL
			&& T2_a->parent()->get_children().size() <= 2)
			return true;
	}
	if (CUT_TWO_B && T1_a->parent()->parent() != NULL) {
		Node *T1_s = T1_a->parent()->get_sibling();
		if (T1_s->is_leaf()) {
			Node *T2_l = T2_a->parent()->parent();
			if (T2_l != NULL && T2_l->get_children().size() <= 2) {
				if (T2_c->parent() != NULL && T2_c->parent()->parent() == T2_l
						&& ((T2_a->parent()->get_children().size() <= 2
						&& T2_c->parent()->get_children().size() <= 2)
						|| T1_s->get_twin()->is_protected())){
					if (T2_l->get_sibling() == T1_s->get_twin()) {
						return true;
					}
					else if (CUT_TWO_B_ROOT && T2_l->parent() == NULL &&
							(T2->contains_rho() ||
							 T2->get_component(0) != T2_l)) {
						return true;
					}
				}
				else if ((T2_l = T2_l->parent()) != NULL
						&& T2_c->parent() == T2_l
						&& ((T2_a->parent()->get_children().size() <= 2
						&& T2_a->parent()->parent()->get_children().size() <= 2
						&& T2_l->get_children().size() <= 2)
						|| T1_s->get_twin()->is_protected())){
					if (T2_l->get_sibling() == T1_s->get_twin()) {
						return true;
					}
					else if (CUT_TWO_B_ROOT && T2_l->parent() == NULL &&
							(T2->contains_rho() ||
							 T2->get_component(0) != T2_l)) {
						return true;
					}
				}
			}
		}
	}
	if (REVERSE_CUT_ONE_B && T1_a->parent()->parent() != NULL) {
		Node *T1_s = T1_a->parent()->get_sibling();
		Node *T2_s = T1_s->get_twin();
		if (T1_s->is_leaf()) {
			if (T2_s->parent() == T2_a->parent()) {
				return true;
			}
			else if (T2_s->parent() == T2_c->parent()
					&& T2_c->parent()->get_children().size() <= 2) {
				return true;
			}
//			else if (REVERSE_CUT_ONE_B_3
//							&& T2_s->is_protected()
//							&& T2_s->parent() != NULL
//							&& T2_s->parent()->parent() == T2_a->parent()
//							&& T2_s->parent()->get_children().size() <= 2) {
//				return true;
//			}
		}
		else if (REVERSE_CUT_ONE_B_2 && T2_c->parent() != NULL
				&& chain_match(T1_s, T2_c->get_sibling(), T2_a))
			return true;
	}
	return false;
}


  __attribute__((always_inline)) inline void strip_whitespace_Inline(string &str) {
    std::string::iterator end_pos = std::remove_if(str.begin(), str.end(), ::isspace);
    str.erase(end_pos, str.end());
  }

  __attribute__((always_inline)) inline void strip_trailing_whitespace_Inline(const string &whitespaces,string &str)
  {
	  size_t start_pos = str.find_first_not_of(whitespaces);
	  str.erase(0, start_pos);
	  size_t end_pos = str.find_last_not_of(whitespaces)+1;
	  str.erase(end_pos, str.size()-end_pos);
  }

  __attribute__((always_inline)) inline void randomize_tree_with_spr_Inline(Forest* T1, Forest* T2, int count)
    {
      for (int spr = 0; spr < count; spr++) {
        vector<Node*> all_nodes = T2->get_component(0)->find_nodes_in_subtree();
        Node* source = all_nodes[rand() % all_nodes.size()];
        Node* target = all_nodes[rand() % all_nodes.size()];
        bool target_in_subtree = false;
        Node* first_leaf = target->find_leaves()[0];
        vector<Node*> source_leaves = source->find_leaves();
        for (int i = 0; i < source_leaves.size(); i++) {
          if (source_leaves[i] == first_leaf) {
        target_in_subtree = true;
        break;
          }
        }

        //Get random node
        //if target is in source's subtree repick
        //spr
        while (source == target ||
           source->is_sibling_of(target) ||
           source->parent() == target ||
           target_in_subtree) {
          source = all_nodes[rand() % all_nodes.size()];
          target = all_nodes[rand() % all_nodes.size()];
          //cout << "Trying: " << source->str_subtree() << " and "<<  target->str_subtree() << endl;
          target_in_subtree = false;
          Node* first_leaf = target->find_leaves()[0];

          int child_count = source->get_children().size();
          if (child_count > 2) {
        int rand_count = rand() % (child_count + 1);
        if (rand_count == child_count || rand_count == 0){
          //cout << "moving whole tree" << endl;
        }
        else if (rand_count == 1) {
          source = source->get_children().front();
          //cout << "moving first child" << endl;
        }
        else {
          list<Node*> to_expand = list<Node*>();
          list<Node*>::iterator c = source->get_children().begin();
          for (int i = 0; i < rand_count; i++) {
            to_expand.push_back(*c);
            c++;
          }
          source = source->expand_children_out(to_expand);
          //cout << "Moving part " << rand_count<< endl;
        }
          }

          vector<Node*> source_leaves = source->find_leaves();
          for (int i = 0; i < source_leaves.size(); i++) {
        if (source_leaves[i] == first_leaf) {
          //cout << "target in subtree" << endl;
          target_in_subtree = true;
          break;
        }
          }
        }

        //cout << "Moving : " << source->str_subtree() << " to " << target->str_subtree() << endl;

        Node* parent = source->parent();
        if (parent != NULL) {
          source->cut_parent();
          if (parent->get_children().size() == 1) {
        parent->contract(true);
          }
        }
        else {
          continue;
        }
        //regraft
        if (target->parent() != NULL) {
          target->parent()->add_child(source);
        }
        else {
          Node* new_parent = target;//new Node();
          Node* replace = new Node(*target);
          new_parent->get_children().clear();
          new_parent->add_child(replace);
          new_parent->add_child(source);
        }

        //T1->print_components();
        //T2->print_components();
      }
    }




}