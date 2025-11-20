#pragma once

/*#include <cstdio>
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
*/

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




	__attribute__((always_inline)) inline bool contains_bipartition_Inline(
	bool (*contains_bipartition)(Node *n, int pre_start, int pre_end,
		int group_1_total, int group_2_total, int *p_group_1_descendants,
		int *p_group_2_descendants),
	Node *n, int pre_start, int pre_end, int group_1_total, int group_2_total, int *p_group_1_descendants, int *p_group_2_descendants) {
		list<Node*>::iterator c;
		int group_1_descendants = 0;
		int group_2_descendants = 0;
		bool found = false;
		bool proper_split = true;
		for(c = n->get_children().begin(); c != n->get_children().end(); c++) {
			int c_group_1_descendants = 0;
			int c_group_2_descendants = 0;
			found = contains_bipartition(*c, pre_start, pre_end, group_1_total,
					group_2_total, &c_group_1_descendants, &c_group_2_descendants);
			if (found)
				return true;
			group_1_descendants += c_group_1_descendants;
			group_2_descendants += c_group_2_descendants;
			if (c_group_1_descendants > 0 && c_group_2_descendants > 0)
				proper_split = false;
		}

		if (n->is_leaf()) {
			int pre = n->get_twin()->get_preorder_number();
			if (pre >= pre_start && pre <= pre_end)
				group_1_descendants++;
			else
				group_2_descendants++;
		}
		else if (proper_split) {
			if (group_1_descendants == group_1_total
					|| group_2_descendants == group_2_total)
				return true;
		}
		if (p_group_1_descendants != NULL)
			*p_group_1_descendants += group_1_descendants;
		if (p_group_2_descendants != NULL)
			*p_group_2_descendants += group_2_descendants;
		return false;
	}


	// Here should be outgroop_root etc.

	// Does not work.
	// __attribute__((always_inline)) inline bool outgroup_reroot_Inline(Node *n, vector<int> &num_in, vector<int> &num_out) {

	__attribute__((always_inline)) inline void count_in_out_Inline(
	void (*count_in_out)(Node *n, vector<int> &num_in, vector<int> &num_out,
	set<string, StringCompare> &outgroup),
		Node *n, vector<int> &num_in, vector<int> &num_out,
		set<string, StringCompare> &outgroup) {
		list<Node *>::iterator c;
		int pre = n->get_preorder_number();
		if (num_in.size() <= pre)
			num_in.resize(pre + 1, 0);
		if (num_out.size() <= pre)
			num_out.resize(pre + 1, 0);
		if (n->is_leaf()) {
			if (outgroup.find(n->get_name()) != outgroup.end()) {
				num_out[pre] = 1;
				num_in[pre] = 0;
			}
			else {
				num_out[pre] = 0;
				num_in[pre] = 1;
			}
		}
		else {
			for(c = n->get_children().begin(); c != n->get_children().end(); c++) {
				count_in_out(*c, num_in, num_out, outgroup);
				num_in[pre] += num_in[(*c)->get_preorder_number()];
				num_out[pre] += num_out[(*c)->get_preorder_number()];
			}
		}
	}







	__attribute__((always_inline)) inline int rf_distance(
		int (*count_differing_bipartitions)(Node* n)
		,Node *T1, Node *T2
		) {
		Forest F1 = Forest(T1);
		Forest F2 = Forest(T2);
		if (!sync_twins(&F1, &F2))
			return 0;
		if (F1.get_component(0)->is_leaf())
			return 0;
		sync_interior_twins(&F1, &F2);
		int rf_d = 0;
		rf_d += count_differing_bipartitions(F1.get_component(0));
		rf_d += count_differing_bipartitions(F2.get_component(0));
		return rf_d;
	}

	__attribute__((always_inline)) inline int count_differing_bipartitions_Inline(
		int (*count_differing_bipartitions)(Node* n),
		Node *n) {
		//cout << "Start: " << n->str_subtree() << endl;
		int count = 0;
		list<Node *>::iterator c;
		for(c = n->get_children().begin(); c != n->get_children().end(); c++) {
			count += count_differing_bipartitions(*c);
		}
		if (n->get_twin() == NULL ||
	//			n->get_depth() > n->get_twin()->get_twin()->get_depth())
				n != n->get_twin()->get_twin()) {
			count++;
				}
		return count;
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