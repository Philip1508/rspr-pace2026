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


	__attribute__((always_inline)) inline string itos_Inline(int i) {
		stringstream ss;
		string a;
		ss << i;
		a = ss.str();
		return a;
	}


	__attribute__((always_inline)) inline Node *find_subtree_of_approx_distance_hlpr_Inline(
	// Shadowed Params
	Node *(*find_subtree_of_approx_distance_hlpr)(Node *n, Forest *F1, Forest *F2, int target_size),

	int (*rSPR_worse_3_approx)(Node *subtree, Forest *T1, Forest *T2),

	//Original Params
	Node *n, Forest *F1, Forest *F2, int target_size) {
	Node *largest_child_subtree = NULL;
	int lcs_size = 0;
	list<Node *>::iterator c;
	for(c = n->get_children().begin(); c != n->get_children().end(); c++) {
		Forest f1 = Forest(F1);
		Forest f2 = Forest(F2);
		Node *subtree = f1.find_by_prenum((*c)->get_preorder_number());
		f1.get_component(0)->disallow_siblings_subtree();
		subtree->allow_siblings_subtree();
//		if (subtree->lchild() != NULL)
//			subtree->lchild()->allow_siblings_subtree();
//		if (subtree->rchild() != NULL)
//			subtree->rchild()->allow_siblings_subtree();

		int cs_size = rSPR_worse_3_approx(subtree, &f1, &f2);
		if (cs_size > lcs_size) {
			largest_child_subtree = *c;
			lcs_size = cs_size;
		}
	}
	if (lcs_size <= 0)
		return n;
	else if (lcs_size < target_size)
		return largest_child_subtree;
	else
		return find_subtree_of_approx_distance_hlpr(largest_child_subtree,
				F1, F2, target_size);
}

	__attribute__((always_inline)) inline Node *find_subtree_of_approx_distance_Inline(
	// SHADOWED PARAMS
	Node *(*find_subtree_of_approx_distance_hlpr)(Node *n, Forest *F1, Forest *F2, int target_size),
	int (*rSPR_worse_3_approx)(Node *subtree, Forest *T1, Forest *T2),
	// ORIGINAL PARAMS
	Node *n, Forest *F1, Forest *F2, int target_size) {
		Forest f1 = Forest(F1);
		Forest f2 = Forest(F2);
		Node *subtree = f1.find_by_prenum(n->get_preorder_number());
		f1.get_component(0)->disallow_siblings_subtree();
		subtree->allow_siblings_subtree();
		//		if (subtree->lchild() != NULL)
		//			subtree->lchild()->allow_siblings_subtree();
		//		if (subtree->rchild() != NULL)
		//			subtree->rchild()->allow_siblings_subtree();
		int size = rSPR_worse_3_approx(subtree, &f1, &f2);
		if (size > target_size)
			return find_subtree_of_approx_distance_hlpr(n, F1, F2, target_size);
		else
			return n;
	}

	__attribute__((always_inline)) inline Node *find_best_root_Inline(
	// SHADOWED PARAMETERS
	void (*find_best_root_hlpr)(Node *T2, int pre_separator, int group_1_total,
		int group_2_total, Node **best_root, double *best_root_b_acc),
	// ORIGINAL PARAMETERS
	Node *T1, Node *T2, double *best_root_b_acc) {
		Forest F1 = Forest(T1);
		Forest F2 = Forest(T2);
		Node *t1 = F1.get_component(0);
		Node *t2 = F2.get_component(0);
		//F1->preorder_number();
		int lchild_pre = t1->lchild()->get_preorder_number();
		int rchild_pre = t1->rchild()->get_preorder_number();
		if (!sync_twins(&F1, &F2)) {
			return NULL;
		}
		//	if (t1->lchild()->get_preorder_number() != lchild_pre ||
		//			t1->rchild()->get_preorder_number() != rchild_pre)
		//		return NULL;
		if (F2.get_component(0)->get_children().size() > 2)
			F2.get_component(0)->fixroot();
		// TODO: maybe stop if F2 is too small?
		int pre_separator = t1->lchild()->get_preorder_number();
		int group_1_total;
		int group_2_total;
		if (t1->rchild()->get_preorder_number() > pre_separator) {
			pre_separator = t1->rchild()->get_preorder_number();
			group_1_total = t1->lchild()->find_leaves().size();
			group_2_total = t1->rchild()->find_leaves().size();
		}
		else {
			group_1_total = t1->rchild()->find_leaves().size();
			group_2_total = t1->lchild()->find_leaves().size();
		}
		//	cout << "g1_total: " << group_1_total << endl;
		//	cout << "g2_total: " << group_2_total << endl;
		Node *best_root = t2->lchild();
		find_best_root_hlpr(t2, pre_separator, group_1_total, group_2_total,
				&best_root, best_root_b_acc);
		//	cout << "t1: " << t1->str_subtree() << endl;
		//	cout << "t2: " << t2->str_subtree() << endl;
		//	cout << "best_root: " << best_root->str_subtree() << endl;
		best_root = T2->find_by_prenum(best_root->get_preorder_number());
		//	cout << "T1: " << T1->str_subtree() << endl;
		//	cout << "best_root: " << best_root->str_subtree() << endl;
		//	T2->reroot(best_root);
		//	cout << "T2: " << T2->str_subtree() << endl;
		return best_root;
	}


	__attribute__((always_inline)) inline Node *find_best_root_Inline(
	Node *(*find_best_root)(Node *T1, Node *T2, double *best_root_b_acc),
	Node *T1, Node *T2) {
		double best_root_b_acc = 0;
		Forest f1 = Forest(T1);
		Forest f2 = Forest(T2);
		sync_twins(&f1, &f2);
		Node *new_root =
			find_best_root(f1.get_component(0), f2.get_component(0), &best_root_b_acc);
		if (new_root != NULL)
			new_root = T2->find_by_prenum(new_root->get_preorder_number());
		return new_root;
	}

	__attribute__((always_inline)) inline double find_best_root_acc_Inline(
	Node *(*find_best_root)(Node *T1, Node *T2, double *best_root_b_acc),
	Node *T1, Node *T2) {
		double best_root_b_acc = -1;
		Forest f1 = Forest(T1);
		Forest f2 = Forest(T2);
		sync_twins(&f1, &f2);
		find_best_root(f1.get_component(0), f2.get_component(0), &best_root_b_acc);
		return best_root_b_acc;
	}


	__attribute__((always_inline)) inline void find_best_root_hlpr_Inline(
	void (*find_best_root_hlpr)(Node *n, int pre_separator, int group_1_total, int group_2_total, Node **best_root, double *best_root_b_acc,int *p_group_1_descendants, int *p_group_2_descendants, int *num_ties),
	Node *T2, int pre_separator, int group_1_total,
		int group_2_total, Node **best_root, double *best_root_b_acc) {
		list<Node*>::iterator c;
		int group_1_descendants = 0;
		int group_2_descendants = 0;
		int num_ties = 2;
		for(c = T2->get_children().begin(); c != T2->get_children().end(); c++) {
			find_best_root_hlpr(*c, pre_separator, group_1_total,
					group_2_total, best_root, best_root_b_acc,
					&group_1_descendants, &group_2_descendants, &num_ties);
		}
	}

	__attribute__((always_inline)) inline void find_best_root_hlpr_Inline(
	void (*find_best_root_hlpr)(Node *n, int pre_separator, int group_1_total,
			int group_2_total, Node **best_root, double *best_root_b_acc,
			int *p_group_1_descendants, int *p_group_2_descendants, int *num_ties),
	Node *n, int pre_separator, int group_1_total,
		int group_2_total, Node **best_root, double *best_root_b_acc,
		int *p_group_1_descendants, int *p_group_2_descendants, int *num_ties) {
	list<Node*>::iterator c;
	int group_1_descendants = 0;
	int group_2_descendants = 0;
//	vector<pair<int,int> > children_splits = vector<pair<int,int>>();
	for(c = n->get_children().begin(); c != n->get_children().end(); c++) {
//		int g_1_desc = 0;
//		int g_2_desc = 0;
		find_best_root_hlpr(*c, pre_separator, group_1_total,
				group_2_total, best_root, best_root_b_acc,
				&group_1_descendants, &group_2_descendants, num_ties);
//				&g_1_desc, &g_2_desc, num_ties);
//		children_splits.push_back(make_pair(g_1_desc,g_2_desc));
//		group_1_descendants += g_1_desc;
//		group_2_descendants += g_2_desc;
	}
	if (n->is_leaf()) {
		int pre = n->get_twin()->get_preorder_number();
		if (pre < pre_separator)
			group_1_descendants++;
		else
			group_2_descendants++;
	}
//	else if (n->get_children.size() > 2) {
//		if (group_1_descendants / (double) group_1_total >= group_2_descendants / (double) group_2_total) {
//			sort(children_splits.begin(),children_splits.end(), g_1_comp(group_1_total, group_2_total, true));
//		}
//	}
	// balanced accuracy
	// don't bother averaging since we only directly compare them
	double tpos = group_1_descendants;
	double fpos = group_2_descendants;
	double fneg = (group_1_total - group_1_descendants);
	double tneg = (group_2_total - group_2_descendants);
	// balanced accuracy for this bipartition
	double b_acc =  tpos / (tpos + fneg)
			+ tneg / (tneg + fpos);
	// balanced accuracy for the opposite bipartition
	double b_acc_opp = fpos / (fpos + tneg)
			+ fneg / (fneg + tpos);
	// use the better bipartition
//	cout << "n: " << n->str_subtree() << endl;
//	cout << "b_acc: " << b_acc << endl;
//	cout << "b_acc_opp: " << b_acc_opp << endl;
//	cout << n->str_subtree() << endl;
//	cout << b_acc << "\t" << b_acc_opp << endl;
	if (b_acc_opp > b_acc)
		b_acc = b_acc_opp;
	if (b_acc > *best_root_b_acc) {
		*best_root = n;
		*best_root_b_acc = b_acc;
		*num_ties = 2;
	}
	else if(b_acc == *best_root_b_acc) {
		int r = rand();
		if (r < RAND_MAX/ *num_ties) {
			*best_root = n;
			*best_root_b_acc = b_acc;
		}
	}
	*p_group_1_descendants += group_1_descendants;
	*p_group_2_descendants += group_2_descendants;
}


	__attribute__((always_inline)) inline Node *find_random_root_Inline(Node *T1, Node *T2) {
		vector<Node *> rroots = T2->find_descendants();
		int r = rand() % rroots.size();
		return rroots[r];
	}

	__attribute__((always_inline)) inline Node *find_best_root_rspr_Inline(
		int (*rSPR_branch_and_bound_simple_clustering)(Node *T1, Node *T2),

		Node *T1, Node *T2) {
		Node *t1 = new Node(*T1);
		//	t1->preorder_number();
		Node *t2  = new Node(*T2);
		int new_prenum = T2->lchild()->get_preorder_number();
		vector<Node *> roots = t2->find_descendants();
		vector<int> root_prenums = vector<int>(roots.size());
		for(int i = 0; i < roots.size(); i++) {
			root_prenums[i] = roots[i]->get_preorder_number();
		}
		int best_distance = INT_MAX;
		int num_ties = 2;
		//	cout << endl;
		//	cout << "find_best_root_rspr" << endl;
		//	cout << "T1: " << T1->str_subtree() << endl;
		//	cout << "T2: " << T2->str_subtree() << endl;
		//	cout << roots.size() << " Rootings" << endl;
		for(int i = 0; i < roots.size(); i++) {
			Node *root = roots[i];
			t2->reroot(root);
			//		cout << "\t" << t2->str_subtree() << endl;
			t2->set_depth(0);
			t2->fix_depths();
			t2->preorder_number();
			int distance = rSPR_branch_and_bound_simple_clustering(t1, t2);
			//		int distance = rf_distance(t1, t2);
			if (distance < best_distance) {
				best_distance = distance;
				new_prenum = root_prenums[i];
				num_ties = 2;
			}
			else if (distance == best_distance) {
				int r = rand();
				if (r < RAND_MAX / num_ties) {
					best_distance = distance;
					new_prenum = root_prenums[i];
				}
				num_ties++;
			}
		}
		Node *new_root = T2->find_by_prenum(new_prenum);
		t1->delete_tree();
		t2->delete_tree();
		return new_root;
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
		//SHADOWED FUNCTION
		int (*count_differing_bipartitions)(Node* n),
		// Original param
		Node *n) {
		//cout << "Start: " << n->str_subtree() << endl;
		int count = 0;

		// Using the Iterator, we iterate over a list of nodes

		list<Node *>::iterator c;
		for(c = n->get_children().begin(); c != n->get_children().end(); c++) {
			count += count_differing_bipartitions(*c);
		}
		if (n->get_twin() == NULL || n != n->get_twin()->get_twin()) {
			count++;
		}
		return count;
	}

   __attribute__((always_inline)) inline bool is_nonbranching_Inline(
   	// SHADOWED PARAMS
 	bool &CUT_ONE_B, bool &CUT_TWO_B, bool &CUT_TWO_B_ROOT, bool &REVERSE_CUT_ONE_B, bool &REVERSE_CUT_ONE_B_2,

 	// ORIGINAL PARAMS
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