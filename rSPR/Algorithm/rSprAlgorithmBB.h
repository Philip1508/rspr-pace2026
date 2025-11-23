



namespace rSprBB {





// rSPR_branch_and_bound recursive helper function
__attribute__((always_inline)) inline int rSPR_branch_and_bound_hlpr_Inline(
// SHADOWED ARGUMENTS:
bool &PREFER_RHO,
bool &LEAF_REDUCTION,
bool &PREFER_NONBRANCHING,
bool &DEEPEST_ORDER,
bool &DEEPEST_PROTECTED_ORDER,
bool &CHECK_MERGE_DEPTH,
bool &EDGE_PROTECTION,
bool &CUT_LOST,
bool &CUT_ALL_B,
bool &CUT_ONE_B,
bool &CUT_ONE_AB,
bool &CUT_TWO_B,
bool &REVERSE_CUT_ONE_B,
bool &REVERSE_CUT_ONE_B_3,
bool &REVERSE_CUT_ONE_B_2,
bool &CUT_TWO_B_ROOT,
bool &BB,
bool &CLUSTER_REDUCTION,
int &NUM_CLUSTERS,
int &MAX_CLUSTERS,
bool &EDGE_PROTECTION_TWO_B,
bool &CUT_AC_SEPARATE_COMPONENTS,
bool &ABORT_AT_FIRST_SOLUTION,
bool &ALL_MAFS,

void (*add_sibling_pair)(set<SiblingPair> *sibling_pairs, Node *a, Node *c, UndoMachine *um),
bool (*is_nonbranching)(Forest *T1, Forest *T2, Node *T1_a, Node *T1_c, Node *T2_a, Node *T2_c),
SiblingPair (*pop_sibling_pair_2)(set<SiblingPair> *sibling_pairs, UndoMachine *um),
SiblingPair (*pop_sibling_pair_3)(set<SiblingPair>::iterator s, set<SiblingPair> *sibling_pairs, UndoMachine *um),

int (*rSPR_worse_3_approx_hlpr)(Forest *T1, Forest *T2, list<Node *> *singletons, list<Node *> *sibling_pairs, Forest **F1, Forest **F2, bool save_forests),
int (*rSPR_branch_and_bound_range)(Forest *T1, Forest *T2, int end_k),

int (*rSPR_branch_and_bound_hlpr_1)(Forest *T1, Forest *T2, int k,
set<SiblingPair> *sibling_pairs, list<Node *> *singletons,
bool cut_b_only, list<pair<Forest,Forest> > *AFs,
list<Node *> *protected_stack, int *num_ties, Node *prev_T1_a, Node *prev_T1_c),

	int (*rSPR_branch_and_bound_hlpr_2)(
	Forest *T1,
	Forest *T2,
	int k,
	set<SiblingPair> *sibling_pairs,
	list<Node *> *singletons,
	bool cut_b_only,
	list<pair<Forest, Forest>> *AFs,
	list<Node *> *protected_stack,
	int *num_ties),

	int (*rSPR_branch_and_bound_hlpr_3)(
	Forest *T1,
	Forest *T2,
	int k,
	set<SiblingPair> *sibling_pairs,
	list<Node *> *singletons,
	bool cut_b_only,
	list<pair<Forest, Forest>> *AFs,
	list<Node *> *protected_stack,
	int *num_ties,
	Node *prev_T1_a,
	Node *prev_T1_c),



// ORIGINAL ARGUMENTS
Forest *T1, Forest *T2, int k,
set<SiblingPair> *sibling_pairs, list<Node *> *singletons,
bool cut_b_only, list<pair<Forest,Forest> > *AFs,
list<Node *> *protected_stack, int *num_ties, Node *prev_T1_a, Node *prev_T1_c)  {
	#ifdef DEBUG
	cout << "rSPR_branch_and_bound_hlpr()" << endl;
	cout << "\tT1: ";
	T1->print_components();
	cout << "\tT2: ";
	T2->print_components();
	cout << "K=" << k << endl;
	cout << "sibling pairs:";
	for (set<SiblingPair>::iterator i = sibling_pairs->begin(); i != sibling_pairs->end(); i++) {
cout << "  ";
(*i).a->print_subtree_hlpr();
cout << ",";
(*i).c->print_subtree_hlpr();
	}
	cout << endl;
	cout << "protected_stack:";
	for (list<Node *>::iterator i = protected_stack->begin(); i != protected_stack->end(); i++) {
cout << "  ";
(*i)->print_subtree_hlpr();
	}
	cout << endl;
	#endif

	UndoMachine um = UndoMachine();


	while(!singletons->empty() || !sibling_pairs->empty()) {
		// Case 1 - Remove singletons
		while(!singletons->empty()) {
			Node *T2_a = singletons->back();
			#ifdef DEBUG
				cout << "Case 1" << endl;
				cout << "a " << T2_a->str_subtree() << endl;
			#endif

			singletons->pop_back();
			// find twin in T1
			Node *T1_a = T2_a->get_twin();

			//if (T1_a->get_sibling_pair_status() > 0)
			//	cout << T1_a->get_sibling(sibling_pairs) << endl;
			// if this is in the first component of T_2 then
			// it is not really a singleton.
			Node *T1_a_parent = T1_a->parent();
			if (T1_a_parent == NULL)
				continue;
			bool potential_new_sibling_pair = T1_a_parent->is_sibling_pair();
			if (T2_a == T2->get_component(0)) {
				// TODO: should we do this when it happens?
				if (!T1->contains_rho()) {
					um.add_event(new AddRho(T1));
					um.add_event(new AddRho(T2));
					T1->add_rho();
					T2->add_rho();
					k--;
					#ifdef DEBUG
					cout << "adding p element, k=" << k << endl;
					#endif
				}
			}

			// cut the edge above T1_a
			um.add_event(new CutParent(T1_a));
			T1_a->cut_parent();

			um.add_event(new AddComponent(T1));
			T1->add_component(T1_a);
			ContractEvent(&um, T1_a_parent);


			Node *node = T1_a_parent->contract();

			if (node != NULL && potential_new_sibling_pair && node->is_sibling_pair()){
				add_sibling_pair(sibling_pairs, node->lchild(), node->rchild(),
						&um);
			}
			#ifdef DEBUG
				cout << "\tT1: ";
				T1->print_components();
				cout << "\tT2: ";
				T2->print_components();
			#endif
		//			sp_i = sibling_pairs->begin();
		}
		if(!sibling_pairs->empty()) {
			Node *T1_a;
			Node *T1_c;
			set<SiblingPair>::iterator deepest_valid = sibling_pairs->end();
			int deepest_depth = INT_MAX;
			int deepest_depth_2 = INT_MAX;
			Node *best_a = NULL;
			Node *best_c = NULL;
			/* pop protected_stack when out of order sibling pairs
			 have already contracted it */
			if(!protected_stack->empty()
					&& protected_stack->back()->get_twin()->parent() == NULL
					&& protected_stack->back()->get_twin() != T1->get_component(0)) {
				um.undo_all();
				return -1;
			}
			while(!protected_stack->empty()
					&& (protected_stack->back()->is_contracted()
					// this shouldn't happen
						|| protected_stack->back()->get_twin()->parent() == NULL)) {
				um.add_event(new ListPopBack(protected_stack));
				protected_stack->pop_back();
			}
			if (LEAF_REDUCTION && !cut_b_only) {
				bool found = false;
				set<SiblingPair>::iterator sp_i = sibling_pairs->begin();
				// correct in case sibling pair involves previous
		/*				if (sp_i != sibling_pairs->begin()) {
					if (check_all_pairs)
						sp_i = sibling_pairs->begin();
					else
						sp_i--;
				}
				*/
				while (sp_i != sibling_pairs->end()) {
					T1_a = (*sp_i).a;
					T1_c = (*sp_i).c;
					if (T1_a->parent() == NULL || T1_a->parent() != T1_c->parent()) {
						um.add_event(new RemoveSetSiblingPairs(sibling_pairs,
									SiblingPair(T1_a, T1_c)));
						set<SiblingPair>::iterator rem = sp_i;
						sp_i++;
						sibling_pairs->erase(rem);
						continue;
					}
					Node *T2_a = T1_a->get_twin();
					Node *T2_c = T1_c->get_twin();
					// select if this is a Case 2 or, optionally, nonbranching
					if ((T2_a->parent() != NULL && T2_a->parent() == T2_c->parent())
							|| (!cut_b_only && PREFER_NONBRANCHING
									&& is_nonbranching(T1, T2, T1_a, T1_c, T2_a, T2_c))) {
						um.add_event(new RemoveSetSiblingPairs(sibling_pairs,
									SiblingPair(T1_a, T1_c)));
						set<SiblingPair>::iterator rem = sp_i;
						sp_i++;
						sibling_pairs->erase(rem);
						found = true;
						break;
					}
					if (DEEPEST_ORDER) {
						int depth;
						int depth2;
						if (T1_a->get_depth() < T1_c->get_depth())
							depth = T1_c->get_depth();
						else
							depth = T1_a->get_depth();
						if (T2_a->get_depth() < T2_c->get_depth())
							depth2 = T2_c->get_depth();
						else
							depth2 = T2_a->get_depth();
						if (deepest_valid == sibling_pairs->end()
								|| deepest_depth < depth
								|| (deepest_depth == depth && deepest_depth_2 < depth2) ) {
							// TODO: this crashes on bigtest2
							// Why can we end up cutting the protected node?
							if (!DEEPEST_PROTECTED_ORDER
									|| protected_stack->empty()
									|| (protected_stack->back()->get_twin()->parent()->get_edge_pre_start()
											<= T1_a->get_preorder_number()
										&& protected_stack->back()->get_twin()->parent()->get_edge_pre_end()
											>= T1_a->get_preorder_number())) {
								deepest_valid = sp_i;
								deepest_depth = depth;
								deepest_depth_2 = depth2;
							}
						}
					}
					/* TODO: create a stack of protected nodes (intervals?) and only
						 accept the deepest sibling pair within the interval
					*/

					/* TODO: remember to pop the stack when we include the protected
					   node */
					sp_i++;
				}
				if (!found) {
					if (sibling_pairs->empty())
						continue;
					else {
						SiblingPair spair;
//						cout << "depth: " << deepest_depth << endl;
						if (DEEPEST_ORDER && deepest_valid != sibling_pairs->end())
							spair = pop_sibling_pair_3(deepest_valid, sibling_pairs, &um);
						else
							spair = pop_sibling_pair_2(sibling_pairs, &um);
						T1_a = spair.a;
						T1_c = spair.c;
					}
				}
			}
			else {
				if (prev_T1_a != NULL && prev_T1_c != NULL) {
					T1_a = prev_T1_a;
					T1_c = prev_T1_c;
					prev_T1_a = NULL;
					prev_T1_c = NULL;
				}
				else {
					SiblingPair spair = pop_sibling_pair_2(sibling_pairs, &um);
					T1_a = spair.a;
					T1_c = spair.c;
				}
			}
//			if (T1_a->parent() != NULL)
//				cout << "a_p: " << T1_a->parent()->str_subtree() << endl;
//			if (T1_c->parent() != NULL)
//				cout << "c_p: " << T1_c->parent()->str_subtree() << endl;
			if (T1_a->parent() == NULL || T1_a->parent() != T1_c->parent()) {
				continue;
			}
			if (!T1_a->can_be_sibling() || !T1_c->can_be_sibling()) {
				continue;
			}
			Node *T1_ac = T1_a->parent();
			// lookup in T2 and determine the case
			Node *T2_a = T1_a->get_twin();
			Node *T2_c = T1_c->get_twin();

			if (T2_a->parent() != NULL && T2_a->parent() == T2_c->parent()) {
				#ifdef DEBUG
					cout << "Case 2" << endl;
					T1_ac->print_subtree();
				#endif
				Node *T2_ac = T2_a->parent();

				if (CHECK_MERGE_DEPTH &&
						(T2_a->get_max_merge_depth() > T2_ac->get_depth()
							|| T2_c->get_max_merge_depth() > T2_ac->get_depth())) {
					um.undo_all();
					return -1;
				}

				if (!protected_stack->empty() &&
						(T2_a == protected_stack->back()
						 	|| T2_c == protected_stack->back())) {
					um.add_event(new ListPopBack(protected_stack));
					protected_stack->pop_back();
				}
				// CAN THIS HAPPEN TWICE?
				if (!protected_stack->empty() &&
						(T2_a == protected_stack->back()
						 	|| T2_c == protected_stack->back())) {
					um.add_event(new ListPopBack(protected_stack));
					protected_stack->pop_back();
				}


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
				//T1->add_deleted_node(T1_a);
				//T1->add_deleted_node(T1_c);
				//T2->add_deleted_node(T2_a);
				//T2->add_deleted_node(T2_c);

				// check if T2_ac is a singleton
				if (T2_ac->is_singleton() && !T1_ac->is_singleton() && T2_ac != T2->get_component(0))
					singletons->push_back(T2_ac);
				// check if T1_ac is part of a sibling pair
				if (T1_ac->parent() != NULL && T1_ac->parent()->is_sibling_pair()) {
				add_sibling_pair(sibling_pairs, T1_ac->parent()->lchild(), T1_ac->parent()->rchild(),
						&um);
				}
				#ifdef DEBUG
					cout << "\tT1: ";
					T1->print_components();
					cout << "\tT2: ";
					T2->print_components();
				#endif
			}
			/* need to copy trees and lists for branching
			 * use forest copy constructor for T1 and T2 giving T1' and T2'
			 * T1' twins are in T2, and same for T2' and T1.
			 * singleton list will be empty except for maybe above the cut,
			 * so this can be created.
			 * fix one set of twins (T2->T1' or T1->T2' not sure)
			 * exploit chained twin relationship to copy sibling pair list
			 * fix other set of twins
			 * swap T2 and T2' root nodes
			 * now do the cut
			 *
			 * note: don't copy for 3rd cut, is a waste
			 */

			// Case 3
			// note: guaranteed that singleton list is empty
			else {
				if (k <= 0) {
					if ((!CUT_LOST || k < 0 ||
								(T1_a->num_lost_children() == 0 &&
								 T1_c->num_lost_children() == 0))
							&& ((T2_c->parent() != NULL && T2_a->parent() != NULL)|| !T2->contains_rho())) {
						singletons->clear();
						um.undo_all();
						return k-1;
					}
				}
				Forest *best_T1;
				Forest *best_T2;
				int best_k = -1;
				int answer_a = -1;
				int answer_b = -1;
				int answer_c = -1;
				bool cut_ab_only = false;
				bool cut_a_only = false;
				bool cut_c_only = false;
				bool cut_a_or_merge_ac = false;
				bool same_component = true;
				int lca_depth = -1;
				int path_length = -1;
				bool cut_b_only_if_not_a_or_c = false;
				bool cob = false;
				int undo_state = um.num_events();
				//  ensure T2_a is below T2_c
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
				Node *T2_b = T2_a->parent()->rchild();
				if (T2_b == T2_a)
					T2_b = T2_a->parent()->lchild();
				bool multi_node = false;
				if (T2_a->parent()->get_children().size() > 2)
					multi_node = true;

			if (CUT_ONE_B) {
				if (T2_a->parent()->parent() == T2_c->parent()
					&& T2_c->parent() != NULL && !cut_b_only)
					cut_b_only=true;
					cob = true;
			}
			else if (CUT_ONE_AB) {
				if (T2_a->parent()->parent() == T2_c->parent()
					&& T2_c->parent() != NULL)
					cut_ab_only=true;
			}
			if (CUT_TWO_B && !cut_b_only && T1_ac->parent() != NULL) {
				Node *T1_s = T1_ac->get_sibling();
				if (T1_s->is_leaf()) {
					Node *T2_l = T2_a->parent()->parent();
					// Note: is this too harsh? If T2_l is nonbinary then can we do cut_b_only_if_not_a_or_c ?
					if (T2_l != NULL && T2_l->get_children().size() <= 2) {
						if (T2_c->parent() != NULL && T2_c->parent()->parent() == T2_l
								&& ((T2_a->parent()->get_children().size() <= 2
										&& T2_c->parent()->get_children().size() <= 2)
									|| T1_s->get_twin()->is_protected())) {
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
								&& ((T2_a->parent()->get_children().size() <= 2
										&& T2_a->parent()->parent()->get_children().size() <= 2
										&& T2_l->get_children().size() <= 2)
									|| T1_s->get_twin()->is_protected())) {
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
			if (REVERSE_CUT_ONE_B && (!cut_b_only || (cob && multi_node)) &&
					T1_ac->parent() != NULL) {
				Node *T1_s = T1_ac->get_sibling();
				if (T1_s->is_leaf()) {
					Node *T2_s = T1_s->get_twin();
					if (T2_s->parent() == T2_a->parent()) {
						cut_c_only=true;
						cut_b_only=false;
						cob=false;
					}
					else if (T2_s->parent() == T2_c->parent()) {
						if (T2_c->parent()->get_children().size() <= 2) {
							cut_a_only=true;
							cut_b_only=false;
							cob=false;
						}
						else {
							cut_a_or_merge_ac=true;
							cut_b_only=false;
							cob=false;
						}
					}
					else if (REVERSE_CUT_ONE_B_3
							// TODO: there is a chance for an additional optimization
							// here. If T2_s is not protected then we can cut c or
							// have to cut a (and s?)
//							&& T2_c->parent()->parent() != NULL
							// TODO: buggy? Do we also have to consider cutting b_1 through the other b's?
							&& T2_s->is_protected()
							&& T2_s->parent() != NULL
							&& T2_s->parent()->parent() == T2_a->parent()
							&& T2_s->parent()->get_children().size() <= 2) {
						//cut_c_only = true;
						cut_b_only=false;
						cob=false;
						if (!T2_a->is_protected()) {
							um.add_event(new ProtectEdge(T2_a));
							T2_a->protect_edge();
						}
					}
				}
				else if (REVERSE_CUT_ONE_B_2 && T2_c->parent() != NULL
						&& rSprUtility::chain_match(T1_s, T2_c->get_sibling(), T2_a)) {
					cut_a_only = true;
					cut_b_only=false;
					cob=false;
				}
			}
			if (REVERSE_CUT_ONE_B_3 && (!cut_b_only || (cob && multi_node)) &&
					T1_ac->parent() != NULL && T1_ac->parent()->parent() != NULL) {
				Node *T1_s = T1_ac->parent()->get_sibling();
				Node *T2_s = T1_s->get_sibling();
				if (T1_s->is_leaf()
						&& T2_s->is_protected()
						&& T2_s->parent() != NULL
						&& T2_s->parent()->parent() == T2_a->parent()
						&& T2_s->parent()->get_children().size() <= 2) {
					//cut_c_only = true;
					cut_b_only=false;
					cob=false;
					if (!T2_a->is_protected()) {
						um.add_event(new ProtectEdge(T2_a));
						T2_a->protect_edge();
					}
				}
			}
			if (CUT_TWO_B_ROOT && cut_a_only == false && cut_c_only == false
					&& cut_b_only_if_not_a_or_c == true) {
				cut_b_only = true;
			}
/*			if (CUT_LOST) {
				if (T1_a->num_lost_children() > 0
						|| T2_a->num_lost_children() > 0) {
					cut_a_only = true;
					cut_b_only = false;
					cut_c_only = false;
					k++;
				}
				else if (T1_c->num_lost_children() > 0
						|| T2_c->num_lost_children() > 0) {
					cut_a_only = false;
					cut_b_only = false;
					cut_c_only = true;
					k++;
				}
				else if (T2_b->is_leaf()) {
					if (T2_b->num_lost_children() > 0
						|| T2_b->get_twin()->num_lost_children() > 0) {
					cut_a_only = false;
					cut_b_only = true;
					cut_c_only = false;
					k++;
					}
				}
			}
*/

			#ifdef DEBUG
					cout << "Case 3" << endl;
					cout << "\tT1: ";
					T1->print_components();
					cout << "\tT2: ";
					T2->print_components();
					cout << "\tK=" << k << endl;
					cout << "\tsibling pairs:";
					for (set<SiblingPair>::iterator i = sibling_pairs->begin(); i != sibling_pairs->end(); i++) {
						cout << "  ";
						(*i).a->print_subtree_hlpr();
						cout << ",";
						(*i).c->print_subtree_hlpr();
					}
					cout << endl;
					cout << "\tprotected_stack:";
					for (list<Node *>::iterator i = protected_stack->begin(); i != protected_stack->end(); i++) {
						cout << "  ";
						(*i)->print_subtree_hlpr();
					}
					cout << endl;
					cout << "\tcut_a_only=" << cut_a_only << endl;
					cout << "\tcut_b_only=" << cut_b_only << endl;
					cout << "\tcut_c_only=" << cut_c_only << endl;
					cout << "\tT2_a " << T2_a->str() << " "
						<< T2_a->get_depth() << endl;
					cout << "\tT2_c " << T2_c->str() << " "
						<< T2_c->get_depth() << endl;
					cout << "\tT2_b " << T2_b->str_subtree() << " "
						<< T2_b->get_depth() << endl;
				#endif



				// copy elements
					/*
				Forest *T1_copy;
				Forest *T2_copy;
				list<Node *> *sibling_pairs_copy;
				Node *T1_a_copy;
				Node *T1_c_copy;
				Node *T2_a_copy;
				Node *T2_c_copy;
				*/
				//list<Node *> *singletons_copy = new list<Node *>();

				// make copies for the approx
				// be careful we do not kill real T1 and T2
				// ie use the copies
				if (BB && !cut_a_only && !cut_b_only && !cut_c_only) {
					list<Node *> *spairs;
					spairs = new list<Node *>();
					spairs->push_back(T1_c);
					spairs->push_back(T1_a);
					for (set<SiblingPair>::iterator i = sibling_pairs->begin(); i != sibling_pairs->end(); i++) {
						spairs->push_back((*i).a);
						spairs->push_back((*i).c);
					}
					int approx_spr = rSPR_worse_3_approx_hlpr(T1, T2,
							singletons, spairs, NULL, NULL, false);
					delete spairs;
					#ifdef DEBUG
						cout << "\tT1: ";
						T1->print_components();
						cout << "\tT2: ";
						T2->print_components();
						cout << "approx =" << approx_spr << endl;
					#endif
					if (approx_spr  >  3*k){
						#ifdef DEBUG
							cout << "approx failed" << endl;
						#endif
						um.undo_all();
						return -1;
					}
					um.undo_to(undo_state);
				}

				if (!cut_a_only && !cut_c_only && !cut_b_only)
					same_component = T2_a->same_component(T2_c, lca_depth, path_length);
				if (cut_b_only)
					same_component = true;

			if (CLUSTER_REDUCTION && (MAX_CLUSTERS < 0 || NUM_CLUSTERS < MAX_CLUSTERS)) {
				// clean up singletons
				// TODO: this is duplication
				/*
			while(!singletons->empty()) {
				Node *T2_a = singletons->back();
				singletons->pop_back();
				// find twin in T1
				Node *T1_a = T2_a->get_twin();
				// if this is in the first component of T_2 then
				// it is not really a singleton.
				Node *T1_a_parent = T1_a->parent();
				if (T1_a_parent == NULL)
					continue;
				bool potential_new_sibling_pair = T1_a_parent->is_sibling_pair();
				if (T2_a == T2->get_component(0)) {
					T1->add_rho();
					T2->add_rho();
					k--;
				}

				// cut the edge above T1_a
				T1_a->cut_parent();
				T1->add_component(T1_a);
				Node *node = T1_a_parent->contract();
				if (node != NULL && potential_new_sibling_pair && node->is_sibling_pair()){
					sibling_pairs->push_front(node->lchild());
					sibling_pairs->push_front(node->rchild());
				}
			}
			*/
	//			cout << "foo" << endl;
	//			cout << "foo2" << endl;
//				cout << "\tT1: ";
//				T1->print_components();
//				cout << "\tT2: ";
//				T2->print_components();
				sync_interior_twins_real(T1, T2);
				list<Node *> *cluster_points = find_cluster_points(T1, T2);
				//cluster_points->erase(++cluster_points->begin(),cluster_points->end());

				// TODO: could this be faster by using the approx to allocate
				// a certain amount of the k to different clusters?
				// TODO: write pseudocode for what we need
				// TODO: then implement it
				// NOTE: need to make a list of ClusterInstances and then
				// solve each.
				// TODO: where should we do this? Just before we would
				// normally branch?

//				cout << "k=" << k << endl;
//				cout << "cp=" << cluster_points->size() << endl;
				if (!cluster_points->empty()) {
					NUM_CLUSTERS++;
					sibling_pairs->clear();
#ifdef DEBUG_CLUSTERS
					cout << "CLUSTERS" << endl;
					for(int j = 0; j < 70; j++) {
						cout << "*";
					}
					cout << endl;
					for(list<Node *>::iterator i = cluster_points->begin();
							i != cluster_points->end(); i++) {
						cout << (*i)->str_subtree() << endl;
						cout << (*i)->get_twin()->str_subtree() << endl;
						for(int j = 0; j < 70; j++) {
							cout << "*";
						}
						cout << endl;
					}
					cout << endl;
#endif

					list<ClusterInstance> clusters =
						cluster_reduction(T1, T2, cluster_points);

					// TODO: make it so we don't need this?
					T1->unsync_interior();
					T2->unsync_interior();
					while(!clusters.empty()) {
						ClusterInstance cluster = clusters.front();
						clusters.pop_front();
						cluster.F1->unsync_interior();
						cluster.F2->unsync_interior();
#ifdef DEBUG_CLUSTERS
						cout << "CLUSTER_START" << endl;
						cout << &(*cluster.F1->get_component(0)) << endl;
						if (!clusters.empty())
							cout << &(*clusters.front().F1->get_component(0)) << endl;
						cout << "\tF1: ";
						cluster.F1->print_components();
						cout << "\tF2: ";
						cluster.F2->print_components();
						cout << "K=" << k << endl;
#endif
						int cluster_spr = -1;
						if (k >= 0) {
							// hack for clusters with no rho
//							cout << __LINE__ << endl;
//							cout << cluster.F2_cluster_node << endl;
//							cout << cluster.F2_has_component_zero << endl;
							if ((cluster.F2_cluster_node == NULL
										|| (cluster.F2_cluster_node->is_leaf()
												&& cluster.F2_cluster_node->parent() == NULL
												&& cluster.F2_cluster_node->
												get_num_clustered_children() <= 1
												&& (cluster.F2_cluster_node !=
														cluster.F2_cluster_node->get_forest()->
															get_component(0))))
										&& cluster.F2_has_component_zero == false) {
//							cout << __LINE__ << endl;
								cluster.F1->add_rho();
								cluster.F2->add_rho();
							}
							cluster_spr = rSPR_branch_and_bound_range(cluster.F1,
									cluster.F2, k);
							if (cluster_spr >= 0) {
//							cout << "cluster k=" << cluster_spr << endl;
//							cout << "\tF1: ";
//							cluster.F1->print_components();
//							cout << "\tF2: ";
//							cluster.F2->print_components();
								k -= cluster_spr;
							}
							else {
								k = -1;
							}
						}
						if (k > -1) {
//							cout << "\tF1: ";
//							T1->print_components();
//							cout << "\tF2: ";
//							T2->print_components();
								if (!cluster.is_original()) {
									int adjustment = cluster.join_cluster(T1, T2);
									k -= adjustment;
									delete cluster.F1;
									delete cluster.F2;

	//						cout << cluster.F1_cluster_node->str_subtree() << endl;
	//						cout << cluster.F2_cluster_node->str_subtree() << endl;
	//						Node *p = cluster.F1_cluster_node;
	//						while (p->parent() != NULL)
	//							p = p->parent();
	//						cout << &(*p) << endl;
	//						cout << p->str_subtree() << endl;


							//cout << "\tF1: ";
							//T1->print_components();
							//cout << "\tF2: ";
							//T2->print_components();
								}
//						else {
//							cout << "original" << endl;
//						}
							}
							else {
								if (!cluster.is_original()) {
									//if (cluster.F1_cluster_node != NULL)
									//	cluster.F1_cluster_node->contract();
									//if (cluster.F2_cluster_node != NULL)
									//	cluster.F2_cluster_node->contract();
									delete cluster.F1;
									delete cluster.F2;
								}
							}
					}
					delete cluster_points;
//					cout << "returning k=" << k << endl;
					NUM_CLUSTERS--;
					return k;
				}
				else {
					T1->unsync_interior();
					T2->unsync_interior();
					// HACK to allow only initial clusters
					// TODO: use UndoMachine in this clustering section
					// and update this clustering to not require copying
					NUM_CLUSTERS++;
				}
				delete cluster_points;

	//			cout << "done" << endl;
			}
				 // make copies for the branching
			/*
				copy_trees(&T1, &T2, &sibling_pairs, &T1_a, &T1_c, &T2_a, &T2_c,
						&T1_copy, &T2_copy, &sibling_pairs_copy,
					&T1_a_copy, &T1_c_copy, &T2_a_copy, &T2_c_copy);
					*/

				Node *node;
				Node *T2_ab = T2_a->parent();

				bool balanced = false;
				bool multi_b1 = T2_ab->get_children().size() > 2;
				bool multi_b2 = false;
				Node *T2_d = NULL;
				if (path_length == 4) {
					if (T2_a->parent()->parent() == T2_c->parent()->parent())
						balanced = true;
					if (balanced)
						multi_b2 = T2_c->parent()->get_children().size() > 2;
					else
						multi_b2 = T2_ab->parent()->get_children().size() > 2;
					if (balanced)
						T2_d = T2_c->get_sibling();
				}

				// cut T2_a
//				if (cut_b_only == false && T2_a->is_protected())
//					cout << "protected k=" << k << endl;
				if (cut_b_only == false && cut_c_only == false &&
						!T2_a->is_protected()
						&& (T2_a->parent()->parent() != NULL
								|| (T2_a->parent() == T2->get_component(0)
										&& !T2->contains_rho())
								|| !T2_b->is_protected()
								|| T2_a->parent()->get_children().size() > 2)) {// &&
//						(!T2_a->parent()->is_protected() ||
//							T2_a->parent()->get_children().size() > 2)) { }
					um.add_event(new CutParent(T2_a));
					T2_a->cut_parent();
					ContractEvent(&um, T2_ab);
					node = T2_ab->contract();
					if (node != NULL && node->is_singleton() &&
							node != T2->get_component(0))
						singletons->push_back(node);
					um.add_event(new AddComponent(T2));
					T2->add_component(T2_a);
					singletons->push_back(T2_a);

					// also require !cut_a_only ?
//					if (EDGE_PROTECTION_TWO_B && T2_c->is_protected() && !cut_a_only) {
//				}

					if (EDGE_PROTECTION_TWO_B && T2_c->is_protected() && !cut_a_only){
						if (path_length == 4) {
							if (!multi_b1 && !multi_b2 && !T2_b->is_protected()) {
								um.add_event(new ProtectEdge(T2_b));
								T2_b->protect_edge();
							}
							if (!multi_b2 && !multi_b1) {
								Node *T2_b2 = T2_b->parent()->get_sibling();
								if (balanced)
									T2_b2 = T2_d;
								if (!T2_b2->is_protected()) {
									um.add_event(new ProtectEdge(T2_b2));
									T2_b2->protect_edge();
								}
						}
						}
					}
					if (cut_a_only) {
						answer_a =
							rSPR_branch_and_bound_hlpr_1(T1, T2, k-1, sibling_pairs,
									singletons, false, AFs, protected_stack, num_ties, T1_c, T1_c->get_sibling());
					}
					else {
						answer_a =
							rSPR_branch_and_bound_hlpr_2(T1, T2, k-1, sibling_pairs,
									singletons, false, AFs, protected_stack, num_ties);
					}
				}
				best_k = answer_a;
				best_T1 = T1;
				best_T2 = T2;

				um.undo_to(undo_state);
/*				#ifdef DEBUG
					cout << "Case 3 CHECK" << endl;
					cout << "\tT1: ";
					T1->print_components();
					cout << "\tT2: ";
					T2->print_components();
					cout << "\tK=" << k << endl;
					cout << "\tsibling pairs:";
					for (list<Node *>::iterator i = sibling_pairs->begin(); i != sibling_pairs->end(); i++) {
						cout << "  ";
						(*i)->print_subtree_hlpr();
					}
					cout << endl;
					cout << "\tcut_b_only=" << cut_b_only << endl;
					cout << "\tT2_a " << T2_a->str() << " "
						<< T2_a->get_depth() << endl;
					cout << "\tT2_c " << T2_c->str() << " "
						<< T2_c->get_depth() << endl;
					cout << "\tT2_b " << T2_b->str_subtree() << " "
						<< T2_b->get_depth() << endl;
				#endif
*/

				//load the copy
				/*
				T1 = T1_copy;
				T2 = T2_copy;
				T1_a = T1_a_copy;
				T1_c = T1_c_copy;
				T2_a = T2_a_copy;
				T2_c = T2_c_copy;
				sibling_pairs = sibling_pairs_copy;
				singletons = new list<Node *>();
				*/


				// make copies for the branching
				/*
				copy_trees(&T1, &T2, &sibling_pairs, &T1_a, &T1_c, &T2_a, &T2_c,
						&T1_copy, &T2_copy, &sibling_pairs_copy,
						&T1_a_copy, &T1_c_copy, &T2_a_copy, &T2_c_copy);
						*/


				// get T2_b
				T2_ab = T2_a->parent();
				T2_b = T2_ab->rchild();

				if (T2_b == T2_a)
					T2_b = T2_ab->lchild();


//				if ((!CUT_AC_SEPARATE_COMPONENTS ||
//							T2_a->find_root() == T2_c->find_root())
//						&& (!multi_node && T2_b->is_protected()))
//					cout << "protected k=" << k << endl;

				// BUG: we need to cut b or c when b is a multifurcating COB
				// TODO: if the LCA of B_1's leaves in T1 is not an ancestor of
				// a then this is safe
				if (multi_node && cob) {
					vector<Node *> B_1 = T2_a->parent()->find_leaves();
					LCA T1_LCA = LCA(T1->get_component(0));
					Node *B_1_lca = NULL;
					for(int i = 0; i < B_1.size(); i++) {
						if (B_1[i] == T2_a) {
							continue;
						}
						if (B_1_lca == NULL) {
							B_1_lca = B_1[i]->get_twin();
						}
						else {
							B_1_lca = T1_LCA.get_lca(B_1_lca, B_1[i]->get_twin());
						}
					}
					Node *T1_a_ancestor = T1_a;
					while(T1_a_ancestor != NULL && T1_a_ancestor != B_1_lca) {
						T1_a_ancestor = T1_a_ancestor->parent();
					}
					if (T1_a_ancestor != NULL) {
						#ifdef DEBUG
							cout << "\tmultifurcating cut_b_only so will also cut c" << endl;
						#endif
						cut_b_only=false;
					}
				}

				// cut T2_b
				if ((!CUT_AC_SEPARATE_COMPONENTS || same_component)
//						&& ((!T2_b->parent()->is_protected()
								&& (((multi_node || !T2_b->is_protected())))
						&& (!ABORT_AT_FIRST_SOLUTION || best_k < 0
							|| !PREFER_RHO || !AFs->front().first.contains_rho() )
						&& !cut_a_only && !cut_c_only
						&& (T2_a->parent()->parent() != NULL
								|| !T2_a->is_protected()
								|| (T2_a->parent() == T2->get_component(0)
										&& !T2->contains_rho()))) {
					if (multi_node) {
						um.add_event(new ChangeEdgePreInterval(T2_a));
						T2_a->copy_edge_pre_interval(T2_ab);
						um.add_event(new CutParent(T2_a));
						T2_a->cut_parent();
						um.add_event(new ChangeEdgePreInterval(T2_ab));
						T2_ab->set_edge_pre_start(-1);
						T2_ab->set_edge_pre_end(-1);
						Node *T2_ab_parent = T2_ab->parent();
						if (T2_ab_parent != NULL) {
							um.add_event(new CutParent(T2_ab));
							T2_ab->cut_parent();
							um.add_event(new AddChild(T2_a));
							T2_ab_parent->add_child(T2_a);
							um.add_event(new AddComponent(T2));
							T2->add_component(T2_ab);
						}
						else {
							if (T2->get_component(0) == T2_ab) {
								um.add_event(new AddComponentToFront(T2));
								T2->add_component(0, T2_a);
							}
							else {
								um.add_event(new AddComponent(T2));
								T2->add_component(T2_a);
								singletons->push_back(T2_a);
							}
						}
					}
					else {
						um.add_event(new CutParent(T2_b));
						T2_b->cut_parent();
						ContractEvent(&um, T2_ab);
						node = T2_ab->contract();
						if (node != NULL && node->is_singleton()
								&& node != T2->get_component(0))
								singletons->push_back(node);
						um.add_event(new AddComponent(T2));
						T2->add_component(T2_b);
						if (T2_b->is_leaf())
							singletons->push_back(T2_b);
					}
				add_sibling_pair(sibling_pairs, T1_a, T1_c,
						&um);

					// TODO: check carefully

					if (cut_a_or_merge_ac) {
						if (!T2_a->is_protected()) {
							um.add_event(new ProtectEdge(T2_a));
							T2_a->protect_edge();
							um.add_event(new ListPushBack(protected_stack));
							protected_stack->push_back(T2_a);
						}
						if (!T2_c->is_protected()) {
							um.add_event(new ProtectEdge(T2_c));
							T2_c->protect_edge();
						}
					}

					if (CUT_ALL_B) {
						answer_b =
							rSPR_branch_and_bound_hlpr_3(T1, T2, k-1,
									sibling_pairs, singletons, true, AFs, protected_stack,
									num_ties, T1_a, T1_c);
					}
					else {
						answer_b =
							rSPR_branch_and_bound_hlpr_3(T1, T2, k-1,
									sibling_pairs, singletons, false, AFs, protected_stack,
									num_ties, T1_a, T1_c);
					}
				}
				if (answer_b > best_k
						|| (answer_b == best_k
							&& PREFER_RHO
							&& T2->contains_rho() )) {
					best_k = answer_b;
					//swap(&best_T1, &T1);
					//swap(&best_T2, &T2);
				}

				um.undo_to(undo_state);

				/*
				delete T1;
				delete T2;
				delete sibling_pairs;
				delete singletons;
				*/


				// load the copy
				/*
				T1 = T1_copy;
				T2 = T2_copy;
				T1_a = T1_a_copy;
				T1_c = T1_c_copy;
				T2_a = T2_a_copy;
				T2_c = T2_c_copy;
				sibling_pairs = sibling_pairs_copy;
				singletons = new list<Node *>();
				*/
//				if (T2_c->is_protected())
//					cout << "protected k=" << k << endl;
				if (!T2_c->is_protected() &&
						!cut_a_or_merge_ac &&
	//					(T2_c->parent() == NULL || !T2_c->parent()->is_protected() ||
	//						T2_c->parent()->get_children().size() > 2) &&
						(!ABORT_AT_FIRST_SOLUTION || best_k < 0
							|| !PREFER_RHO || !AFs->front().first.contains_rho() )
						&& cut_b_only == false && cut_ab_only == false
						&& cut_a_only == false
						// TODO: do we allow this if T2_c has no parent?
						// it has to be under rho, right?
						&& (T2_c->parent() == NULL
								|| T2_c->parent()->parent() != NULL
								|| (T2_c->parent() == T2->get_component(0)
										&& !T2->contains_rho())
								|| !T2_c->get_sibling()->is_protected()
								|| T2_c->parent()->get_children().size() > 2)) {// &&


					if (T2_c->parent() != NULL) {
						Node *T2_c_parent = T2_c->parent();
						um.add_event(new CutParent(T2_c));
						T2_c->cut_parent();
						ContractEvent(&um, T2_c_parent);
						node = T2_c_parent->contract();
						if (node != NULL && node->is_singleton()
								&& node != T2->get_component(0))
							singletons->push_back(node);
						um.add_event(new AddComponent(T2));
						T2->add_component(T2_c);
					}
					else {
						// don't decrease k
						k++;
					}
					if (EDGE_PROTECTION && !cut_c_only) {
						if (!T2_a->is_protected()) {
							um.add_event(new ProtectEdge(T2_a));
							T2_a->protect_edge();
//							if (DEEPEST_PROTECTED_ORDER && !cut_c_only) {
							if (DEEPEST_PROTECTED_ORDER) {
								um.add_event(new ListPushBack(protected_stack));
								protected_stack->push_back(T2_a);
							}
							// TODO: add to protected list
						}
//						if (EDGE_PROTECTION_TWO_B && !cut_c_only) {
//					}
						// TODO: problem here :(
						if (EDGE_PROTECTION_TWO_B) {
							if (path_length == 4) {
								if (!multi_b1 && !multi_b2 && !T2_b->is_protected()) {
									um.add_event(new ProtectEdge(T2_b));
									T2_b->protect_edge();
								}
								if (!multi_b2 && !multi_b1) {
									Node *T2_b2 = T2_b->parent()->get_sibling();
									if (balanced)
										T2_b2 = T2_d;
									if (!T2_b2->is_protected()) {
										um.add_event(new ProtectEdge(T2_b2));
										T2_b2->protect_edge();
									}
								}
							}
						}
						if (path_length == 5)
							T2_a->set_max_merge_depth(lca_depth);
					}
						singletons->push_back(T2_c);
						if (cut_c_only) {
							answer_c =
								rSPR_branch_and_bound_hlpr_1(T1, T2, k-1, sibling_pairs,
										singletons, false, AFs, protected_stack, num_ties, T1_a, T1_a->get_sibling());
						}
						else {
							answer_c =
								rSPR_branch_and_bound_hlpr_2(T1, T2, k-1, sibling_pairs,
										singletons, false, AFs, protected_stack, num_ties);
						}
						if (answer_c > best_k
									|| (answer_c == best_k
									&& PREFER_RHO
									&& T2->contains_rho() )) {
							best_k = answer_c;
							//swap(&best_T1, &T1);
							//swap(&best_T2, &T2);
						}
				}
				/*
				delete T1;
				delete T2;
				delete sibling_pairs;
				delete singletons;
				*/

				um.undo_to(undo_state);

				//T1 = best_T1;
				//T2 = best_T2;

#ifdef DEBUG_UNDO
		 while(um.num_events() > 0) {
				cout << "Undo step " << um.num_events() << endl;
				cout << "T1: ";
				T1->print_components();
				cout << "T2: ";
				T2->print_components();
					cout << "sibling pairs:";
					for (set<SiblingPair>::iterator i = sibling_pairs->begin(); i != sibling_pairs->end(); i++) {
						cout << "  ";
						(*i).a->print_subtree_hlpr();
						cout << ",";
						(*i).c->print_subtree_hlpr();
					}
					cout << endl;
			 um.undo();
			cout << endl;
		 }
#else
		 um.undo_all();
#endif
				singletons->clear();
				return best_k;
			}
			cut_b_only = false;
		}
	}

	if (k >= 0) {
		if (PREFER_RHO && !AFs->empty() && !AFs->front().first.contains_rho() && T1->contains_rho()) {
			if (!ALL_MAFS)
				AFs->clear();
			AFs->push_front(make_pair(Forest(T1),Forest(T2)));
			*num_ties = 2;
		}
		else if (ALL_MAFS || AFs->empty()) {
			AFs->push_back(make_pair(Forest(T1),Forest(T2)));
		}
		else if (!PREFER_RHO || AFs->front().first.contains_rho() == T1->contains_rho()) {
			if (rand() < RAND_MAX/ *num_ties) {
				AFs->clear();
				AFs->push_back(make_pair(Forest(T1),Forest(T2)));
			}
			(*num_ties)++;
		}
	}

#ifdef DEBUG_UNDO
		 while(um.num_events() > 0) {
				cout << "Undo step " << um.num_events() << endl;
				cout << "T1: ";
				T1->print_components();
				cout << "T2: ";
				T2->print_components();
					cout << "sibling pairs:";
					for (set<SiblingPair>::iterator i = sibling_pairs->begin(); i != sibling_pairs->end(); i++) {
						cout << "  ";
						(*i).a->print_subtree_hlpr();
						cout << ",";
						(*i).c->print_subtree_hlpr();
					}
					cout << endl;
			 um.undo();
			cout << endl;
		 }
#else
		 um.undo_all();
#endif

	return k;
}




	// ---------------------------------------------------------------
	// SEPARATOR
	// ---------------------------------------------------------------



	__attribute__((always_inline)) inline int rSPR_branch_and_bound_simple_clustering_Inline(
	// SHADOWED ARGUMENTS
	int &MAX_SPR,
	bool &PREFER_RHO,
	bool &MULTIFURCATING,
	int &CLUSTER_TUNE,
	bool &CLAMP,

	int (*rSPR_worse_3_mult_approx)(Forest *T1, Forest *T2),
	int (*rSPR_worse_3_approx_2)(Forest *T1, Forest *T2),

	int (*rSPR_branch_and_bound_mult_range)(Forest *T1, Forest *T2, int start_k, int end_k),
	int (*rSPR_branch_and_bound_range)(Forest *T1, Forest *T2, int start_k, int end_k),


	// ORIGINAL ARGUMENTS
	Forest *T1, Forest *T2, bool verbose, map<string, int> *label_map, map<int, string> *reverse_label_map) {
	Forest F1 = *T1;//Forest(T1);
	Forest F2 = *T2;//Forest(T2);
	Forest F3 = Forest(F1);
	Forest F4 = Forest(F2);

	bool do_cluster = true;

//	bool old_rho = PREFER_RHO;
	PREFER_RHO = true;
	if (verbose) {
		cout << "T1: ";
		F1.print_components();
		cout << "T2: ";
		F2.print_components();
	}

	int full_approx_spr;
	if (MULTIFURCATING) {
	  full_approx_spr = rSPR_worse_3_mult_approx(&F3, &F4);
	}
	else {
	  full_approx_spr = rSPR_worse_3_approx_2(&F3, &F4);
	}
	if (full_approx_spr <= CLUSTER_TUNE) {
		do_cluster = false;
	}
	if (verbose) {

		cout << "approx F1: ";
		F3.print_components();
		cout << "approx F2: ";
		F4.print_components();
		// what the AF shows
		cout << "approx drSPR=" << F4.num_components()-1 << endl;
		/* what we use to get the lower bound: 3 * the number of cutting rounds in
			 the approx algorithm
		*/
		//cout << "approx drSPR=" << full_approx_spr << endl;
		cout << "\n";
	}

	if (!sync_twins(&F1, &F2))
		return 0;
	if (F1.get_component(0)->is_leaf())
		return 0;
	list<Node *> *cluster_points;
	sync_interior_twins_real(&F1, &F2);
	if (do_cluster) {
		cluster_points = find_cluster_points(&F1, &F2);
	}

		int total_k = 0;

		if (do_cluster && !cluster_points->empty()) {

			list<ClusterInstance> clusters =
				cluster_reduction(&F1, &F2, cluster_points);

			F1.unsync_interior();
			F2.unsync_interior();
		int k = 0;
		int i = 0;
		while(!clusters.empty()) {
			i++;
			ClusterInstance cluster = clusters.front();
			clusters.pop_front();
			cluster.F1->unsync_interior();
			cluster.F2->unsync_interior();
			if (verbose) {
				cout << "C" << i << "_1: ";
				cluster.F1->print_components();
				cout << "C" << i << "_2: ";
				cluster.F2->print_components();
			}
			Forest f1 = Forest(cluster.F1);
			Forest f2 = Forest(cluster.F2);

			int min_spr;
			if (MULTIFURCATING) {
			  min_spr = rSPR_worse_3_mult_approx(&f1, &f2);
			}
			else {
			  min_spr = rSPR_worse_3_approx_2(&f1, &f2);
			}
			min_spr /= 3;

			if (verbose) {
				cout << "cluster approx drSPR=" << f1.num_components()-1 << endl;
				//cout << "cluster approx drSPR=" << min_spr << endl;
				cout << endl;
			}

			int cluster_spr = -1;
			k = MAX_SPR - total_k;
			if (k >= 0) {
				// hack for clusters with no rho
				if ((cluster.F2_cluster_node == NULL
							|| (cluster.F2_cluster_node->is_leaf()
									&& cluster.F2_cluster_node->parent() == NULL
									&& cluster.F2_cluster_node->
									get_num_clustered_children() <= 1
									&& (cluster.F2_cluster_node !=
											cluster.F2_cluster_node->get_forest()->
												get_component(0))))
							&& cluster.F2_has_component_zero == false) {
					cluster.F1->add_rho();
					cluster.F2->add_rho();
				}
				if (MULTIFURCATING) {
				  cluster_spr = rSPR_branch_and_bound_mult_range(cluster.F1,
						cluster.F2, min_spr, MAX_SPR - total_k);
				}
				else {
				  cluster_spr = rSPR_branch_and_bound_range(cluster.F1,
						cluster.F2, min_spr, MAX_SPR - total_k);
				}
				if (cluster_spr >= 0) {
					if (verbose) {
	  				cout << endl;
	  				cout << "F" << i << "_1: ";
	  				cluster.F1->print_components();
	  				cout << "F" << i << "_2: ";
	  				cluster.F2->print_components();
	  				cout << "cluster exact drSPR=" << cluster_spr << endl;
	  				cout << endl;
					}
					total_k += cluster_spr;
					k -= cluster_spr;
				}
				else {
					k = -1;
					if (verbose) {
						cout << "cluster exact drSPR=?  " << "k=" << k << " too large"
							<< endl;
						cout << "\n";
					}
					if (CLAMP) {
						total_k = MAX_SPR + 1;
					}
					else {
							total_k += min_spr;
					}
				}
			}
			if (k > -1) {
				if (!cluster.is_original()) {
					int adjustment = cluster.join_cluster(&F1, &F2);
					total_k += adjustment;
					delete cluster.F1;
					delete cluster.F2;
				}
			}
			else {
				if (!cluster.is_original()) {
					//if (cluster.F1_cluster_node != NULL)
					//	cluster.F1_cluster_node->contract();
					//if (cluster.F2_cluster_node != NULL)
					//	cluster.F2_cluster_node->contract();
					delete cluster.F1;
					delete cluster.F2;
				}
			}
		}
		delete cluster_points;
	}
	else {
		if (do_cluster) {
			delete cluster_points;
		}
		full_approx_spr /= 3;
		if (MULTIFURCATING) {
		  total_k = rSPR_branch_and_bound_mult_range(&F1, &F2, full_approx_spr, MAX_SPR);
		}
		else {
		  total_k = rSPR_branch_and_bound_range(&F1, &F2, full_approx_spr, MAX_SPR);
		}
		int i = 1;
		if (total_k < 0) {
			if (CLAMP) {
				total_k = MAX_SPR;
			}
			else {
		 	 total_k = full_approx_spr;
			}
		}

	}

	if (verbose) {
		F1.numbers_to_labels(reverse_label_map);
		F2.numbers_to_labels(reverse_label_map);
	 	cout << endl;
	 	cout << "F1: ";
	 	F1.print_components();
	 	cout << "F2: ";
	 	F2.print_components();
	 	cout << "total exact drSPR=" << total_k << endl;
	 	cout << endl;
	}
	return total_k;
}

// ---------------------------------------------------------------
// SEPARATOR
// ---------------------------------------------------------------


// KEY A
__attribute__((always_inline)) inline int rSPR_branch_and_bound_simple_clustering_Inline(
//SHADOWED ARGUMENTS
int &MAX_SPR,
bool &PREFER_RHO,
bool &MULTIFURCATING,
int &CLUSTER_TUNE,
bool &COUNT_LOSSES,
bool &SHOW_CLUSTERS,
int &MIN_SPR,
int &CLUSTER_MAX_SPR,
float &INITIAL_TREE_FRACTION,
bool &SPLIT_APPROX,
int &SPLIT_APPROX_THRESHOLD,
bool &CLAMP,

int (*rSPR_worse_3_mult_approx)(Forest *T1, Forest *T2),
int (*rSPR_worse_3_approx_2)(Forest *T1, Forest *T2),
int (*rSPR_worse_3_approx_3)(Node *subtree, Forest *T1, Forest *T2),
void (*reduction_leaf_mult)(Forest *T1, Forest* T2),
void (*reduction_leaf)(Forest *T1, Forest *T2),
int (*rSPR_branch_and_bound_mult)(Forest *T1, Forest *T2, int k),
int (*rSPR_branch_and_bound)(Forest *T1, Forest *T2, int k),
Node *(*find_subtree_of_approx_distance)(Node *n, Forest *F1, Forest *F2, int target_size),
int (*rSPR_branch_and_bound_hlpr)(Forest *T1, Forest *T2, int k,
		set<SiblingPair> *sibling_pairs, list<Node *> *singletons, bool cut_b_only,
		list<pair<Forest,Forest> > *AFs, list<Node *> *protected_stack,
		int *num_ties),



// ORIGINAL ARGUMENTS
Node *T1, Node *T2, bool verbose, map<string, int> *label_map, map<int, string> *reverse_label_map, int min_k, int max_k, Forest **out_F1, Forest **out_F2) {
	bool do_cluster = true;
	if (max_k > MAX_SPR)
		max_k = MAX_SPR;
	else if (max_k == -1)
		max_k = INT_MAX;

	if (T2->get_preorder_number() == -1) {
	  T2->preorder_number();
	}

	ClusterForest F1 = ClusterForest(T1);
	ClusterForest F2 = ClusterForest(T2);
	Forest F3 = Forest(F1);
	Forest F4 = Forest(F2);


//	bool old_rho = PREFER_RHO;
	PREFER_RHO = true;
	if (verbose) {
		cout << "T1: ";
		F1.print_components();
		cout << "T2: ";
		F2.print_components();
	}
	int full_approx_spr;
	if (MULTIFURCATING) {
	  full_approx_spr = rSPR_worse_3_mult_approx(&F3, &F4);
	}
	else {
	  full_approx_spr = rSPR_worse_3_approx_2(&F3, &F4);
	}
	if (full_approx_spr < CLUSTER_TUNE) {
		do_cluster = false;
	}
	if (verbose) {

		cout << "approx F1: ";
		F3.print_components();
		cout << "approx F2: ";
		F4.print_components();
		// what the AF shows
		cout << "approx drSPR=" << F4.num_components()-1 << endl;
		/* what we use to get the lower bound: 3 * the number of cutting rounds in
			 the approx algorithm
		*/
		//cout << "approx drSPR=" << full_approx_spr << endl;
		cout << "\n";
	}
//	if (F1.get_component(0)->get_preorder_number() == -1)
//		F1.get_component(0)->preorder_number();
//	if (F2.get_component(0)->get_preorder_number() == -1)
//		F2.get_component(0)->preorder_number();

	if (!sync_twins(&F1, &F2))
		return 0;
	if (F1.get_component(0)->is_leaf())
		return 0;
	if (F1.get_component(0)->get_preorder_number() == -1) {
		F1.get_component(0)->preorder_number();
		F2.get_component(0)->preorder_number();
	}
	int loss = 0;
	list<Node *> *cluster_points;
	if (F1.get_component(0)->get_edge_pre_start() == -1) {
		F1.get_component(0)->edge_preorder_interval();
		F2.get_component(0)->edge_preorder_interval();
	}
	if (LEAF_REDUCTION2) {

	  if (MULTIFURCATING) {
	    reduction_leaf_mult(&F1, &F2);
	  }
	  else {
	    reduction_leaf(&F1, &F2);
	    }
//		F1.get_component(0)->preorder_number();
//		F2.get_component(0)->preorder_number();
//		F1.get_component(0)->edge_preorder_interval();
//		F2.get_component(0)->edge_preorder_interval();
	}
	if (COUNT_LOSSES) {
		loss += F1.get_component(0)->count_lost_subtree();
		loss += F2.get_component(0)->count_lost_subtree();
	}
		//F1.print_components();
		//F2.print_components();
	if (do_cluster) {
		sync_interior_twins(&F1, &F2);
		cluster_points = find_cluster_points(&F1, &F2);
		//	list<Node *> *cluster_points = new list<Node *>();
		for(list<Node *>::iterator i = cluster_points->begin();
				i != cluster_points->end(); i++) {
			string cluster_name = "X";
			/*
			if (verbose) {
					stringstream ss;
					ss << F1.size();
					cluster_name += ss.str();
					//int num_labels = label_map.size();
					//label_map.insert(make_pair(cluster_name,num_labels));
					//reverse_label_map.insert(
					//		make_pair(num_labels,cluster_name));
					//ss.str("");
					//ss << num_labels;
					//cluster_name = ss.str();
			}
			*/

			Node *n = *i;
			if (n->parent()->parent() == NULL
					&& n->get_sibling() != NULL &&
					n->get_sibling()->get_name() == "X")
				continue;
			Node *n_parent = n->parent();
			Node *twin = n->get_twin();
			Node *twin_parent = twin->parent();

			if (twin_parent == NULL)
				continue;

			F1.add_cluster(n,cluster_name);

			F2.add_cluster(twin,cluster_name);

			Node *n_cluster =
					F1.get_cluster_node(F1.num_clusters()-1);
			Node *twin_cluster =
					F2.get_cluster_node(F2.num_clusters()-1);
			n_cluster->set_twin(twin_cluster);
			twin_cluster->set_twin(n_cluster);

		}
		if (verbose)
			cout << endl << "CLUSTERS" << endl;

	}

	// component 0 needs to be done last
	F1.add_component(F1.get_component(0));
	F2.add_component(F2.get_component(0));

	int k;
	int num_clusters = F1.num_components();
	int total_k = 0;

	if(SHOW_CLUSTERS){
		cout << "Clusters start" << endl;
		for(int i = 1; i < num_clusters; i++) {
			Forest f1 = Forest(F1.get_component(i));
			f1.print_components();
		}
		cout << "Clusters end" << endl;
	}


	for(int i = 1; i < num_clusters; i++) {
		if (i == num_clusters - 1) {
			PREFER_RHO = false;
		}
		int exact_spr = -1;
		//vector<Node *> comps = vector<Node *>();
		//comps.push_back(F1.get_component(i));
		Forest f1 = Forest(F1.get_component(i));
		//Forest f1 = Forest(comps);
		//comps.clear();

		//comps.push_back(F1.get_component(i));
		Forest f2 = Forest(F2.get_component(i));
		//Forest f2 = Forest(comps);
		//comps.clear();
		Forest f1a = Forest(f1);
		Forest f2a = Forest(f2);
		Forest *f1_cluster;
		Forest *f2_cluster;

		if (verbose) {
			cout << "C" << i << "_1: ";
			f1.print_components();
			cout << "C" << i << "_2: ";
			f2.print_components();
		}
		int approx_spr;
		if (MULTIFURCATING) {
		  approx_spr = rSPR_worse_3_mult_approx(&f1a, &f2a);
		}
		else {
		  approx_spr = rSPR_worse_3_approx_2(&f1a, &f2a);
		}
		if (verbose) {
			cout << "cluster approx drSPR=" << f2a.num_components()-1 << endl;
			//cout << "cluster approx drSPR=" << approx_spr << endl;

			cout << endl;
		}

		int min_spr = approx_spr / 3;
		if (min_spr < MIN_SPR - total_k)
			min_spr = MIN_SPR - total_k;
		int total_split_k = 0;

		bool done_cluster = false;
		bool done_split = false;

		double tree_fraction = INITIAL_TREE_FRACTION;

		if (min_spr < min_k)
			min_spr = min_k;

		while(!done_cluster) {
			done_cluster = true;

			for(k = min_spr - total_split_k; true; k++) {
				if (k < 0)
					k = 0;
				if (SPLIT_APPROX && !done_split && k >= SPLIT_APPROX_THRESHOLD) {
					done_cluster = false;
					break;
				}
				Forest f1t = Forest(f1);
//				Forest f1t = f1;
				Forest f2t = Forest(f2);
//				Forest f2t = f2;
				f1t.unsync();
				f2t.unsync();
				exact_spr = -1;
				if (verbose) {
					cout << k << " ";
  				cout.flush();
				}
				if (k + total_k <= max_k && k <= CLUSTER_MAX_SPR) {
					if (f1t.get_component(0)->get_name() == DEAD_COMPONENT) {
						f1t.add_rho();
						f2t.add_rho();
					}
					if (MULTIFURCATING) {
					  exact_spr = rSPR_branch_and_bound_mult(&f1t, &f2t, k);
					}
					else {
					  exact_spr = rSPR_branch_and_bound(&f1t, &f2t, k);
					}
				}
				if (exact_spr >= 0 || k + total_k > max_k ||
						k > CLUSTER_MAX_SPR) {
					if (k > CLUSTER_MAX_SPR) {
						f1t.swap(&f1a);
						f2t.swap(&f2a);
//						cout << "foo" << endl;
					}
					if (exact_spr >= 0) {
						exact_spr += total_split_k;
						if (verbose) {
	  					cout << endl;
	  					cout << "F" << i << "_1: ";
	  					f1t.print_components();
	  					cout << "F" << i << "_2: ";
	  					f2t.print_components();
	  					cout << "cluster exact drSPR=" << exact_spr << endl;
	  					cout << endl;
						}
						total_k += exact_spr;
					}
					else {
						// TODO: don't just the MAX_SPR here
						// incorporate extra information
						// toggle?
						if (verbose) {
							cout << "cluster exact drSPR=?  " << "k=" << k << " too large"
								<< endl;
							cout << "\n";
						}
						if (false && k > CLUSTER_MAX_SPR) {
							// TODO: this should be an approx of the remaining forest
//							total_k += approx_spr;
						}
						else if (CLAMP) {
							total_k = max_k;
						}
						else {
							Forest f1a = Forest(f1);
							Forest f2a = Forest(f2);

							int approx_spr;
							if (MULTIFURCATING) {
							  approx_spr = rSPR_worse_3_mult_approx(&f1a, &f2a);
							}
							else{
							  approx_spr = rSPR_worse_3_approx_2(&f1a, &f2a);
							}
								//total_k += min_spr;
								total_k += approx_spr / 3;
						}
					}
					if ( i < num_clusters - 1) {
						F1.join_cluster(i,&f1t);
						F2.join_cluster(i,&f2t);
					}
					else {
						F1.join_cluster(&f1t);
						F2.join_cluster(&f2t);
					}
					break;
				}
			}
			done_split = done_cluster;
			int num_splits = 0;
			while (SPLIT_APPROX && !done_split) {
				//IN_SPLIT_APPROX = true;
				Node *original_split_node = find_subtree_of_approx_distance(
						f1.get_component(0), &f1, &f2, SPLIT_APPROX_THRESHOLD*2);
				if (original_split_node == f1.get_component(0) &&
						num_splits > 0)
					done_split = true;
				else {
					Forest f1a = Forest(f1);
					Forest f2a = Forest(f2);
					Node *a_split_node =
					f1a.find_by_prenum(original_split_node->get_preorder_number());
					f1a.get_component(0)->disallow_siblings_subtree();
						a_split_node->allow_siblings_subtree();
//					if (a_split_node->lchild() != NULL)
//						a_split_node->lchild()->allow_siblings_subtree();
//					if (a_split_node->rchild() != NULL)
//						a_split_node->rchild()->allow_siblings_subtree();
					// something odd going on here
					int start = rSPR_worse_3_approx_3(a_split_node, &f1a, &f2a);
					if (start == INT_MAX)
						start = 0;
					start /= 3;
					int end = f1.get_component(0)->size();
					for(k = start; true; k++) {
						// TODO: figure out the bug here
						if (k > end) {
							k = 0;
							done_split = true;
							break;
						}
				/*	if (k > SPLIT_APPROX_THRESHOLD) {
						k = 0;
						tree_fraction *= 0.75;
						if (verbose)
							cout << "tree_fraction: " << tree_fraction << endl;
						continue;
					}*/
						Forest f1s = Forest(f1);
						Forest f2s = Forest(f2);
						if (!sync_twins(&f1s, &f2s)) {
							k = 0;
							done_split = true;
							break;
						}
						if (verbose) {
							cout << k << " ";
		  				cout.flush();
						}
						Node *split_node = f1s.find_by_prenum(original_split_node->get_preorder_number());
						f1s.get_component(0)->disallow_siblings_subtree();
							split_node->allow_siblings_subtree();
//						if (split_node->lchild() != NULL)
//							split_node->lchild()->allow_siblings_subtree();
//						if (split_node->rchild() != NULL)
//							split_node->rchild()->allow_siblings_subtree();
							//f1s.get_component(0)->find_subtree_of_size(tree_fraction);
							set<SiblingPair > *sibling_pairs =
								find_sibling_pairs_set(split_node);
							list<Node *> singletons = f2s.find_singletons();
							list<pair<Forest,Forest> > AFs = list<pair<Forest,Forest> >();
							list<Node *> protected_stack = list<Node *>();

							int num_ties = 2;

							int split_k = rSPR_branch_and_bound_hlpr(&f1s, &f2s, k,
									sibling_pairs, &singletons, false, &AFs,
									&protected_stack, &num_ties);
							delete sibling_pairs;
							if (!AFs.empty()) {
								AFs.front().first.swap(&f1);
								AFs.front().second.swap(&f2);
								f2.unprotect_edges();
								f1.get_component(0)->allow_siblings_subtree();
								AFs.clear();
								total_split_k += k - split_k;
		//						if (k < SPLIT_APPROX_THRESHOLD * 0.75) {
		//							tree_fraction *= 2;
		//							if (tree_fraction > INITIAL_TREE_FRACTION)
		//								tree_fraction = INITIAL_TREE_FRACTION;
		//						}
								if (verbose)
									cout << "split_k: " << k << endl;
								break;
							}
					}
				}
				//IN_SPLIT_APPROX = false;
				num_splits++;
			}

			// TODO: approx again? seperate approxes ?
		}
	}

		if (F1.contains_rho()) {
			F1.get_component(0)->delete_tree();
			F2.get_component(0)->delete_tree();
			F1.erase_components(0, num_clusters);
			F2.erase_components(0, num_clusters);
		}
		else {
			F1.get_component(num_clusters)->delete_tree();
			F2.get_component(num_clusters)->delete_tree();
			F1.erase_components(1, num_clusters+1);
			F2.erase_components(1, num_clusters+1);
		}
		// fix hanging roots
		for(int i = 0; i < F1.num_components(); i++) {
			F1.get_component(i)->contract(true);
			F2.get_component(i)->contract(true);
		}
		if (verbose) {
			F1.numbers_to_labels(reverse_label_map);
			F2.numbers_to_labels(reverse_label_map);
			cout << "F1: ";
			F1.print_components();
			cout << "F2: ";
			F2.print_components();
			cout << "total exact drSPR=" << total_k << endl;
		}
		if (out_F1 != NULL)
			*out_F1 = new Forest(&F1);
			//F1.swap(out_F1);
		if (out_F2 != NULL)
			*out_F2 = new Forest(&F2);
//			F2.swap(out_F2);
		if (out_F1 != NULL && out_F2 != NULL) {
//			out_F1->resync();
			sync_twins(*out_F1, *out_F2);
	//		sync_interior_twins_real(out_F1, out_F2);
		}

	if (do_cluster) {
		delete cluster_points;
	}
//	PREFER_RHO = old_rho;
	total_k += loss;
/*	cout << "F1: ";
	for (int i = 0; i < F1.num_components(); i++) {
		if (i > 0)
			cout << " ";
		F1.get_component(i)->expand_contracted_nodes();
		cout << F1.get_component(i)->str_edge_pre_interval_subtree();
	}
	cout << endl;
	cout << "F2: ";
	for (int i = 0; i < F2.num_components(); i++) {
		if (i > 0)
			cout << " ";
		F2.get_component(i)->expand_contracted_nodes();
		cout << F2.get_component(i)->str_edge_pre_interval_subtree();
	}
	cout << endl << endl;
*/
	#ifdef DEBUG_CASE_COUNTER
	print_mult_case_count();
        #endif

	return total_k;
}










}


