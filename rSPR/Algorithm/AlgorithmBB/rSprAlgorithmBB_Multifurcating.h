//
// Created by Philip Kail on 11/22/25.
//







class bb_mult_recurse_data {
public:
	Forest *T1;
	Forest *T2;
	list<Node*> *sibling_groups;
	list<Node*> *singletons;
	map<Node *, Node*> node_map;
};



namespace rSprBB_Multifurcating {


/* Generates a copy of the data used to recurse on rSPR_branch_and_bound_mult_hlpr so
   that original forests aren't clobbered,
   updates pointers accordingly */
//
	__attribute__((always_inline)) inline bb_mult_recurse_data *generate_mult_recurse_data_Inline(
		Forest *T1, Forest *T2, list<Node*> *sibling_groups, list<Node*> *singletons) {

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

	__attribute__((always_inline)) inline void mult_cut_and_cleanup_Inline(
		Node* to_cut, Forest *T2, list<Node*> *singletons) {
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
	__attribute__((always_inline)) inline void mult_cut_all_except_and_cleanup_Inline(
		Node* T2_a1, Forest *T2, list<Node*> *singletons) {
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
	//__attribute__((always_inline)) inline
	void MULT_RHO_CUT_AND_RESOLVE(

	bb_mult_recurse_data *(*generate_mult_recurse_data)(Forest *T1, Forest *T2, list<Node*> *sibling_groups, list<Node*> *singletons),
	int (*rSPR_branch_and_bound_mult_hlpr)(Forest *T1, Forest *T2, int k, list<Node*> *sibling_groups, list<Node*> *singletons, Node *protected_node, list<pair<Forest,Forest>> *AFs, int* num_ties),

	
		Forest *T1, Forest *T2, int &k ,int &best_k, list<Node*> *sibling_groups,
		list<Node*> *singletons, list<pair<Forest,Forest>> *AFs,
						int* num_ties,
		Node *rho_to_singleton, Node *node_to_protect) {
		{
			bb_mult_recurse_data *next_data = generate_mult_recurse_data(T1, T2, sibling_groups, singletons);
			Node* protect = NULL;
			next_data->T1->add_rho();
			next_data->T2->add_rho();
			next_data->singletons->push_front(next_data->node_map[rho_to_singleton]);
			int result_k = rSPR_branch_and_bound_mult_hlpr(next_data->T1, next_data->T2, k - 1, next_data->sibling_groups, next_data->singletons, protect, AFs, num_ties);
			delete next_data->T1;
			delete next_data->T2;
			delete next_data->sibling_groups;
			delete next_data->singletons;
			delete next_data;
			if (result_k > best_k) {
				best_k = result_k;
			}
		}
	}


	//TODO (Ben): try using a routine instead of a macro, I suspect it will slow down because of
	//            stack and parameter passing though
	//__attribute__((always_inline)) inline
	void MULT_BB_CUT_AND_RESOLVE(
			bool &ALL_MAFS,

			int (*rSPR_branch_and_bound_mult_hlpr)(Forest *T1, Forest *T2, int k, list<Node*> *sibling_groups, list<Node*> *singletons, Node *protected_node, list<pair<Forest,Forest>> *AFs, int* num_ties),

			bb_mult_recurse_data *(*generate_mult_recurse_data)(Forest *T1, Forest *T2, list<Node*> *sibling_groups, list<Node*> *singletons),

			void (*mult_cut_and_cleanup)(Node* to_cut, Forest *T2, list<Node*> *singletons),

			void (*mult_cut_all_except_and_cleanup)(Node* T2_a1, Forest *T2, list<Node*> *singletons),


	Forest *T1, Forest *T2, int &k,int &best_k, list<Node*> *sibling_groups,
		list<Node*> *singletons, list<pair<Forest,Forest>> *AFs, int* num_ties,

		
	vector<Node*> &nodes_to_cut,  vector<Node*> &nodes_to_exclude_cutting, Node *node_to_protect) {
		bb_mult_recurse_data *next_data = generate_mult_recurse_data(T1, T2, sibling_groups, singletons);
		int num_cuts = 0;
		for (int i = 0; i < nodes_to_cut.size(); i++) {
			Node* T2_ax_next = next_data->node_map[nodes_to_cut[i]];
			mult_cut_and_cleanup(T2_ax_next, next_data->T2, next_data->singletons);
			num_cuts++;
		}
		for (int i = 0; i < nodes_to_exclude_cutting.size(); i++) {
			Node* T2_ax_next = next_data->node_map[nodes_to_exclude_cutting[i]];
			mult_cut_all_except_and_cleanup(T2_ax_next, next_data->T2, next_data->singletons);
			num_cuts++;
		}
		Node* protect = next_data->node_map[node_to_protect];
		int result_k = rSPR_branch_and_bound_mult_hlpr(
			next_data->T1, next_data->T2, k - num_cuts,
			next_data->sibling_groups, next_data->singletons,
			protect, AFs, num_ties);
		delete next_data->T1;
		delete next_data->T2;
		delete next_data->sibling_groups;
		delete next_data->singletons;
		delete next_data;
		if (result_k > best_k) {
			best_k = result_k;
			if (!ALL_MAFS && best_k > -1) {

			}
		}
	}


	
	//TODO: UndoMachine, then cleanup all constructors, bb_mult_recurse_data relying on copies of trees
    //__attribute__((always_inline)) inline
	int rSPR_branch_and_bound_mult_hlpr_Inline(
    // Shadowed Params
    bool &BB,
    bool &MULT_4_BRANCH,
    bool &ALL_MAFS,

    int (*rSPR_worse_3_mult_approx_hlpr)(Forest *T1, Forest *T2, list<Node *> *singletons, list<Node *> *sibling_pairs, Forest **F1, Forest **F2, bool save_forests),

    int (*rSPR_branch_and_bound_mult_hlpr)(Forest *T1, Forest *T2, int k, list<Node*> *sibling_groups, list<Node*> *singletons, Node *protected_node, list<pair<Forest,Forest>> *AFs, int* num_ties),
    bb_mult_recurse_data *(*generate_mult_recurse_data)(Forest *T1, Forest *T2, list<Node*> *sibling_groups, list<Node*> *singletons),
    void (*mult_cut_and_cleanup)(Node* to_cut, Forest *T2, list<Node*> *singletons),
    void (*mult_cut_all_except_and_cleanup)(Node* T2_a1, Forest *T2, list<Node*> *singletons),

	// Original Params
    Forest *T1, Forest *T2,
				    int &k,
				    list<Node*> *sibling_groups, list<Node*> *singletons,
				    Node *protected_node, list<pair<Forest,Forest>> *AFs,
				    int* num_ties){
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
	      MULT_BB_CUT_AND_RESOLVE(ALL_MAFS,rSPR_branch_and_bound_mult_hlpr,generate_mult_recurse_data,mult_cut_and_cleanup,mult_cut_all_except_and_cleanup, T1, T2, k, best_k, sibling_groups, singletons, AFs, num_ties, to_cut, to_cut_except, protected_node);
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
	      MULT_BB_CUT_AND_RESOLVE(ALL_MAFS,rSPR_branch_and_bound_mult_hlpr,generate_mult_recurse_data,mult_cut_and_cleanup,mult_cut_all_except_and_cleanup, T1, T2, k, best_k, sibling_groups, singletons, AFs, num_ties,to_cut, to_cut_except, protected_node);
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
		MULT_BB_CUT_AND_RESOLVE(ALL_MAFS,rSPR_branch_and_bound_mult_hlpr,generate_mult_recurse_data,mult_cut_and_cleanup,mult_cut_all_except_and_cleanup, T1, T2, k, best_k, sibling_groups, singletons, AFs, num_ties,to_cut, to_cut_except, protected_node);
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
	      MULT_BB_CUT_AND_RESOLVE(ALL_MAFS,rSPR_branch_and_bound_mult_hlpr,generate_mult_recurse_data,mult_cut_and_cleanup,mult_cut_all_except_and_cleanup, T1, T2, k, best_k, sibling_groups, singletons, AFs, num_ties,to_cut, to_cut_except, protected_node);	      
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
	      MULT_BB_CUT_AND_RESOLVE(ALL_MAFS,rSPR_branch_and_bound_mult_hlpr,generate_mult_recurse_data,mult_cut_and_cleanup,mult_cut_all_except_and_cleanup, T1, T2, k, best_k, sibling_groups, singletons, AFs, num_ties,to_cut, to_cut_except, protected_node);
	    }
	    /*

	      This is NOT in the cases. Figure out whats happening with 8.2 with r = 2

	     */
	    /*
	    {
	      cout << "Case 7.3 Cut T2_ax: " << T2_ax->str() << " Protected Node: " << protected_node->str() << endl;
	      vector<Node*> to_cut = {T2_ax};
	      vector<Node*> to_cut_except = {};
	      MULT_BB_CUT_AND_RESOLVE(ALL_MAFS,rSPR_branch_and_bound_mult_hlpr,generate_mult_recurse_data,mult_cut_and_cleanup,mult_cut_all_except_and_cleanup, T1, T2, k, best_k, sibling_groups, singletons, AFs, num_ties,to_cut, to_cut_except, protected_node);
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
	      MULT_BB_CUT_AND_RESOLVE(ALL_MAFS,rSPR_branch_and_bound_mult_hlpr,generate_mult_recurse_data,mult_cut_and_cleanup,mult_cut_all_except_and_cleanup, T1, T2, k, best_k, sibling_groups, singletons, AFs, num_ties,to_cut, to_cut_except, protected_node);
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
	      MULT_BB_CUT_AND_RESOLVE(ALL_MAFS,rSPR_branch_and_bound_mult_hlpr,generate_mult_recurse_data,mult_cut_and_cleanup,mult_cut_all_except_and_cleanup, T1, T2, k, best_k, sibling_groups, singletons, AFs, num_ties,to_cut, to_cut_except, protected_node);	      
	    }
	    /*

	      Not part of the outlined special cases

	     */
	    {
	      vector<Node*> to_cut = {T2_a2};
	      vector<Node*> to_cut_except = {};
	      MULT_BB_CUT_AND_RESOLVE(ALL_MAFS,rSPR_branch_and_bound_mult_hlpr,generate_mult_recurse_data,mult_cut_and_cleanup,mult_cut_all_except_and_cleanup, T1, T2, k, best_k, sibling_groups, singletons, AFs, num_ties,to_cut, to_cut_except, protected_node);	     
	    }
	    {
	      vector<Node*> to_cut = {};
	      vector<Node*> to_cut_except = {T2_a1};
	      MULT_BB_CUT_AND_RESOLVE(ALL_MAFS,rSPR_branch_and_bound_mult_hlpr,generate_mult_recurse_data,mult_cut_and_cleanup,mult_cut_all_except_and_cleanup, T1, T2, k, best_k, sibling_groups, singletons, AFs, num_ties,to_cut, to_cut_except, protected_node);	     
	    }

	    /*
	    if (T1_sibling_group->get_children().size() == 2 &&
		deepest_siblings.size() == 2 &&
		T2_a1 == T2_ax &&
		T2_a2 != protected_node) {
	      vector<Node*> to_cut = {T2_a2};
	      vector<Node*> to_cut_except = {};
	      MULT_BB_CUT_AND_RESOLVE(ALL_MAFS,rSPR_branch_and_bound_mult_hlpr,generate_mult_recurse_data,mult_cut_and_cleanup,mult_cut_all_except_and_cleanup, T1, T2, k, best_k, sibling_groups, singletons, AFs, num_ties,to_cut, to_cut_except, protected_node);
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
	      MULT_BB_CUT_AND_RESOLVE(ALL_MAFS,rSPR_branch_and_bound_mult_hlpr,generate_mult_recurse_data,mult_cut_and_cleanup,mult_cut_all_except_and_cleanup, T1, T2, k, best_k, sibling_groups, singletons, AFs, num_ties,to_cut, to_cut_except, NULL);
	    }
	    else if (T2_a1 == T2->get_component(0)) {
	      if (!T1->contains_rho()) {
		MULT_RHO_CUT_AND_RESOLVE(generate_mult_recurse_data,rSPR_branch_and_bound_mult_hlpr,T1, T2, k, best_k, sibling_groups, singletons, AFs, num_ties,T2_a1, NULL);
	      }
	    }
#ifdef DEBUG
	    cout << "Case 8.1b cut a2" << endl;
#endif
	    //8.1 : Cut a2
	    if (T2_a2->parent() != NULL) {
	      vector<Node*> to_cut = {T2_a2};
	      vector<Node*> to_cut_except = {};
	      MULT_BB_CUT_AND_RESOLVE(ALL_MAFS,rSPR_branch_and_bound_mult_hlpr,generate_mult_recurse_data,mult_cut_and_cleanup,mult_cut_all_except_and_cleanup, T1, T2, k, best_k, sibling_groups, singletons, AFs, num_ties,to_cut, to_cut_except, NULL);
	    }
	    else if (T2_a2 == T2->get_component(0)) {
	      if (!T1->contains_rho()) {
		MULT_RHO_CUT_AND_RESOLVE(generate_mult_recurse_data,rSPR_branch_and_bound_mult_hlpr,T1, T2, k, best_k, sibling_groups, singletons, AFs, num_ties,T2_a2, NULL);
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
	      MULT_BB_CUT_AND_RESOLVE(ALL_MAFS,rSPR_branch_and_bound_mult_hlpr,generate_mult_recurse_data,mult_cut_and_cleanup,mult_cut_all_except_and_cleanup, T1, T2, k, best_k, sibling_groups, singletons, AFs, num_ties,to_cut, to_cut_except, NULL);
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
		MULT_BB_CUT_AND_RESOLVE(ALL_MAFS,rSPR_branch_and_bound_mult_hlpr,generate_mult_recurse_data,mult_cut_and_cleanup,mult_cut_all_except_and_cleanup, T1, T2, k, best_k, sibling_groups, singletons, AFs, num_ties,to_cut, to_cut_except, deepest_siblings[i]);
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
		    MULT_BB_CUT_AND_RESOLVE(ALL_MAFS,rSPR_branch_and_bound_mult_hlpr,generate_mult_recurse_data,mult_cut_and_cleanup,mult_cut_all_except_and_cleanup, T1, T2, k, best_k, sibling_groups, singletons, AFs, num_ties,to_cut, to_cut_except, NULL);
		    break;
		  }
		}

	      }
	      }*/


	      else if (deepest_siblings.size() == 2) {
		  vector<Node*> to_cut = {T2_a2};
		  vector<Node*> to_cut_except = {};
		  MULT_BB_CUT_AND_RESOLVE(ALL_MAFS,rSPR_branch_and_bound_mult_hlpr,generate_mult_recurse_data,mult_cut_and_cleanup,mult_cut_all_except_and_cleanup, T1, T2, k, best_k, sibling_groups, singletons, AFs, num_ties,to_cut, to_cut_except, NULL);
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
	      MULT_BB_CUT_AND_RESOLVE(ALL_MAFS,rSPR_branch_and_bound_mult_hlpr,generate_mult_recurse_data,mult_cut_and_cleanup,mult_cut_all_except_and_cleanup, T1, T2, k, best_k, sibling_groups, singletons, AFs, num_ties,to_cut, to_cut_except, NULL);
	    }
	    //8.3 : Cut a2
	    if (T2_a2->parent() != NULL) {
#ifdef DEBUG
	    cout << "Case 8.3b cut a2" << endl;
#endif
	      vector<Node*> to_cut = {T2_a2};
	      vector<Node*> to_cut_except = {};
	      MULT_BB_CUT_AND_RESOLVE(ALL_MAFS,rSPR_branch_and_bound_mult_hlpr,generate_mult_recurse_data,mult_cut_and_cleanup,mult_cut_all_except_and_cleanup, T1, T2, k, best_k, sibling_groups, singletons, AFs, num_ties,to_cut, to_cut_except, NULL);
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
	      MULT_BB_CUT_AND_RESOLVE(ALL_MAFS,rSPR_branch_and_bound_mult_hlpr,generate_mult_recurse_data,mult_cut_and_cleanup,mult_cut_all_except_and_cleanup, T1, T2, k, best_k, sibling_groups, singletons, AFs, num_ties,to_cut, to_cut_except, NULL);
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
	      MULT_BB_CUT_AND_RESOLVE(ALL_MAFS,rSPR_branch_and_bound_mult_hlpr,generate_mult_recurse_data,mult_cut_and_cleanup,mult_cut_all_except_and_cleanup, T1, T2, k, best_k, sibling_groups, singletons, AFs, num_ties,to_cut, to_cut_except, NULL);
	    }
	    //8.4 : Cut a1 and a2's B's
	    {
#ifdef DEBUG
	    cout << "Case 8.4b cut B1 & B2" << endl;
#endif
	      vector<Node*> to_cut = {};
	      vector<Node*> to_cut_except = {T2_a1, T2_a2};
	      MULT_BB_CUT_AND_RESOLVE(ALL_MAFS,rSPR_branch_and_bound_mult_hlpr,generate_mult_recurse_data,mult_cut_and_cleanup,mult_cut_all_except_and_cleanup, T1, T2, k, best_k, sibling_groups, singletons, AFs, num_ties,to_cut, to_cut_except, NULL);
	    }

	    //8.4 : All except a1->p, then b1
	    {
#ifdef DEBUG
	    cout << "Case 8.4c cut B1 & B`1" << endl;
#endif
	      vector<Node*> to_cut = {};
	      vector<Node*> to_cut_except = {T2_a1->parent(), T2_a1};
	      MULT_BB_CUT_AND_RESOLVE(ALL_MAFS,rSPR_branch_and_bound_mult_hlpr,generate_mult_recurse_data,mult_cut_and_cleanup,mult_cut_all_except_and_cleanup, T1, T2, k, best_k, sibling_groups, singletons, AFs, num_ties,to_cut, to_cut_except, T2_a2);
	    }
	    //8.4 : All except a2->p, then b2
	    {
#ifdef DEBUG
	    cout << "Case 8.4d cut B2 & B`2" << endl;
#endif
	      vector<Node*> to_cut = {};
	      vector<Node*> to_cut_except = {T2_a2->parent(), T2_a2};
	      MULT_BB_CUT_AND_RESOLVE(ALL_MAFS,rSPR_branch_and_bound_mult_hlpr,generate_mult_recurse_data,mult_cut_and_cleanup,mult_cut_all_except_and_cleanup, T1, T2, k, best_k, sibling_groups, singletons, AFs, num_ties,to_cut, to_cut_except, T2_a1);
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
		MULT_BB_CUT_AND_RESOLVE(ALL_MAFS,rSPR_branch_and_bound_mult_hlpr,generate_mult_recurse_data,mult_cut_and_cleanup,mult_cut_all_except_and_cleanup, T1, T2, k, best_k, sibling_groups, singletons, AFs, num_ties,to_cut, to_cut_except, T2_a2);	
	      }
	      //8.4 Cut a2 prot a1
	      {
#ifdef DEBUG
	      cout << "Case 8.4f a2" << endl;
#endif
		vector<Node*> to_cut = {T2_a2};
		vector<Node*> to_cut_except = {};
		MULT_BB_CUT_AND_RESOLVE(ALL_MAFS,rSPR_branch_and_bound_mult_hlpr,generate_mult_recurse_data,mult_cut_and_cleanup,mult_cut_all_except_and_cleanup, T1, T2, k, best_k, sibling_groups, singletons, AFs, num_ties,to_cut, to_cut_except, T2_a1);	
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
	      MULT_BB_CUT_AND_RESOLVE(ALL_MAFS,rSPR_branch_and_bound_mult_hlpr,generate_mult_recurse_data,mult_cut_and_cleanup,mult_cut_all_except_and_cleanup, T1, T2, k, best_k, sibling_groups, singletons, AFs, num_ties,to_cut, to_cut_except, NULL);
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
	      MULT_BB_CUT_AND_RESOLVE(ALL_MAFS,rSPR_branch_and_bound_mult_hlpr,generate_mult_recurse_data,mult_cut_and_cleanup,mult_cut_all_except_and_cleanup, T1, T2, k, best_k, sibling_groups, singletons, AFs, num_ties,to_cut, to_cut_except, NULL);
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
		MULT_BB_CUT_AND_RESOLVE(ALL_MAFS,rSPR_branch_and_bound_mult_hlpr,generate_mult_recurse_data,mult_cut_and_cleanup,mult_cut_all_except_and_cleanup, T1, T2, k, best_k, sibling_groups, singletons, AFs, num_ties,to_cut, to_cut_except, deepest_siblings[i]);//TODO protect ai
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
		MULT_BB_CUT_AND_RESOLVE(ALL_MAFS,rSPR_branch_and_bound_mult_hlpr,generate_mult_recurse_data,mult_cut_and_cleanup,mult_cut_all_except_and_cleanup, T1, T2, k, best_k, sibling_groups, singletons, AFs, num_ties,to_cut, to_cut_except, deepest_siblings[i]);//TODO protect ai
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
	      MULT_BB_CUT_AND_RESOLVE(ALL_MAFS,rSPR_branch_and_bound_mult_hlpr,generate_mult_recurse_data,mult_cut_and_cleanup,mult_cut_all_except_and_cleanup, T1, T2, k, best_k, sibling_groups, singletons, AFs, num_ties,to_cut, to_cut_except, NULL);
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

	      MULT_BB_CUT_AND_RESOLVE(ALL_MAFS,rSPR_branch_and_bound_mult_hlpr,generate_mult_recurse_data,mult_cut_and_cleanup,mult_cut_all_except_and_cleanup, T1, T2, k, best_k, sibling_groups, singletons, AFs, num_ties,to_cut, to_cut_except, NULL);
	    }
	    //8.6 cut B2
	    {
#ifdef DEBUG	     
	      cout << "Case 8.6c cut B2" << endl;
#endif
	      vector<Node*> to_cut = {};
	      vector<Node*> to_cut_except = {T2_a2};
	      MULT_BB_CUT_AND_RESOLVE(ALL_MAFS,rSPR_branch_and_bound_mult_hlpr,generate_mult_recurse_data,mult_cut_and_cleanup,mult_cut_all_except_and_cleanup, T1, T2, k, best_k, sibling_groups, singletons, AFs, num_ties,to_cut, to_cut_except, NULL);
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
	      MULT_BB_CUT_AND_RESOLVE(ALL_MAFS,rSPR_branch_and_bound_mult_hlpr,generate_mult_recurse_data,mult_cut_and_cleanup,mult_cut_all_except_and_cleanup, T1, T2, k, best_k, sibling_groups, singletons, AFs, num_ties,to_cut, to_cut_except, NULL);
	    }	      
	    //8.7 cut a2
	    {
#ifdef DEBUG	     
	      cout << "Case 8.7b cut a2" << endl;
#endif
	      vector<Node*> to_cut = {T2_a2};
	      vector<Node*> to_cut_except = {};
	      MULT_BB_CUT_AND_RESOLVE(ALL_MAFS,rSPR_branch_and_bound_mult_hlpr,generate_mult_recurse_data,mult_cut_and_cleanup,mult_cut_all_except_and_cleanup, T1, T2, k, best_k, sibling_groups, singletons, AFs, num_ties,to_cut, to_cut_except, T2_a1);
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

	      MULT_BB_CUT_AND_RESOLVE(ALL_MAFS,rSPR_branch_and_bound_mult_hlpr,generate_mult_recurse_data,mult_cut_and_cleanup,mult_cut_all_except_and_cleanup, T1, T2, k, best_k, sibling_groups, singletons, AFs, num_ties,to_cut, to_cut_except, NULL);		
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

	      MULT_BB_CUT_AND_RESOLVE(ALL_MAFS,rSPR_branch_and_bound_mult_hlpr,generate_mult_recurse_data,mult_cut_and_cleanup,mult_cut_all_except_and_cleanup, T1, T2, k, best_k, sibling_groups, singletons, AFs, num_ties,to_cut, to_cut_except, T2_a1);	
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
	      MULT_BB_CUT_AND_RESOLVE(ALL_MAFS,rSPR_branch_and_bound_mult_hlpr,generate_mult_recurse_data,mult_cut_and_cleanup,mult_cut_all_except_and_cleanup, T1, T2, k, best_k, sibling_groups, singletons, AFs, num_ties,to_cut, to_cut_except, NULL);
	    }
	    //mimics case 8.1, should add rho
	    else if (T2_a1 == T2->get_component(0) && arbitrary_lca == NULL) {
	      if (!T1->contains_rho()) {
		MULT_RHO_CUT_AND_RESOLVE(generate_mult_recurse_data,rSPR_branch_and_bound_mult_hlpr,T1, T2, k, best_k, sibling_groups, singletons, AFs, num_ties,T2_a1, NULL);
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
	      MULT_BB_CUT_AND_RESOLVE(ALL_MAFS,rSPR_branch_and_bound_mult_hlpr,generate_mult_recurse_data,mult_cut_and_cleanup,mult_cut_all_except_and_cleanup, T1, T2, k, best_k, sibling_groups, singletons, AFs, num_ties,to_cut, to_cut_except, NULL);
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
	      MULT_BB_CUT_AND_RESOLVE(ALL_MAFS,rSPR_branch_and_bound_mult_hlpr,generate_mult_recurse_data,mult_cut_and_cleanup,mult_cut_all_except_and_cleanup, T1, T2, k, best_k, sibling_groups, singletons, AFs, num_ties,to_cut, to_cut_except, NULL);
	    }
	    else if (T2_a2 == T2->get_component(0) && arbitrary_lca == NULL) {
	      
	      if (!T1->contains_rho()) {
		MULT_RHO_CUT_AND_RESOLVE(generate_mult_recurse_data,rSPR_branch_and_bound_mult_hlpr,T1, T2, k, best_k, sibling_groups, singletons, AFs, num_ties,T2_a2, NULL);
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
	      MULT_BB_CUT_AND_RESOLVE(ALL_MAFS,rSPR_branch_and_bound_mult_hlpr,generate_mult_recurse_data,mult_cut_and_cleanup,mult_cut_all_except_and_cleanup, T1, T2, k, best_k, sibling_groups, singletons, AFs, num_ties,to_cut, to_cut_except, NULL);
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



	

	__attribute__((always_inline)) inline int rSPR_branch_and_bound_mult_range_Inline(
	bool &LEAF_REDUCTION,

	void (*reduction_leaf_mult)(Forest *T1, Forest* T2),
	int (*rSPR_branch_and_bound_mult)(Forest *T1, Forest *T2, int start_k),

	Forest *T1, Forest *T2, int start_k, int end_k){
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


	__attribute__((always_inline)) inline int rSPR_branch_and_bound_mult_Inline(
	bool &ALL_MAFS,
	int (*rSPR_branch_and_bound_mult_hlpr)(Forest *T1, Forest *T2, int k, list<Node*> *sibling_groups, list<Node*> *singletons, Node *protected_node, list<pair<Forest,Forest>> *AFs, int* num_ties),


	Forest *T1, Forest *T2, int k){

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


}



