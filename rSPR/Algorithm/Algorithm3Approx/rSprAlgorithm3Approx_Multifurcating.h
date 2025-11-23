//
// Created by ghillie-lich on 11/23/25.
//

namespace rSpr3Approx_Multifurcating {
    __attribute__((always_inline)) inline void reduction_leaf_mult_Inline(Forest *T1, Forest* T2) {

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
}

