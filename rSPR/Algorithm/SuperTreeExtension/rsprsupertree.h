//
// Created by ghillie-lich on 11/23/25.
//

#ifndef RSPR_PACE2026_RSPRSUPERTREE_H
#define RSPR_PACE2026_RSPRSUPERTREE_H

namespace rsprSupertree {

    __attribute__((always_inline)) inline int rSPR_total_distance(

    bool &MAIN_CALL,
    bool &PREFER_RHO,
    bool &VERBOSE,
    bool &FIND_RATE,
    bool &MULTIFURCATING,

    int (*rSPR_branch_and_bound_simple_clustering)(Node *T1, Node *T2, bool verbose),

    // ORIGINAL PARAMS
    Node *T1, vector<Node *> &gene_trees, vector<int> *original_scores) {
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



}



#endif //RSPR_PACE2026_RSPRSUPERTREE_H