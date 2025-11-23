//
// Created by user on 11/23/25.
//

#include "Forest.h"

#ifndef RSPR_PACE2026_PROBLEMSOLUTION_H
#define RSPR_PACE2026_PROBLEMSOLUTION_H

#endif //RSPR_PACE2026_PROBLEMSOLUTION_H


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