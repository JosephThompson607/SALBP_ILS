//
// Created by Joseph Thompson on 2025-05-16.
//

#ifndef ILS_H
#define ILS_H
#include "albp_solution.h"
#include "ALBP.h"

ALBPSolution iterated_local_search(const ALBP& albp, int max_iter, float op_probs, bool verbose=false);

#endif //ILS_H
