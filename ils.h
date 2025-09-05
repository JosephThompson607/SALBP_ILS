//
// Created by Joseph Thompson on 2025-05-16.
//

#ifndef ILS_H
#define ILS_H
#include "albp_solution.h"
#include "ALBP.h"

ALBPSolution iterated_local_search(const ALBP& albp, int max_iter, int time_limit, float op_probs, bool verbose=false,const std::vector<int> &initial_solution = std::vector<int>() );
ALBPSolution ils_solve_SALBP1(int C, int N, const std::vector<int>& task_times, const std::vector<std::vector<int>>& raw_precedence, int max_iter, std::optional<int> time_limit, float op_probs, bool verbose=false,  const std::vector<int> &initial_solution = std::vector<int>());

#endif //ILS_H
