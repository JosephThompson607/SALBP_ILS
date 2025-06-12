//
// Created by Joseph Thompson on 2025-06-12.
//

#ifndef MHH_H
#define MHH_H
#include "albp_solution.h"
#include "ALBP.h"


ALBPSolution modified_hoffman(const ALBP& albp , int max_loads = 5000);
ALBPSolution mhh_solve_salbp1(int C, int N, const std::vector<int>& task_times, const std::vector<std::vector<int>>& raw_precedence, bool verbose=false);



#endif //MHH_H
