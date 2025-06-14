//
// Created by Joseph Thompson on 2025-06-12.
//

#ifndef MHH_H
#define MHH_H


#include "albp_solution.h"
#include "ALBP.h"

class MHH {
    public:
        explicit MHH(const ALBP& albp, float alpha=0.0f, float beta=0.0f, int max_attempts = 5000);
        ALBPSolution solve( );

private:
    const ALBP& albp_;
    void gen_load(int depth, int remaining_time, int start, int n_eligible, float cost);
    std::vector<int> eligible_tasks_{}; //tasks currently available for assignment
    std::vector<int> s_task_assign_{}; //task assignments to a given station
    std::vector<int> best_s_task_assign_{}; //task assignments to a given station
    std::vector<int> n_prec_{}; //number of predecessor unassigned. 0 if available, 1 if not
    int n_attempts_;
    int max_attempts_;
    float alpha_, beta_, min_cost_=0;

    std::vector<int>  pw_{};
};
ALBPSolution mhh_solve_salbp1(int C, int N, const std::vector<int>& task_times, const std::vector<std::vector<int>>& raw_precedence);
#endif //MHH_H
