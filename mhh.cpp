//
// Created by Joseph Thompson on 2025-06-12.
//
#include "albp_solution.h"
#include "ALBP.h"
#include "mhh.h"

void gen_load(const ALBP&,std::vector<int> eligible_tasks,std::vector<int> task_assignments, int n_eligible, float cost) {

}
ALBPSolution modified_hoffman(const ALBP& albp, int max_loads ) {
    ALBPSolution mhh_sol = ALBPSolution(albp.N);
    int n_assigned = 0;

    int n_stations = 0;
    std::vector<int> n_prec(albp.N, 0);
    std::vector<int> eligible_tasks(albp.N, -1);
    std::vector<int> task_assignments(albp.N, -1);
    while (n_assigned < albp.N) {
        int n_eligible = 0;
        for(int i = 0; i < albp.N; i++) {
            n_prec[i] = albp.dir_pred[i].size();
            if (n_prec[i] == 0) {
                eligible_tasks[n_eligible++] = i;
            }
        }
        gen_load(albp, eligible_tasks, task_assignments, n_eligible, cost );


    }


    return mhh_sol;
}
ALBPSolution mhh_solve_salbp1(int C, int N, const std::vector<int>& task_times, const std::vector<std::vector<int>>& raw_precedence, bool verbose) {
    ALBP albp(C, N, task_times, raw_precedence);
    ALBPSolution result =modified_hoffman(albp);
    return result;
}

