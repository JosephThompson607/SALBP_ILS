//
// Created by Joseph Thompson on 2025-06-12.
//
#include "albp_solution.h"
#include "ALBP.h"
#include "mhh.h"

#include <cfloat>

MHH::MHH(const ALBP& albp, const float alpha, const float beta, const int max_attempts):
    albp_(albp),


     n_prec_(albp.N, 0),                    // Size = num_tasks, init to 0
     alpha_(alpha),                                        // Initialize to 0
     beta_(beta),
    max_attempts_(max_attempts),
    n_attempts_(0),
    // Initialize to 0
     min_cost_(FLT_MAX)
    {
    s_task_assign_.reserve(albp.N),
    best_s_task_assign_.reserve(albp.N),
    eligible_tasks_.resize(albp.N,-1);
    pw_ = get_positional_weight(albp_) ;
    //initilize counts of unassigned predecessors
    for(int i = 0; i < albp_.N; i++) {
        n_prec_[i] = albp_.dir_pred[i].size();
    }

}



    ALBPSolution MHH::solve() {
        ALBPSolution mhh_sol = ALBPSolution(albp_.N);
        mhh_sol.n_stations = 0;





        int total_assigned = 0;
        while (total_assigned < albp_.N) {
            s_task_assign_.clear();
            best_s_task_assign_.clear();
            int n_eligible = 0;
            n_attempts_ = 0;
            for(int i = 0; i < albp_.N; i++) {
                if (n_prec_[i] == 0) {
                    eligible_tasks_[n_eligible++] = i;
                }
            }
            min_cost_ = FLT_MAX;
            gen_load(0, albp_.C, 0, n_eligible, albp_.C );

            for (const int i : best_s_task_assign_) {
                mhh_sol.task_assignment[i] = mhh_sol.n_stations;
                n_prec_[i] = -1;
                for (const int j : albp_.dir_suc[i]) {
                    n_prec_[j] --;
                }
            }
            total_assigned +=  best_s_task_assign_.size();
            mhh_sol.n_stations++;
            }
        mhh_sol.task_to_station();
        mhh_sol.station_to_ranking();
        return mhh_sol;
    }


    ALBPSolution mhh_solve_salbp1(int C, int N, const std::vector<int>& task_times, const std::vector<std::vector<int>>& raw_precedence) {
        ALBP albp(C, N, task_times, raw_precedence);
        for(float i = 0; i < albp.N; i++) {}
        auto mhh= MHH(albp);
        ALBPSolution result =mhh.solve();
        return result;
    }

    void MHH::gen_load( int depth, int remaining_time,const int start, int n_eligible, float cost) {

        int full_load = 1;
        for(int i=start;i<n_eligible;i++) {
            if (int task = eligible_tasks_[i]; albp_.task_time[task] <= remaining_time) {
                full_load = 0;
                int n_sub_eligible = n_eligible;
                s_task_assign_.push_back(task);
                n_prec_[task] = -1;
                for (const int j : albp_.dir_suc[task]) {
                    n_prec_[j] --;
                    if (n_prec_[j] == 0) {
                        eligible_tasks_[n_sub_eligible++] = j;
                    }
                }
                int sub_remaining_time = remaining_time - albp_.task_time[task];
                float sub_cost = cost - albp_.task_time[task] - alpha_ * pw_[task] - beta_ * albp_.suc[task].size();
                if (sub_cost < min_cost_) {
                    min_cost_ = sub_cost;
                    best_s_task_assign_ = s_task_assign_;
                }
                gen_load(depth+1, sub_remaining_time, i+1, n_sub_eligible, sub_cost );

                //undo the changes
                s_task_assign_.pop_back();
                n_prec_[task] = 0;
                for (const int j : albp_.dir_suc[task]) {
                    n_prec_[j] ++;
                }
                if (n_attempts_ >= max_attempts_) return;
            }



        }


    n_attempts_+=full_load;

};

