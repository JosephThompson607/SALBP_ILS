//
// Created by Joseph Thompson on 2025-06-12.
//
#include "albp_solution.h"
#include "ALBP.h"
#include "Hoff.h"
#include <iostream>
#include <cfloat>
#include <algorithm>
#include "salbp_basics.h"
Hoff::Hoff(const ALBP& albp, const float alpha, const float beta, const int max_attempts):
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



    ALBPSolution Hoff::solve() {
        ALBPSolution hoff_sol = ALBPSolution(albp_.N);
        hoff_sol.n_stations = 0;
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
                hoff_sol.task_assignment[i] = hoff_sol.n_stations;
                n_prec_[i] = -1;
                for (const int j : albp_.dir_suc[i]) {
                    n_prec_[j] --;
                }
            }
            total_assigned +=  best_s_task_assign_.size();
            hoff_sol.n_stations++;
            }

        hoff_sol.task_to_station();
        hoff_sol.station_to_ranking();
        hoff_sol.station_to_load(albp_);
        hoff_sol.find_windows(albp_);
        hoff_sol.cycle_time = *std::max_element(hoff_sol.loads.begin(), hoff_sol.loads.end());
        return hoff_sol;
    }


ALBPSolution hoff_solve(const ALBP &albp, int alpha_iter=4, int beta_iter =-1, float alpha_size = 0.005, float beta_size=0.005, bool reverse=true) {
    if (beta_iter < 0) {
        beta_iter = albp.N;
    }
    auto mhh= Hoff(albp);
    ALBPSolution best_result =mhh.solve();
    best_result.method = "mhh";
    for(int i = 0; i <=alpha_iter; i++) {
        float alpha = alpha_size * i;
        for (int j = 0; j < beta_iter; j++) {
            float beta = beta_size * j;
            auto mhh1= Hoff(albp, alpha, beta);

            if (ALBPSolution result =mhh1.solve(); result.n_stations < best_result.n_stations){
                std::cout << "best stations: " << best_result.n_stations << " new: "<<result.n_stations << std::endl;
                best_result = result;
            }
        }
    }
    //Reverse pass
    if (reverse) {
        ALBP rev_albp = albp.reverse();
        for(int i = 0; i <=alpha_iter; i++) {
            float alpha = alpha_size * i;
            for (int j = 0; j < beta_iter; j++) {
                float beta = beta_size* j;
                auto mhh1= Hoff(rev_albp, alpha, beta);

                if (ALBPSolution result =mhh1.solve(); result.n_stations < best_result.n_stations){
                    std::cout << "best stations: " << best_result.n_stations << " new: "<<result.n_stations << std::endl;
                    best_result = result;
                    best_result.reverse();
                    best_result.method = "reversed_hoff";
                }
            }
        }
    }


    return best_result;
}
ALBPSolution hoff_solve_salbp1(const ALBP &albp, int alpha_iter, int beta_iter, float alpha_size, float beta_size, bool reverse) {
    if (alpha_iter < 0) {
        alpha_iter = albp.N;
    }
    ALBPSolution best_result = hoff_solve(albp, alpha_iter, beta_iter, alpha_size, beta_size, reverse);
    return best_result;
}

ALBPSolution hoff_solve_salbp1(const int C,const int N, const std::vector<int>& task_times, const std::vector<std::vector<int>>& raw_precedence, int alpha_iter, int beta_iter, float alpha_size, float beta_size, bool reverse) {
    if (alpha_iter < 0) {
        alpha_iter = N;
    }
    ALBP albp = ALBP::type_1(C, N, task_times, raw_precedence);
    ALBPSolution best_result =hoff_solve(albp, alpha_iter, beta_iter, alpha_size, beta_size, reverse);
    return best_result;
}

void Hoff::gen_load( int depth, int remaining_time,const int start, int n_eligible, float cost) {

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

