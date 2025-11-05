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
Hoff::Hoff(const ALBP& albp,  int alpha_iter, int beta_iter , float alpha_size, float beta_size, bool reverse, int max_attempts):
    albp_(albp),
    alpha_iter_(alpha_iter),
    beta_iter_(beta_iter),
    alpha_size_(alpha_size),
    beta_size_(beta_size),
    reverse_(reverse),
     n_prec_(albp.N, 0),                    // Size = num_tasks, init to 0
     alpha_(0),                                        // Initialize to 0
     beta_(0),
    max_attempts_(max_attempts),
    n_attempts_(0),
    // Initialize to 0
     min_cost_(FLT_MAX)
    {
    lb_ = calc_salbp_1_bin_lbs(albp.task_time, albp.C);
    s_task_assign_.reserve(albp.N),
    best_s_task_assign_.reserve(albp.N),
    pw_ = get_positional_weight(albp_) ;
    //initilize counts of unassigned predecessors
    reinitialize();

}

    void Hoff::reinitialize() {
    eligible_tasks_.resize(albp_.N,-1);
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

ALBPSolution Hoff::multi_solve() {
    std::cout <<"starting multi solve" << std::endl;
    ALBPSolution best_result =solve();
    best_result.method = "hoff";
    std::cout << "changing alpha and beta current best " << best_result.n_stations << std::endl;
    for(int i = 0; i <=alpha_iter_; i++) {
        alpha_ = alpha_size_ * i;
        for (int j = 0; j < beta_iter_; j++) {
            beta_ = beta_size_ * j;
            //Check for optimality
            if (best_result.n_stations == lb_) {
                std::cout << "optimal solution !" << std::endl;
                best_result.optimal = true;
                return best_result;
            }
            //Resets original data and re-solves
            reinitialize();
            if (ALBPSolution result =solve(); result.n_stations < best_result.n_stations){
                //std::cout << "best stations: " << best_result.n_stations << " new: "<<result.n_stations << std::endl;
                best_result = result;
                if (best_result.n_stations == lb_) {

                }
            }
        }
    }
    std::cout <<"ending multi solve best result has " << best_result.n_stations <<std::endl;

    return best_result;
}

ALBPSolution hoff_solve(const ALBP &albp, int alpha_iter=4, int beta_iter =-1, float alpha_size = 0.005, float beta_size=0.005, bool reverse=true){
    if (beta_iter < 0) {
        beta_iter = albp.N;
    }
    auto hoff = Hoff(albp, alpha_iter, beta_iter, alpha_size, beta_size);
    ALBPSolution current_best =  hoff.multi_solve();

    if (reverse && current_best.optimal == false){
        ALBP rev_albp = albp.reverse();
        auto hoff = Hoff(rev_albp, alpha_iter, beta_iter, alpha_size, beta_size);
        ALBPSolution rev_sol = hoff.solve();
        if (rev_sol.n_stations < current_best.n_stations){

        current_best = rev_sol;}
    }


    return current_best;
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
            if ((n_attempts_ >= max_attempts_) || (remaining_time==0)) return;
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

            }
        }
    n_attempts_+=full_load;
};

