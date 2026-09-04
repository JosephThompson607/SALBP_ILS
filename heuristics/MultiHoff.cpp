//
// Created by Joseph Thompson on 2025-06-12.
//
#include "../albp_solution.h"
#include "../ALBP.h"
#include "MultiHoff.h"
#include <iostream>
#include <cfloat>
#include <algorithm>
#include "salbp_basics.h"
#include <map>
#include <optional>
#include <numeric>
MultiHoff::MultiHoff(const ALBP& albp, const int max_attempts,
                    const std::optional<std::vector<float>>& alpha_schedule,
                const std::optional<std::vector<float>>& beta_schedule,
                const std::optional<float> gamma,
                const std::optional<std::vector<int>>& rankings,
                const std::optional<unsigned> seed):

    albp_(albp),
    mhh_sol_(albp.N),
    remaining_task_times_(albp.task_time),
    dir_pred_(albp.dir_pred),
    dir_suc_(albp.dir_suc),
    n_prec_(albp.N, 0),
    n_suc_(albp.N, 0),
    n_prec_orig_(albp.N, 0),
    n_suc_orig_(albp.N, 0),
    rng_(seed ? *seed : std::random_device{}()),
    n_attempts_(0),
    max_attempts_(max_attempts),
    min_cost_(FLT_MAX),


    alpha_sched_(alpha_schedule.value_or(std::vector<float>{0.0f})),
    beta_sched_(beta_schedule.value_or(std::vector<float>{0.0f})),
    gamma_(gamma.value_or(0.0)),
    //Making gamma_ perturbation scale with smallest task time
    pert_size_(static_cast<float>(*std::min_element(albp.task_time.begin(), albp.task_time.end()))),

    ub_(albp.N),
     lb_(calc_salbp_1_bin_lbs(albp.task_time, albp.C))
    {


    s_task_assign_.reserve(albp.N);

    if (!rankings){
    forw_ranking_ = pw_ranking(albp_) ; //Ranking of tasks in forward/backward direction for hoffman
    back_ranking_ = rpw_ranking(albp_) ;
        }
    else {
        forw_ranking_ = *rankings;
        back_ranking_ = *rankings;
    }
    initialize_current_s_assignments();

    //initialize counts of unassigned predecessors

    
    for(int i = 0; i < albp_.N; i++) {
        n_prec_orig_[i] = dir_pred_[i].size();
        if (n_prec_orig_[i] == 0) {
            no_prec_tasks_.push_back(i);
        }
    }
    for (int i = 0; i < albp.N; i++) {
        n_suc_orig_[i] = dir_suc_[i].size();
        if (n_suc_orig_[i] == 0) {
            no_suc_tasks_.push_back(i);
        }

    }

}
    void MultiHoff::initialize_current_s_assignments() {
    // Initialize all needed keys
    reverse_s_assignment_.clear();
    for (int i = 0; i < albp_.N; ++i) {
        reverse_s_assignment_[i] = { -1};
    }
}


    std::vector<int> MultiHoff::filter_eligible(std::vector<int>& elig) {
        std::vector<int> result;
        result.reserve(elig.size());
        for (int i =0; i < elig.size(); i++) {
            int task = elig[i];
            if (n_prec_[task] != -1) {
                result.push_back(task);
            }
        }
        return result;
}



bool MultiHoff::check_ub() const {
    //Get unassigned tasks

    int lb = calc_salbp_1_bin_lbs(remaining_task_times_, albp_.C);
    //std::cout << "forward station " << forward_station_ << " backward station " << backward_station_<<"The lb is" << lb << std::endl;
    if (ub_ < forward_station_ + backward_station_ + lb +1 ) {
        //std::cout <<"skipping" << std::endl;
        return true;
    }
    return false;

}

void MultiHoff::mark_task_assigned(const int task, std::vector<int>& elig, const bool remove_old) {
    if (remove_old) {
        auto it = std::find(elig.begin(), elig.end(), task);
        if (it != elig.end()) {
            *it = elig.back();
            elig.pop_back();
        }
    }
    n_prec_[task] = -1;
    n_suc_[task] = -1;
    remaining_task_times_[task] = 0;
}


int MultiHoff::one_packing_search( std::vector<int>&elig, const int station) {
    s_task_assign_.clear();
    best_s_task_assign_.clear();
    n_attempts_ = 0;
    min_cost_ = FLT_MAX;
    gen_load(0, albp_.C, 0, albp_.C, elig );
    if (back_pass_) {
        s_backwards_.push_front(best_s_task_assign_);
    }
    else {
        s_forwards_.push_back(best_s_task_assign_);
    }
    int last_task = std::numeric_limits<int>::max();
    for (const int task : best_s_task_assign_) {
        //Remove task from consideration
        mark_task_assigned(task, elig);
        //Adds task to eligible list
        add_new_available(elig, task);
        //Updates records for caching
        if (back_pass_){
            reverse_s_assignment_[task] = {station};
        }
        else {
            if (reverse_s_assignment_[task] > -1) {
                last_task = std::min(last_task, reverse_s_assignment_[task]);
            }
        }

    }
    if (check_ub()) {
        return -1;
    }
    return last_task;
}


    ALBPSolution MultiHoff::solve_one_pass() {
    bool improved=false;

    for (int mf =0; mf <ub_; mf++) {
        back_pass_=false;
        std::vector<int> eligible_tasks_forward = no_prec_tasks_;
        std::vector<int> eligible_tasks_backward = no_suc_tasks_;
        forward_station_ = 0;
        backward_station_ = 0;
        n_prec_ = n_prec_orig_;
        n_suc_ = n_suc_orig_;
        remaining_task_times_ = albp_.task_time;
        int n_forward_assigned_tasks = 0;
        if (mf > 1 ) { //reloading cached data for forward pass
            for (int s = 0; s < mf-1; s++) {
                for (const int task : s_forwards_[s]){
                    mark_task_assigned(task, eligible_tasks_forward);
                    add_new_available(eligible_tasks_forward, task);
                    n_forward_assigned_tasks++;
                }
            }
            forward_station_ = mf -1;
        }
        //Hoffman pass forward
        int last_task = std::numeric_limits<int>::max(); //-1 if ub violation, or last cached station to remove
        while (forward_station_ < mf && !eligible_tasks_forward.empty() && last_task != -1) {
            sort_by_ranking(eligible_tasks_forward, forw_ranking_);
            last_task =one_packing_search(eligible_tasks_forward, forward_station_);
            ++forward_station_;
        }
        //BACKWARDS PASS
        back_pass_ = true;
        int n_cached = 0;
        if (mf > 0 && last_task != -1 ) {
            //Remove tasks from eligible forward that were already assigned
            eligible_tasks_backward = filter_eligible(eligible_tasks_backward);
            //using catched partial solutions from previous passes
            //Remove all stations from cache that were impacted by forward pass
            if (s_backwards_.size() > last_task) {
                s_backwards_.erase(s_backwards_.begin(), s_backwards_.end() - last_task);
            }
            //Add remaining station assignments from previous solutions
            for (int i=s_backwards_.size()-1; i >= 0; i--) {
                for (int task: s_backwards_[i]) {
                    n_cached++;
                    mark_task_assigned(task, eligible_tasks_backward);
                    add_new_available(eligible_tasks_backward, task);

                }
            }

            backward_station_ = s_backwards_.size();
        }
        //hoffman pass packwards
        while (!eligible_tasks_backward.empty()&&last_task != -1) {
            sort_by_ranking(eligible_tasks_backward, back_ranking_);
            last_task = one_packing_search(eligible_tasks_backward, backward_station_);
            ++backward_station_;
        }
        int n_stations = forward_station_ + backward_station_;
        // std::cout << "combined has solution has "<< n_stations << "last task "<< last_task <<std::endl;
        // std::cout << "Solution printout forward" << std::endl;
        // for (int s = 0; s < mf; s++) {
        //     for (const int task : s_forwards_[s]) {
        //         std::cout << "station task" << s << " " << task << std::endl;
        //     }
        // }
        // std::cout << "Solution printout backward" << std::endl;
        // for (int s = 0; s<backward_station_; s++) {
        //     for (const int task : s_backwards_[s]) {
        //         std::cout << "station task" << s + mf << " " << task << std::endl;
        //     }
        // }


        if (n_stations < ub_ && last_task != -1) {
            improved = true;
            ub_ = n_stations;
            for (int s = 0; s < mf; s++) {
                for (const int task : s_forwards_[s]) {
                    mhh_sol_.task_assignment[task] = s;
                }
            }
            for (int s = 0; s<backward_station_; s++) {
                for (const int task : s_backwards_[s]) {
                    mhh_sol_.task_assignment[task] = mf + s;
                }
            }
        }

        if (ub_==lb_) {
            mhh_sol_.optimal =true;
            break;
        }
    }
    if (improved){ //Save new solution if we have an improvement
        if (reverse_) {
            std::reverse(mhh_sol_.task_assignment.begin(), mhh_sol_.task_assignment.end());
        }
        mhh_sol_.n_stations = ub_;
        mhh_sol_.task_to_station_and_load(albp_);
        mhh_sol_.station_to_ranking(false);
        return mhh_sol_;
    }
        return mhh_sol_;

}


void MultiHoff::reverse_solve_order() {
    /* This function reverses the elements necessary for a backwards solve. Note that it modifies some internal state*/
    reverse_ = !reverse_;
    back_pass_ = false;
    std::swap(dir_pred_, dir_suc_);
    std::swap(n_prec_orig_, n_suc_orig_);
    std::swap(forw_ranking_, back_ranking_);
    std::swap(no_prec_tasks_, no_suc_tasks_);
    initialize_current_s_assignments();
    s_forwards_.clear();
    s_backwards_.clear();
    remaining_task_times_ = albp_.task_time;



}

ALBPSolution MultiHoff::solve() {
    float first_alpha = alpha_sched_.front();
    float first_beta = beta_sched_.front();

    alpha_ = first_alpha;
    beta_ = first_beta;
    ALBPSolution best_result = solve_one_pass();

    if (ub_ != lb_) {

        reverse_solve_order();
        ALBPSolution best_backward = solve_one_pass();
        if (best_backward.n_stations < best_result.n_stations) {
            best_result = best_backward;
        }
    }


    if (ub_ != lb_) {
        for (float alpha:alpha_sched_) {
            for (float beta:beta_sched_) {
                if (alpha == first_alpha && beta == first_beta) {
                    continue;
                }
                reverse_solve_order();//Need to reset the direction
                alpha_ = alpha;
                beta_ = beta;
                ALBPSolution try_forward = solve_one_pass();
                if (try_forward.n_stations < best_result.n_stations) {
                    best_result = try_forward;
                    if (ub_ == lb_) {
                        best_result.optimal = true;
                        break;
                    }
                }
                reverse_solve_order();
                ALBPSolution try_backward = solve_one_pass();
                if (try_backward.n_stations < best_result.n_stations) {
                    best_result = try_backward;
                }
                if (ub_ == lb_) {
                    best_result.optimal = true;
                    break;
                }
            }

        }

    }
    else {
        best_result.optimal = true;
    }

    best_result.method = "MultiHoff";

    return best_result;
}



void MultiHoff::add_new_available(std::vector<int> &eligible_tasks, const int task) {
    if (back_pass_) {
        for (const int j : dir_pred_[task]) {
            n_suc_[j] --;
            if (n_suc_[j] == 0) {
                eligible_tasks.push_back(j);
            }
        }
    }
    else {
        for (const int j : dir_suc_[task]) {
            n_prec_[j] --;
            if (n_prec_[j] == 0) {
                eligible_tasks.push_back(j);
            }
        }
    }


}

void MultiHoff::remove_new_available(std::vector<int> &eligible_tasks, const int task) {
    if (back_pass_) {
        for (const int j : dir_pred_[task]) {
            if (n_suc_[j] == 0) {
                //Start search from the back because it should be in the back
                auto it = std::find(eligible_tasks.rbegin(), eligible_tasks.rend(), j);
                if (it != eligible_tasks.rend()) {
                    // Found it! Convert to forward iterator if needed:
                    auto forward_it = it.base() - 1;
                    std::swap(*forward_it, eligible_tasks.back());
                    eligible_tasks.pop_back();
                }
            }
            n_suc_[j]++;
        }
    }
    else {
        for (const int j : dir_suc_[task]) {
            if (n_prec_[j] == 0) {
                auto it = std::find(eligible_tasks.rbegin(), eligible_tasks.rend(), j);

                if (it != eligible_tasks.rend()) {
                    // Found it! Convert to forward iterator if needed:
                    auto forward_it = it.base() - 1;
                    std::swap(*forward_it, eligible_tasks.back());
                    eligible_tasks.pop_back();
                }

            }
            n_prec_[j]++;
        }
    }

}
void MultiHoff::gen_load( int depth, int remaining_capacity,const int start, float cost, std::vector<int>eligible_tasks) {
    std::uniform_real_distribution<double> dis(0.0, 1.0);

    // 4. Generate the random number

        int full_load = 1;
        for(int i=start;i<eligible_tasks.size();i++) {
            if ((n_attempts_ >= max_attempts_) || (remaining_capacity==0)) return;
            if (int task = eligible_tasks[i]; albp_.task_time[task] <= remaining_capacity) {
                full_load = 0;
                s_task_assign_.push_back(task);
                n_prec_[task] = -1;
                add_new_available(eligible_tasks, task);
                int sub_remaining_capacity = remaining_capacity - albp_.task_time[task];
                float sub_cost = cost - albp_.task_time[task];
                auto random_num = static_cast<float>((dis(rng_)));
                if (back_pass_) { //Update costs with weighted values alpha_ and beta_. Alpha beta are determined at start of heuristic pass

                    sub_cost = sub_cost - alpha_ * static_cast<float>(back_ranking_[task]) - beta_ * static_cast<float>(albp_.pred[task].size())  + gamma_ *random_num * pert_size_ ;
                }
                else {
                    sub_cost = sub_cost - alpha_ * static_cast<float>(forw_ranking_[task]) - beta_ * static_cast<float>(albp_.suc[task].size()) + gamma_ *random_num * pert_size_;

                }
                if (sub_cost < min_cost_) {
                    min_cost_ = sub_cost;
                    best_s_task_assign_ = s_task_assign_;
                }
                gen_load(depth+1, sub_remaining_capacity, i+1, sub_cost,eligible_tasks );

                //undo the changes
                s_task_assign_.pop_back();
                n_prec_[task] = 0;
                remove_new_available(eligible_tasks, task);

            }
        }
    n_attempts_+=full_load;
};

inline void remove_tasks_unordered(std::vector<int>& vec, const std::vector<int>& to_remove) {
    for (int task : to_remove) {
        auto it = std::find(vec.begin(), vec.end(), task);
        if (it != vec.end()) {
            *it = vec.back();
            vec.pop_back();
        }
    }
}

ALBPSolution mhh_solve(const ALBP &albp, const std::optional<std::vector<float>> &alpha_schedule, const std::optional<std::vector<float>> &
                       beta_schedule,const std::optional<float> gamma, const std::optional<std::vector<int>> &task_priorities, const std::optional<unsigned int> seed) {

    MultiHoff mhh= MultiHoff(albp, 5000, alpha_schedule, beta_schedule,gamma, task_priorities, seed);
    ALBPSolution best_result =mhh.solve();






    return best_result;
}
ALBPSolution mhh_solve_salbp1(const ALBP &albp, const std::optional<std::vector<float>> &alpha_schedule, const std::optional<std::vector<float>> &
                              beta_schedule, const std::optional<float> gamma , const std::optional<std::vector<int>> &task_priorities, const std::optional<unsigned> seed) {

    ALBPSolution best_result  = mhh_solve(albp, alpha_schedule, beta_schedule,gamma, task_priorities, seed);
    return best_result;
}

ALBPSolution mhh_solve_salbp1(const int C, const int N, const std::vector<int>& task_times, const std::vector<std::vector<int>>& raw_precedence, const
                              std::optional<std::vector<float>> &alpha_schedule, const std::optional<std::vector<float>> &beta_schedule, const std::optional<float> gamma , const std::
                              optional<std::vector<int>> &task_priorities, const std::optional<unsigned> seed) {
    // print_vector(alpha_schedule.value_or(std::vector<float>{0.0f}));
    // print_vector(alpha_schedule.value_or(std::vector<float>{0.0f}));
    ALBP albp = ALBP::type_1(C, N, task_times, raw_precedence);
    ALBPSolution best_result =mhh_solve(albp, alpha_schedule, beta_schedule,gamma, task_priorities, seed);
    return best_result;
}

void swap_and_pop(int item, std::vector<int> &vec) {
    auto it = std::find(vec.begin(), vec.end(), item);
    if (it != vec.end()) {
        *it = vec.back();      // Replace with last element
        vec.pop_back();        // Remove last element - O(1)
    }
}