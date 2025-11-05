//
// Created by Joseph Thompson on 2025-06-12.
//
#include "albp_solution.h"
#include "ALBP.h"
#include "MultiHoff.h"
#include <iostream>
#include <cfloat>
#include <algorithm>
#include "salbp_basics.h"
#include <map>
#include <numeric>
MultiHoff::MultiHoff(const ALBP& albp, const int max_attempts):
    albp_(albp),
    n_prec_(albp.N, 0),
    n_suc_(albp.N, 0),
    n_prec_orig_(albp.N, 0),
    n_suc_orig_(albp.N, 0),
    n_attempts_(0),
    max_attempts_(max_attempts),// Size = num_tasks, init to 0
    min_cost_(FLT_MAX),
    remaining_task_times_(albp.task_time), //Will zero out task times if task is assigned
    dir_pred_(albp.dir_pred),
    dir_suc_(albp.dir_suc),
    mhh_sol_(albp.N),
    // Fill with 0 to N-1


    ub_(albp.N),
    // Initialize to 0
     lb_(calc_salbp_1_bin_lbs(albp.task_time, albp.C))
    {


    s_task_assign_.reserve(albp.N),


    pw_ranking_ = pw_ranking(albp_) ;
    rpw_ranking_ = rpw_ranking(albp_) ;
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
            sort_by_ranking(eligible_tasks_forward, pw_ranking_);
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
            sort_by_ranking(eligible_tasks_backward, rpw_ranking_);
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
        std::cout << "saving best solution, n stations " <<ub_ << std::endl;
        if (reverse_) {
            std::reverse(mhh_sol_.task_assignment.begin(), mhh_sol_.task_assignment.end());
        }
        mhh_sol_.n_stations = ub_;
        mhh_sol_.task_to_station_and_load(albp_);
        mhh_sol_.station_to_ranking();
        return mhh_sol_;
    }
        return mhh_sol_;

}


void MultiHoff::reverse_solve_order() {
    /* This function reverses the elements necessary for a backwards solve. Note that it modifies some internal state*/
    reverse_ = true;
    back_pass_ = false;
    std::swap(dir_pred_, dir_suc_);
    std::swap(n_prec_orig_, n_suc_orig_);
    std::swap(pw_ranking_, rpw_ranking_);
    std::swap(no_prec_tasks_, no_suc_tasks_);
    initialize_current_s_assignments();
    s_forwards_.clear();
    s_backwards_.clear();
    remaining_task_times_ = albp_.task_time;



}


ALBPSolution MultiHoff::solve() {
    std::cout << "\n Starting forward solve..." << std::endl;

    ALBPSolution best_result = solve_one_pass();
    std::cout << " \n \n Forward solve complete. UB=" << ub_ << " LB=" << lb_ << std::endl;

    if (ub_ != lb_) {
        std::cout << "Starting reverse solve..." << std::endl;

        reverse_solve_order();
        std::cout << "Reverse order setup complete" << std::endl;

        ALBPSolution best_backward =solve_one_pass();
        std::cout << "Reverse solve complete" << std::endl;

        if (best_backward.n_stations < best_result.n_stations) {
            best_result = best_backward;
        }
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

        int full_load = 1;
        for(int i=start;i<eligible_tasks.size();i++) {
            if ((n_attempts_ >= max_attempts_) || (remaining_capacity==0)) return;
            if (int task = eligible_tasks[i]; albp_.task_time[task] <= remaining_capacity) {
                full_load = 0;
                s_task_assign_.push_back(task);
                n_prec_[task] = -1;
                add_new_available(eligible_tasks, task);
                int sub_remaining_capacity = remaining_capacity - albp_.task_time[task];
                float sub_cost = cost - albp_.task_time[task]; //- alpha_ * pw_[task] - beta_ * albp_.suc[task].size();
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

ALBPSolution mhh_solve(const ALBP &albp ) {

    MultiHoff mhh= MultiHoff(albp);
    ALBPSolution best_result =mhh.solve();






    return best_result;
}
ALBPSolution mhh_solve_salbp1(const ALBP &albp ) {

    ALBPSolution best_result  = mhh_solve(albp);
    return best_result;
}

ALBPSolution mhh_solve_salbp1(const int C,const int N, const std::vector<int>& task_times, const std::vector<std::vector<int>>& raw_precedence) {

    ALBP albp = ALBP::type_1(C, N, task_times, raw_precedence);
    ALBPSolution best_result =mhh_solve(albp);
    return best_result;
}

inline void swap_and_pop(int item, std::vector<int> &vec) {
    auto it = std::find(vec.begin(), vec.end(), item);
    if (it != vec.end()) {
        *it = vec.back();      // Replace with last element
        vec.pop_back();        // Remove last element - O(1)
    }
}