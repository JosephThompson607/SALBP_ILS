//
// Created by Joseph Thompson on 2026-06-30.
//

#include "tabu.h"

#include <algorithm>
#include <deque>
#include <filesystem>

#include "../ALBP.h"
#include "../albp_solution.h"
#include <iostream>
#include <random>

#include "Hoff.h"
#include "salbp_basics.h"
//TABU LIST MANAGEMENT
TabuList::TabuList(int n, int max_size) {
    n_ = n;
    max_size_ = max_size;
    is_tabu_.resize(n *n , -1);
}


int TabuList::tabu_status(const std::pair<int, int> &move) const {
    return is_tabu_[index(move)];
}

bool TabuList::del_move(const std::pair<int, int> &move) {
    is_tabu_[index(move)] = false;
    auto it = std::find(tabu_list_.begin(), tabu_list_.end(), move);
    if (it != tabu_list_.end()) {
        tabu_list_.erase(it);
        return true; // found and removed
    }
    return false; // not found
}

void TabuList::check_tenure() {
    if (tabu_list_.size() ==0) return;
    std::pair first = tabu_list_.front();
    while (tabu_status(first) == iteration_ % max_size_) {
        is_tabu_[index(first)] = -1;
        tabu_list_.pop_front();
        first = tabu_list_.front();
    }

}

void TabuList::step() {
    increment_counter();
    check_tenure();
}

void TabuList::increment_counter() {
    iteration_ = (iteration_ + 1) % max_size_;
}

int TabuList::get_counter() const {
return iteration_;
}


void TabuList::print() const {
    std::cout << "iteration " << iteration_<< ": Tabu list (" << tabu_list_.size() << "/" << max_size_ << "): ";
    for (const auto &mv : tabu_list_) {
        std::cout << "(" << mv.first << "->" << mv.second <<", t: " << is_tabu_[index(mv)]<< ") " ;
    }
    std::cout << std::endl;
}


void TabuList::insert(const std::pair<int, int> &move) {
    int idx = index(move);
    if (is_tabu_[idx]>=0) return;//If is already tabu, stop from inserting
    //Add move to tabu list
    tabu_list_.push_back(move);
    is_tabu_[idx] = iteration_ % max_size_;
}
void TabuList::reset() {
    std::fill(is_tabu_.begin(), is_tabu_.end(), -1);
    tabu_list_.clear();
    iteration_ = 0;
}

int TabuList::index(const std::pair<int, int> &move) const {
    return static_cast<int>(move.first) * n_ + move.second;
}


//TABU Search




void Tabu::add_init_solution(std::vector<int>init_solution) {
    std::cout<<"Initial solution added, processing"<<std::endl;
    best_ = process_init_solution(albp_, init_solution);

}



int Tabu::try_shift(const ALBPSolution &sol, int task,  int old_station,int new_station) const {
    int obj = sol.tot_ovlo2_;
    //Get rid of overloads that will change
    obj -= sol.sq_overloads_[old_station];
    obj -= sol.sq_overloads_[new_station];
    //recalculate overloads
    int ovlo1 =  std::max(sol.loads[new_station] + albp_.task_time[task] - albp_.C,0);
    int ovlo2 =  std::max(sol.loads[old_station] - albp_.task_time[task] - albp_.C,0);
    obj += ovlo1*ovlo1 + ovlo2*ovlo2;
    return obj;
}

int Tabu::try_swap(const ALBPSolution &sol, int task_1,  int task_2,  int station_1,
    int station_2) const {
    int obj = sol.tot_ovlo2_;
    //Get rid of overloads that will change
    obj -= sol.sq_overloads_[station_1];
    obj -= sol.sq_overloads_[station_2];
    //recalculate overloads
    int ovlo1 =  std::max(sol.loads[station_2] + albp_.task_time[task_1] - albp_.task_time[task_2] -albp_.C,0);
    int ovlo2 =  std::max(sol.loads[station_1] - albp_.task_time[task_1] + albp_.task_time[task_2] -albp_.C,0);
    obj += ovlo1*ovlo1 + ovlo2*ovlo2;
    return obj;
}

void Tabu::perform_shift(ALBPSolution &sol, int task, int task_idx, int old_station,int new_station) {
    //Changes the task assignments
    sol.task_assignment[task] = new_station;
    //changes stations assignments
    sol.station_assignments[old_station][task_idx] = sol.station_assignments[old_station].back();  // Move last element to the one being removed
    sol.station_assignments[old_station].pop_back();// Remove last element
    sol.station_assignments[new_station].push_back(task);
    //Changes the station loads
    sol.loads[old_station] -= albp_.task_time[task];
    sol.loads[new_station] += albp_.task_time[task];
    sol.tot_ovlo2_ -= sol.sq_overloads_[old_station];
    sol.tot_ovlo2_ -= sol.sq_overloads_[new_station];
    sol.overloads_[new_station] =  std::max(sol.loads[new_station] - albp_.C,0);
    sol.overloads_[old_station] =  std::max(sol.loads[old_station]- albp_.C,0);
    sol.sq_overloads_[new_station] = sol.overloads_[new_station]*sol.overloads_[new_station];
    sol.sq_overloads_[old_station] = sol.overloads_[old_station]*sol.overloads_[old_station];
    sol.tot_ovlo2_ += sol.sq_overloads_[new_station];
    sol.tot_ovlo2_ += sol.sq_overloads_[old_station];
    //Update the cycle time (if applicable)
    auto it = std::max_element(sol.loads.begin(), sol.loads.end());
    if (it != sol.loads.end()) {
        size_t index = std::distance(sol.loads.begin(), it);
        auto value = *it;
        sol.cycle_time = value;
        sol.max_ovlo_ = std::max(value-albp_.C,0);
        max_ovlo_station_ = index;
    }
    //Changes the earliest and latest for parents and children

    sol.update_window(albp_,task);




}

// void Tabu::perform_swap(ALBPSolution &sol, int task_1, int task_idx_1, int task_2, int task_idx_2, int station_1,
//     int station_2) {
//     perform_shift(sol, task_1, task_idx_1, station_1, station_2);
//     perform_shift(sol, task_2, task_idx_2, station_2, station_1);
// }

void Tabu::quad_overloads(ALBPSolution &sol, const int C) {
    sol.overloads_.resize(sol.loads.size());
    sol.sq_overloads_.resize(sol.loads.size());
    max_ovlo_station_ = 0;
    sol.max_ovlo_ = 0;
    sol.tot_ovlo2_ = 0;
    for (int i = 0; i < sol.loads.size(); ++i) {
        int ove1 = std::max( sol.loads[i]-C, 0);
        int sq_ove = ove1 * ove1;
        sol.sq_overloads_[i] = sq_ove;
        sol.overloads_[i] = ove1;
        sol.tot_ovlo2_ += sol.sq_overloads_[i];
        if (sol.overloads_[i] > sol.max_ovlo_) {
            max_ovlo_station_ = i;
            sol.max_ovlo_ = sol.overloads_[i];
        }
    }
}

void Tabu::elim_station(ALBPSolution &sol) {
    //eliminates smallest staiton, not respecting cycle time
    auto it = std::min_element(sol.loads.begin(), sol.loads.end());
    const long s = std::distance(sol.loads.begin(), it);
    ALBPSolution current_sol = sol;
    //Make sure that the station gets emptied
    while (current_sol.station_assignments[s].size() >0) {
        for (int i=0; i < sol.station_assignments[s].size(); i++) {
            int task = sol.station_assignments[s][i];
            if (current_sol.task_assignment[task] != s) continue; //If task was reassigned already, move on
            int target_c = std::numeric_limits<int>::max(); //Create a large objective value for minimization
            auto &tasks_at_s = current_sol.station_assignments[s];
            auto it = std::find(tasks_at_s.begin(), tasks_at_s.end(), task);
            int cur_idx = static_cast<int>(std::distance(tasks_at_s.begin(), it));


            auto local_best = TabuMove(current_sol, target_c, std::pair(-1,-1), std::pair(-1,-1), -1);
            for (int j = sol.earliest[task]; j <= sol.latest[task]; j++ ) {
                if( j != (s)){
                    int obj = try_shift(current_sol, task, s, j);
                    if (obj < target_c) {
                        target_c = obj;
                        local_best = TabuMove(current_sol, obj, std::pair(task,s), std::pair(task,j), cur_idx);
                    }
                }
            }
            do_move(current_sol, local_best, false);
        }
    }
    current_sol.n_stations = sol.n_stations-1;
    current_sol.station_assignments.erase(current_sol.station_assignments.begin() + s);
    //removes station
    current_sol.loads.erase(current_sol.loads.begin() + s);
    //Caculates initial overloads
    quad_overloads(current_sol, albp_.C);
    //Fixes task assignment
    current_sol.station_to_task();
    current_sol.find_windows(albp_);
    sol = current_sol;
}

void Tabu::print_move(const TabuMove &move) const {
    std::cout << "=== Move ===" << std::endl;
    std::cout << "Station loads (before move): ";
    for (size_t i = 0; i < move.sol.loads.size(); ++i) {
        std::cout << "[" << i << "]=" << move.sol.loads[i] << " ";
    }
    std::cout << std::endl;

    if (move.move == std::pair<int,int>(-1,-1)) {
        std::cout << "No move selected (local_best empty)" << std::endl;
    } else {
        std::cout << "Shift: task " << move.move.first
                   << " station " << move.reversal.second
                   << " -> station " << move.move.second << std::endl;
        if (move.move2.has_value()) {
            std::cout << "Swap with: task " << move.move2->first
                       << " station " << move.reversal2->second
                       << " -> station " << move.move2->second << std::endl;
        }
        std::cout << "Predicted objective after move: " << move.obj << std::endl;
    }

    tabu_.print();
    std::cout << "============" << std::endl;
}

int Tabu::roulette_select(ALBPSolution&sol) {
    std::discrete_distribution<int> dist(sol.sq_overloads_.begin(), sol.sq_overloads_.end());
    return dist(rng_);
}

bool Tabu::is_optimal(ALBPSolution & sol) {
    if (sol.cycle_time<= albp_.C && sol.n_stations == lb_) {
        sol.optimal= true;
        return true;
    }
    return false;

}

TabuMove Tabu::shift(int s, int & local_obj, const ALBPSolution & current) {
    TabuMove local_best = TabuMove(current, local_obj, std::pair(-1,-1), std::pair(-1,-1), -1);
    for (int i=0; i < current.station_assignments[s].size(); i++) {
        int task = current.station_assignments[s][i];
        for (int j = current.earliest[task]; j <= current.latest[task]; j++ ) {
            if( j != (s)){
                std::pair move = std::pair(task, j);
                std::pair<int,int> reversal = std::pair(task, s);
                int obj = try_shift(current, task, s, j);
                if (obj < current.tot_ovlo2_) {//Aspiration criteria
                    local_obj = obj;
                    return TabuMove(current, obj, reversal,move, i);
                }
                if ( tabu_.tabu_status(move) == -1  && obj <=  local_obj) {
                    local_obj = obj;
                    local_best = TabuMove(current, obj, reversal, move, i);
                }

            }
        }
        }

    return local_best;
}

TabuMove Tabu::swap(int s, int &local_obj, const ALBPSolution &current) {
    TabuMove local_best = TabuMove(current, local_obj, std::pair(-1,-1), std::pair(-1,-1), -1);
    for (int task_idx=0; task_idx < current.station_assignments[s].size(); task_idx++) {
        int task = current.station_assignments[s][task_idx];
        for (int new_station = current.earliest[task]; new_station <= current.latest[task]; new_station++ ) {
            for (int new_task_idx=0; new_task_idx < current.station_assignments[new_station].size(); new_task_idx++){
                int new_task = current.station_assignments[new_station][new_task_idx];
                if ((new_station == (s)) || (new_station<s && albp_.prec_mat[task * albp_.N + new_task]==1) || (new_station>s && albp_.prec_mat[new_task * albp_.N + task]==1) ) {
                    continue;
                }
                if(current.earliest[new_task]<= s && current.latest[new_task]>= s) {
                    std::pair move1 = std::pair(task, new_station);
                    std::pair move2 = std::pair(new_task, s);
                    std::pair<int,int> reversal1 = std::pair(task, s);
                    std::pair<int,int> reversal2 = std::pair(new_task, new_station);
                    int obj = try_swap(current, task, new_task, s, new_station);
                    if (obj < current.tot_ovlo2_) {//Aspiration criteria
                        local_obj = obj;
                        return TabuMove(current, obj, reversal1, move1, task_idx, reversal2, move2,  new_task_idx);
                    }
                    if ( (tabu_.tabu_status(move1)   == -1) && (tabu_.tabu_status(move2)== -1) && obj <=  local_obj) {
                        local_obj = obj;
                        local_best = TabuMove(current, obj, reversal1, move1, task_idx, reversal2,move2,  new_task_idx);
                    }
                    }

                }
            }
        }
    return local_best;
}

void Tabu::do_move(ALBPSolution & sol, const TabuMove &move, const bool update_tabu= true) {

    if (move.move != std::pair<int,int>(-1,-1)) {
        sol = move.sol;
        if (update_tabu) {
            tabu_.insert(move.reversal);

        }
        int task_1 = move.move.first;
        int task_1_idx = move.task_idx;
        int station_2 = move.move.second;
        int station_1 = move.reversal.second;
        perform_shift(sol, task_1, task_1_idx, station_1, station_2);
        //If move was swap, must update info
        if (move.move2.has_value()) {
            if (update_tabu) tabu_.insert(move.reversal2.value());
            int task_2 = move.move2.value().first;
            int task_2_idx = move.task2_idx.value();
            station_1 = move.move2.value().second;
            station_2 = move.reversal2.value().second;
            perform_shift(sol, task_2, task_2_idx, station_2, station_1);
        }


    }



}

void Tabu::shift_and_swap(ALBPSolution& sol, int  s) {
    //Makes the current solution look bad so it gets changed
    int loc_obj = std::numeric_limits<int>::max();
    TabuMove loc_move = TabuMove(sol, loc_obj,  std::pair(-1,-1), std::pair(-1,-1), -1);
    //Check shift
    TabuMove shift_move = shift(s, loc_obj, sol);
    if (shift_move.obj < loc_move.obj){
        loc_move =shift_move;
        }
    //Check swap
    TabuMove swap_move = swap(s, loc_obj, sol);
    if (swap_move.obj < loc_move.obj) {
        loc_move =swap_move;
    }
    //Do the best move
    print_move(loc_move);
    do_move(sol, loc_move, true);
    tabu_.step(); // Clear out tabu list if old moves are being stored



}

ALBPSolution Tabu::solve(){
    start_time_ = std::chrono::steady_clock::now();
    //Get initial SALBP-1 solution if one is not given
    if (best_.station_assignments.empty()) {
        // best_= hoff_solve_salbp1(albp_);
        // best_.method = "hoff(tabu start)";
        best_= priority_solve_salbp_1(albp_, 4, 42);
        best_.method = "priority ranking(tabu start)";
    }
    best_.overloads_.resize(best_.station_assignments.size(), 0);
    best_.sq_overloads_.resize(best_.station_assignments.size(), 0);
    current_ = best_;
    //Makes sure best has not been obtained from initial solution
    if (is_optimal(best_)) return best_;

    //Remove station with lowest load until we reach a cycle time violation
    tabu_.reset();
    while (current_.cycle_time <= albp_.C && !time_exceeded()){
        elim_station(current_);
        if (current_.cycle_time <= albp_.C) {
            best_ = current_;
            best_.method = "tabu: pre-elimination";
            if (is_optimal(best_)) return best_;
        }
    }




    while (!time_exceeded() && !(is_optimal(best_))) {
        std::bernoulli_distribution coin(0.9);
        int s = 0;
        if (coin(rng_)) {
            s= roulette_select(current_);
        } else {            // skip it (10% of calls)
            std::uniform_int_distribution<size_t> dist(0, current_.n_stations-1);
            s = dist(rng_);

        }
        shift_and_swap(current_, s);
        if (current_.cycle_time <= albp_.C){
            best_ = current_;
            best_.method = "tabu: main";
            if (is_optimal(best_)) return best_;
            elim_station(current_);
            }
    }
    return best_;




}

bool Tabu::time_exceeded() const {
    const auto now = std::chrono::steady_clock::now();
    return (now - start_time_) >= time_limit_;
}



//Main function calls
ALBPSolution tabu_solve_salbp1(const ALBP &albp , std::optional<double> time_limit) {
    double limit = time_limit.value_or(3.0);
    //creates tabu instance, (with calculating BPP lbs)
    Tabu tab= Tabu(albp, limit);
    ALBPSolution best_result =tab.solve();
    return best_result;
}

ALBPSolution tabu_solve_salbp1(int C, int N,
             const std::vector<int> &task_times,
             const std::vector<std::vector<int> > &raw_precedence, std::optional<double> time_limit) {
    //creates tabu instance, (with calculating BPP lbs)
    ALBP albp = ALBP::type_1(C, N, task_times, raw_precedence);
    return tabu_solve_salbp1(albp, time_limit);
}
