//
// Created by Joseph Thompson on 2025-07-02.
//

#include "vdls.h"

#include <iostream>
#include <random>
#include <optional>
#include <algorithm>
#include <ranges>
#include "Hoff.h"
#include "MultiHoff.h"
#include "salbp_basics.h"


struct Profile {
        std::chrono::milliseconds local_search_time{0};
        std::chrono::milliseconds perturbation_time{0};
        std::chrono::milliseconds hoff_search_time{0};
        int local_search_calls = 0;
        int perturbations = 0;
        int improvements = 0;
};
Profile PROFILE;

bool VDLS::time_exceeded() const {
        auto now = std::chrono::steady_clock::now();
        bool time_exceeded = (now - start_time_) >= time_limit_;
        return time_exceeded;
}
void VDLS::add_init_solution(std::vector<int>init_solution) {
        std::cout<<"Initial solution added, processing"<<std::endl;
        best_ = process_init_solution(albp_, init_solution);
}
ALBPSolution VDLS::solve_type_1(  ) {
        if (best_.station_assignments.empty()) {
                best_= mhh_solve_salbp1(albp_); //Get initial SALBP-1 solution
                best_.method = "VDLS: hoff start";

        }
        int n_stations = best_.n_stations;
        const int salbp_1_lb = calc_salbp_1_lbs(albp_);
        std::cout<<"best hoff solution "<< best_.n_stations<< " lb: "<< salbp_1_lb<<std::endl;
        while (best_.n_stations > salbp_1_lb && !time_exceeded() ) {

                n_stations --;//Try again with one fewer station
                //Check to see if even possible in allotted cycle time
                lb_ = calc_salbp_2_lbs(albp_.task_time, n_stations);
                if (lb_ > albp_.C) {
                        std::cout << "Calculated lb: " << lb_ << " Impossible to fit " << n_stations << std::endl;
                        break;
                }


                ALBPSolution local_best = vdls_heuristic(n_stations, albp_.C);
                local_best.method = "VDLS: local search";
                if (local_best.cycle_time <= albp_.C) {
                        best_ = local_best;


                }
                else break; //Couldn't find an SALBP-2 solution <C for the given number of stations. Give up.
        }
        return best_;
}

ALBPSolution VDLS::solve_type_2(  ) {
        int lb = calc_salbp_2_lbs(albp_.task_time, albp_.S);
        lb_ = lb; //Had to feed in separate lower bound for SALBP-1, that is why it looks goofy on this one
        best_ = vdls_heuristic(albp_.S, lb);
        best_.method = "VDLS_SALBP2";
        return best_;
}

ALBPSolution VDLS::hoff_search(int n_stations) const {
        ALBP albp = albp_;//Copying problem so we can modify the cycle time
        int lb = calc_salbp_2_lbs(albp_.task_time, n_stations);
        int ub = calc_salbp_2_ub(albp_.task_time, n_stations);
        int a0 = 0;
        int a1 = 1;
        ALBPSolution test_sol{albp.N};
        ALBPSolution best_sol = test_sol;
        do  {
                albp.C = lb +a0 ;
                test_sol = hoff_solve_salbp1(albp);
                if (test_sol.n_stations > n_stations) {
                        lb += a0 +1 ;
                        const int a2 = a0+a1;
                        a0 = a1;
                        a1 = a2;
                }
                else {
                        ub =lb + a0;
                        a0 = 0;
                        a1 = 1;
                        best_sol = test_sol;
                        // if (lb+1==ub) {
                        //         lb=ub;
                        // }
                }



        }
        while ((lb < ub));

        return best_sol;
}

/* Moves a task to a new station and updates the solution**/
void VDLS::perform_shift(ALBPSolution &sol, int task, int task_idx, int old_station, int new_station) {
        //Changes the task assignments
        sol.task_assignment[task] = new_station;
        //changes stations assignments
        sol.station_assignments[old_station][task_idx] = sol.station_assignments[old_station].back();  // Move last element to the one being removed
        sol.station_assignments[old_station].pop_back();// Remove last element
        sol.station_assignments[new_station].push_back(task);
        //Changes the station loads
        int old_load = sol.loads[old_station];
        sol.loads[old_station] -= albp_.task_time[task];
        sol.loads[new_station] += albp_.task_time[task];
        //Update the cycle time (if applicable)

        if (sol.loads[new_station] > sol.cycle_time) {
                sol.cycle_time = sol.loads[new_station];
                sol.critical_stations = {new_station};
        }
        else if (sol.loads[new_station] == sol.cycle_time){
                std::replace(sol.critical_stations.begin(), sol.critical_stations.end(), old_station, new_station);

        }
        //Old station was the maximum, and new station did not reach the maximum, improvement
        else if (old_load >= sol.cycle_time) {
                if (sol.critical_stations.size() == 1) {
                        auto updates = get_critical_stations(sol.loads);
                        sol.cycle_time =updates.first;
                        sol.critical_stations = updates.second;
                }
                else {
                       //There are other critical stations, remove old one
                        swap_and_pop(old_station, sol.critical_stations);
                }

        }
        //Changes the earliest and latest for parents and children
        sol.update_window(albp_,task, old_station);

}
/*Recursive DFS algorithm for exploring different task shifts up to a given depth. Returns true if there was a local
 * improvement to the solution
 */
bool VDLS::local_search(const ALBPSolution &incumbent, ALBPSolution &local_best, int depth, int previous_task) {
        if (depth ==max_depth_ || incumbent.cycle_time <= lb_ || time_exceeded()) {
                return false;
        }
        //Reassign tasks to other stations, provided it wasn't previous task (DFS)
        std::vector<int> to_check = {incumbent.critical_stations.back()}; //Only look at previous move for depths greater than 1
        if (depth==0) to_check = incumbent.critical_stations;
        //Start with fresher critical stations at back
        for (int old_station : std::ranges::reverse_view(to_check)) {
                for (int task_idx=0; task_idx < incumbent.station_assignments[old_station].size(); task_idx++){
                        int task = incumbent.station_assignments[old_station][task_idx];
                        if (task != previous_task) {
                                for (int new_station = incumbent.earliest[task]; new_station <= incumbent.latest[task]; new_station++ ) {
                                        if( new_station != (old_station)){
                                                if (time_exceeded()) return false;
                                                ALBPSolution new_sol = incumbent;
                                                perform_shift(new_sol, task, task_idx, old_station, new_station);
                                                //std::cout<<"Before I moved "<< task<< " from " << old_station << "  to station: "<<  new_station <<  " depth " << depth << std::endl;
                                                if ((new_sol.critical_stations.size() < local_best.critical_stations.size() && new_sol.cycle_time == local_best.cycle_time) || new_sol.cycle_time < local_best.cycle_time) {
                                                        // std::cout<<"IMPROVEMENT moved "<< task<< " from " << old_station << "  to station: "<<  new_station <<  " depth " << depth << " new station load  " << new_sol.loads[new_station] <<  " old loads " << new_sol.loads[old_station] << " new sol cycle time " << new_sol.cycle_time  << " best cycle time " <<  local_best.cycle_time<<std::endl;
                                                        local_best = new_sol;
                                                        return true;

                                                }
                                                //Go to next level, check for improvement
                                                if (local_search(new_sol, local_best,depth + 1, task )) {
                                                        return true;
                                                }


                                                //         if (improved ==false && time_exceeded()) return false;
                                        }
                                }
                        }
                }
        }

        return false;
}

std::vector<int> get_max_indices(const std::vector<int> & station_load, const std::vector<int> &can_change) {

        std::vector<int> max_indices;
        int max_val = -1;
        for (int i = 0; i < station_load.size(); ++i) {
                if (can_change[i] == 0) continue;
                if (station_load[i] == max_val) {

                        max_indices.push_back(i);
                }
                else if (station_load[i] > max_val) {
                        max_val = station_load[i];
                        max_indices = {static_cast<int>(i)};
                }
        }
        return max_indices;
}

/*randomly selects from a vector of integers */
int VDLS::random_selection(const std::vector<int>& int_vec) {
        std::uniform_int_distribution<> dist(0, int_vec.size() - 1);
        return int_vec[dist(rng_)];
}


std::vector<int> range_excluding(int start, int end, int exclude) {
        std::vector<int> result;
        for (int i = start; i <= end; ++i) {
                if (i != exclude) {
                        result.push_back(i);
                }
        }
        return result;
}

int select_new_station(ALBPSolution &sol, int task, int station_idx) {
        int earliest = sol.earliest[task];
        int latest = sol.latest[task];
        //Move task to a station that is not the current station
        std::vector<int> station_range = range_excluding(earliest, latest, station_idx);
        if (station_range.empty()) {
                return -1;
        }
        int selected_station = station_range[0];
        for (int i = 1; i < station_range.size(); ++i) {
                if (sol.loads[station_range[i]] < sol.loads[selected_station]) {
                        selected_station = station_range[i];
                }

        }
        return selected_station;
}

std::vector<int> check_windows(ALBPSolution &sol, int station) {
        std::vector<int> task_indices;
        for (int i = 0; i < sol.station_assignments[station].size(); i++) {
                int task = sol.station_assignments[station][i];
                int left = sol.earliest[task];
                int right = sol.latest[task];
                if (left != right) {
                        task_indices.push_back(i);
                }
        }
        return task_indices;
}

/* Task based perturbation */
void VDLS::perturbation(ALBPSolution &current_sol) {
        std:: vector<int> can_change(current_sol.n_stations, 1 );//keeps track of perturbed stations
        for (int i=0; i<n_perts_; i++) {
                std::vector<int> max_indices = get_max_indices(current_sol.loads, can_change);;
                //if there are no more max indices, that means that we have ran out of changeable stations
                if (max_indices.empty()) {
                        return;
                }
                //Select a station to remove tasks from
                const int station_idx = random_selection(max_indices);
                std::vector<int> filtered_tasks = check_windows(current_sol, station_idx);
                //Check for station feasibility:If no tasks can be moved, try different station
                if (filtered_tasks.empty()) {
                        can_change[station_idx] = 0;
                        continue;
                }
                //Select task
                const int task_idx = random_selection(filtered_tasks);
                const int task = current_sol.station_assignments[station_idx][task_idx];
                //select new station to move to
                const int new_station = select_new_station(current_sol, task, station_idx);
                if (new_station != -1) {
                        perform_shift(current_sol, task, task_idx, station_idx, new_station);
                        can_change[new_station] = 0;
                }

        }

}


ALBPSolution VDLS::vdls_heuristic( int n_stations,  int lb) {
        auto start_hoff = std::chrono::steady_clock::now();
        ALBPSolution current_solution = hoff_search(n_stations);
        ALBPSolution best_sol = current_solution;
        while ((best_sol.cycle_time > lb) &&(!time_exceeded())) {
                //Keep running local search until it stops improving
                ALBPSolution local_best = current_solution;
                while (local_search(current_solution, local_best,0, -1)){
                       // current_solution.print();
                        current_solution = local_best;
                        if ( local_best.cycle_time < best_sol.cycle_time ) {
                                best_sol = local_best;
                                best_sol.method = "VDLS (local search)";
                                auto now = std::chrono::steady_clock::now();
                                best_sol.elapsed_ms = std::chrono::duration_cast<std::chrono::milliseconds>(now - start_time_);
                                if (best_sol.cycle_time == lb) {
                                        best_sol.optimal = true;
                                        return best_sol;
                                }
                        }
                }
                std::cout <<" perturbation! "<< best_sol.elapsed_ms.count() << " ms" << std::endl;
                perturbation(current_solution);
        }
        return best_sol;
}


ALBPSolution vdls_solve_salbp1(const ALBP &albp, const std::vector<int> &initial_solution, std::optional<int> max_attempts, std::optional<double> time_limit,
        std::optional<unsigned> seed ) {
        int attempts = max_attempts.value_or(100000);  // default if not passed
        double limit = time_limit.value_or(7200.);

        auto vdls= VDLS(albp, limit, seed);
        if (!initial_solution.empty()) {
                vdls.add_init_solution(initial_solution);
        }
        ALBPSolution result =vdls.solve_type_1();
        return result;
}

ALBPSolution vdls_solve_salbp2(const ALBP &albp, const std::vector<int> &initial_solution, std::optional<int> max_attempts, std::optional<double>
                               time_limit, std::optional<unsigned> seed ) {
        int attempts = max_attempts.value_or(5000);  // default if not passed
        double limit = time_limit.value_or(7200.0);

        auto vdls= VDLS(albp,  limit, seed);
        if (!initial_solution.empty()) {
                vdls.add_init_solution(initial_solution);
        }
        ALBPSolution result =vdls.solve_type_2();
        return result;
}
ALBPSolution vdls_solve_salbp2(const int S, const int N, const std::vector<int>& task_times, const std::vector<std::vector<int>>& raw_precedence, const std::vector<int> &initial_solution, std::optional<int> max_attempts,
                               std::optional<double> time_limit, std::optional<unsigned> seed ) {

        ALBP albp = ALBP::type_2(S, N, task_times, raw_precedence);

        return vdls_solve_salbp2(albp, initial_solution, max_attempts, time_limit, seed);
}

ALBPSolution vdls_solve_salbp1(const int C, const int N, const std::vector<int>& task_times, const std::vector<std::vector<int>>& raw_precedence, const std::vector<int> &initial_solution, std::optional<int> max_attempts,
                               std::optional<double> time_limit, std::optional<unsigned> seed ) {
        ALBP albp = ALBP::type_1(C, N, task_times, raw_precedence);
        return vdls_solve_salbp1(albp, initial_solution, max_attempts, time_limit, seed);
}
