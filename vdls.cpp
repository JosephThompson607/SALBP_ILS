//
// Created by Joseph Thompson on 2025-07-02.
//

#include "vdls.h"

#include <iostream>
#include <random>
#include <__ostream/basic_ostream.h>

#include "mhh.h"

bool VDLS::time_exceeded() const {
        auto now = std::chrono::steady_clock::now();
        return (now - start_time_) >= time_limit_;
}
void VDLS::add_init_solution(std::vector<int>init_solution) {
        std::cout<<"Initial solution added, processing"<<std::endl;
        best_ = process_init_solution(albp_, init_solution);
}
ALBPSolution VDLS::solve_type_1(  ) {
        if (best_.station_assignments.empty()) {
                std::cout<<"no initial solution, calculating new solution"<<std::endl;
                best_= mhh_solve_salbp1(albp_); //Get initial SALBP-1 solution
        }
        int n_stations = best_.n_stations;
        const int salbp_1_lb = calc_lb_1(albp_.task_time, albp_.C); //TODO: add in other lower bounds
        while (best_.n_stations > salbp_1_lb && !time_exceeded() ) {
                ALBPSolution candidate = best_;
                n_stations --;//Try again with one fewer station
                //Check to see if even possible in allotted cycle time
                lb_ = calc_salbp_2_lbs(albp_.task_time, candidate.n_stations);
                if (lb_ > albp_.C) {
                        std::cout << "Calculated lb: " << lb_ << " Impossible to fit " << candidate.n_stations << std::endl;
                        break;
                }


                ALBPSolution local_best = vdls_heuristic(n_stations, albp_.C);
                if (local_best.cycle_time <= albp_.C) {
                        best_ = local_best;


                }
                else break; //Couldn't find an SALBP-2 solution <C for the given number of stations. Give up.
        }
        return best_;
}

ALBPSolution VDLS::solve_type_2(  ) {
        best_ = vdls_heuristic(albp_.S, albp_.C);
        return best_;
}

ALBPSolution VDLS::hoff_search(int n_stations) {
        ALBP albp = albp_;//Copying problem so we can modify the cycle time
        int lb = calc_salbp_2_lbs(albp_.task_time, n_stations);
        int ub = calc_salbp_2_ub(albp_.task_time, n_stations);
        int a0 = 0;
        int a1 = 1;
        ALBPSolution test_sol{albp.N};
        while (lb < ub) {
                albp.C = lb;
                test_sol = hoff_solve_salbp1(albp);
                if (test_sol.n_stations > n_stations) {
                        lb += a0;
                        const int a2 = a0+a1;
                        a0 = a1;

                        a1 = a2;
                }
                else {


                        ub =lb + a0;
                        a0 = 0;
                        a1 = 1;
                }



        }


        return test_sol;
}
/* Moves a task to a new station and updates the solution**/
void VDLS::perform_shift(ALBPSolution &sol, int task, int task_idx, int old_station,int new_station) {
        //Changes the task assignments
        sol.task_assignment[task] = new_station;
        //changes stations assignments
        sol.station_assignments[old_station][task_idx] = sol.station_assignments[old_station].back();  // Move last element to the one being removed
        sol.station_assignments[old_station].pop_back();// Remove last element
        sol.station_assignments[new_station].push_back(task);
        //Changes the station loads
        sol.load[old_station] += albp_.task_time[task];
        sol.load[new_station] -= albp_.task_time[task];
        //Changes the earliest and latest for parents and children
        sol.update_window(albp_,task);




}
/*Recursive DFS algorithm for exploring different task shifts up to a given depth. Returns true if there was a local
 * improvement to the solution
 */
bool VDLS::local_search(ALBPSolution &local_best, ALBPSolution &incumbent, int depth, int last_task, bool improved) {
        if (depth ==max_depth_ || local_best.cycle_time <= lb_) {
                return false;
        }
        //Get the station with the highest load
        auto max_it = std::max_element(incumbent.load.begin(),incumbent.load.end());
        int station_idx = std::distance(incumbent.load.begin(),max_it);
        //Reassign tasks to other stations, provided it wasn't previous task (DFS)
        for (int i=0; i < incumbent.station_assignments[station_idx].size(); i++){
                int task = incumbent.station_assignments[station_idx][i];
                if (task != last_task) {
                        for (int j = incumbent.earliest[task]; j < incumbent.latest[task]; j++ ) {
                                if( j != (station_idx)){
                                        ALBPSolution new_sol = incumbent;
                                        perform_shift(new_sol, task, i, station_idx, j);
                                        if (new_sol.cycle_time < local_best.cycle_time) {
                                                local_best = new_sol;
                                                return true;
                                        }
                                        improved = local_search(local_best,new_sol, depth + 1, task, improved);
                                        if (improved == true) return true;
                                }
                        }
                }
        }
        return improved;
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
int random_selection(const std::vector<int>& int_vec, std::optional<unsigned int> seed = std::nullopt) {
        std::mt19937 gen;

        if (seed.has_value()) {
                gen.seed(seed.value());
        } else {
                std::random_device rd;
                gen.seed(rd());
        }

        std::uniform_int_distribution<> dist(0, int_vec.size() - 1);
        return int_vec[dist(gen)];
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
        std::vector<int> station_range = range_excluding(earliest, latest, station_idx);
        if (station_range.empty()) {
                return -1;
        }
        int selected_station = station_range[0];
        for (int i = 1; i < station_range.size(); ++i) {
                if (sol.load[station_range[i]] < sol.load[selected_station]) {
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
void VDLS::perturbation(ALBPSolution &new_sol) {
        std:: vector<int> can_change;
        can_change.resize(new_sol.n_stations, 1 );//Keeps track of stations that were already perturbed
        for (int i=0; i<n_perts_; i++) {
                std::vector<int> max_indices = get_max_indices(new_sol.load, can_change);
                //if there are no more max indices, that means that we have ran out of changeable stations
                if (max_indices.empty()) {
                        return;
                }
                const int station_idx = random_selection(max_indices);
                std::vector<int> filtered_tasks = check_windows(new_sol, station_idx);
                if (filtered_tasks.empty()) { //Stations tasks are all forced into that station, try different station
                        can_change[station_idx] = 0;
                        continue;
                }
                const int task_idx = random_selection(filtered_tasks);
                const int task = new_sol.station_assignments[station_idx][task_idx];
                const int new_station = select_new_station(new_sol, task, station_idx);
                perform_shift(new_sol, task, task_idx, station_idx, new_station);
                can_change[new_station] = 0;




        }

}


ALBPSolution VDLS::vdls_heuristic( int n_stations,  int lb) {
        ALBPSolution local_best = hoff_search(n_stations);
        ALBPSolution new_sol = local_best;
        while ((n_attempts_ < max_attempts_ ) && (local_best.cycle_time > lb) &&(!time_exceeded())) {
                bool improved = false;
                do {//Perform local search, restarting if we find an improved solution with depth 0
                        improved = local_search(local_best,new_sol, 0, -1, improved);
                        if (improved == true) new_sol = local_best;
                } while (improved == true);
                perturbation(new_sol);
                n_attempts_ ++;
        }
        return local_best;
}


ALBPSolution vdls_solve_salbp1(const ALBP &albp, std::optional<int> max_attempts , std::optional<int> time_limit   ) {
        int attempts = max_attempts.value_or(5000);  // default if not passed
        int limit = time_limit.value_or(7200);

        auto vdls= VDLS(albp, attempts, limit);
        ALBPSolution result =vdls.solve_type_1();
        return result;
}

ALBPSolution vdls_solve_salbp2(const ALBP &albp, std::optional<int> max_attempts , std::optional<int> time_limit   ) {
        int attempts = max_attempts.value_or(5000);  // default if not passed
        int limit = time_limit.value_or(7200);

        auto vdls= VDLS(albp, attempts, limit);
        ALBPSolution result =vdls.solve_type_2();
        return result;
}
ALBPSolution vdls_solve_salbp2(const int S,const int N, const std::vector<int>& task_times, const std::vector<std::vector<int>>& raw_precedence, const std::vector<int> &initial_solution,std::optional<int> max_attempts,
                                std::optional<int> time_limit) {
        int attempts = max_attempts.value_or(5000);  // default if not passed
        int limit = time_limit.value_or(7200);      // default if not passed
        ALBP albp = ALBP::type_2(S, N, task_times, raw_precedence);
        auto vdls= VDLS(albp, attempts, limit);
        if (!initial_solution.empty()) {
                vdls.add_init_solution(initial_solution);
        }
        ALBPSolution best = vdls.solve_type_2();
        return best;
}

ALBPSolution vdls_solve_salbp1(const int C,const int N, const std::vector<int>& task_times, const std::vector<std::vector<int>>& raw_precedence, const std::vector<int> &initial_solution,std::optional<int> max_attempts,
                                std::optional<int> time_limit) {
        int attempts = max_attempts.value_or(5000);  // default if not passed
        int limit = time_limit.value_or(7200);      // default if not passed
        ALBP albp = ALBP::type_1(C, N, task_times, raw_precedence);
        auto vdls= VDLS(albp, attempts, limit);
        if (!initial_solution.empty()) {
                vdls.add_init_solution(initial_solution);
        }
        ALBPSolution best = vdls.solve_type_1();
        return best;
}
