
#include "albp_solution.h"
#include "ALBP.h"
#include "mhh.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <filesystem>
#include <numeric>
#include <algorithm>
#include <climits>
#include <random>
#include <ctime>
#include <vector>
#include <tuple>
#include <functional>
#include <stdexcept>  // For std::runtime_error

void ALBPSolution::print() const {
    std::cout << "ALBP Solution with " << n_stations << " stations: " <<std::endl;
    for (int i = 0; i < n_stations; ++i) {
        std::cout << "Station " << i + 1 << ": ";
        for (int j : station_assignments[i]) {
            std::cout << j + 1 << " ";
        }
        std::cout << std::endl;
    }
    std::cout << "Here is the ranking :" << std::endl;
    for (int i = 0; i < n_tasks; ++i) {
        std::cout << ranking[i] << " ";
    }
    std::cout << std::endl;
}
void ALBPSolution::task_to_station(){
    // Convert task assignment to station assignment
    station_assignments.clear();
    station_assignments.resize(n_stations);
    for (int i = 0; i < n_tasks; ++i) {

        station_assignments[task_assignment[i]].push_back(i);
    }
}

void ALBPSolution::station_to_load(const ALBP &albp) {
    load.clear();
    load.resize(n_stations);
    for (int i = 0; i < n_stations; ++i) {

        for ( const int task: station_assignments[i]) {
            load[i] += (albp.task_time[task]);
        }
    }
}
/* Finds earliest and latest stations for all task **/
void ALBPSolution::find_windows(const ALBP &albp) {
    find_all_earliest(albp);
    find_all_latest(albp);
}

/* Finds earliest stations for a task **/
void ALBPSolution::find_earliest(const ALBP &albp, int i) {
    for (int j: albp.dir_pred[i]){
        if (task_assignment[j] > earliest[i]){
            earliest[i] = task_assignment[j];
        }
    }
}

void ALBPSolution::find_all_earliest(const ALBP &albp) {
    earliest.resize(albp.N,0);
    for (int i =0; i < albp.N; i++){
        find_earliest(albp, i);

    }
}

/* Finds latest stations for a task **/
void ALBPSolution::find_latest(const ALBP &albp, const int i) {
    for (int j: albp.dir_suc[i]) {
        if (task_assignment[j] < latest[i]){
            latest[i] = task_assignment[j];
        }
    }
}
/* updates the latest station of the predecessors of a task, in the case where the task is before all other successors*/
void ALBPSolution::update_pred_latest(const ALBP &albp, const int i) {
    for (int j: albp.dir_pred[i]) {
        if (task_assignment[i] < latest[j]) {
            latest[j] = task_assignment[i];
        }
    }
}

void ALBPSolution::update_suc_earliest(const ALBP &albp, const int i) {
    for (int j: albp.dir_suc[i]) {
        if (task_assignment[i] > earliest[j]) {
            earliest[j] = task_assignment[i];
        }
    }
}

void ALBPSolution::update_window(const ALBP &albp, const int i) {
    update_pred_latest(albp,i);
    update_suc_earliest(albp,i);
}
 void ALBPSolution::find_all_latest(const ALBP &albp){
    latest.resize(albp.N,n_stations-1);
    for (int i =0; i < albp.N; i++){
        find_latest(albp, i);

        }
    }



void ALBPSolution::station_to_task(){
    // Convert station assignment to task assignment
    task_assignment.clear();
    task_assignment.resize(n_tasks, -1);
    for (int i = 0; i < n_stations; ++i) {
        for (int j = 0; j < station_assignments[i].size(); ++j) {
            int task = station_assignments[i][j];
            if (task >= 0 && task < n_tasks) {
                task_assignment[task] = i;
            }
        }
    }
}

void ALBPSolution::station_to_ranking() {
    // Convert task assignment to ranking
    ranking.clear();
    task_ranking.resize(n_tasks);
    int ranking_counter = 0;
    for (int i = 0; i < n_stations; ++i) {
        for (int j : station_assignments[i]) {
            task_ranking[j] = ranking_counter;
            ranking.push_back(j);
            ++ranking_counter;
        }
    }
}

void ALBPSolution::ranking_to_task_ranking() {
    //Converts list of tasks ordered by rank to vector of ranks for each task ordered sequentially
    task_ranking.clear();
    task_ranking.resize(n_tasks, 0);
    for (int i = 0; i < n_tasks; ++i) {
        task_ranking[ranking[i]] = i;
    }
}


int calc_lb_1(const std::vector<int>& task_time, const int C) {
    int lb_1 = std::accumulate(task_time.begin(), task_time.end(), 0);
    lb_1 = (lb_1 + C - 1) / C; // ceil(lb_1 / C)
    return lb_1;
}

int calc_lb_2(const std::vector<int>& task_time, const int C) {
    int red_count = 0; //tasks over C/2
    int blue_count = 0; // tasks under or equal to C/2
    for (const auto& task : task_time) {
        if (static_cast<double>(task) > static_cast<double>(C)/2) {
            ++red_count;
        }
        else {
            ++blue_count;
        }
    }
    return red_count + (blue_count+1)/2;
}

int calc_salbp_2_lb_1(const std::vector<int>& task_time, const int S) {
    return (std::accumulate(task_time.begin(), task_time.end(), 0)+ S-1)/ S;
}

int calc_salbp_2_lb_2(const std::vector<int>& task_time, const int S) {
    return ( *std::max_element(task_time.begin(),task_time.end()));
}

int calc_salb_2_lb_3(const std::vector<int>& task_time, const int S) {
    std::vector<int> sorted = task_time;
    std::sort(sorted.begin(), sorted.end());
    int hole_size = size(task_time)/S;
    return std::accumulate(sorted.begin(), sorted.begin()+hole_size, 0) ;
}
int calc_salbp_2_lbs(const std::vector<int>& task_time, const int S) {
    std::vector<int> lbs = {calc_salbp_2_lb_1(task_time,S), calc_salbp_2_lb_2(task_time,S), calc_salb_2_lb_3(task_time,S)};
    return ( *std::max_element(lbs.begin(),lbs.end()));
}

int calc_salbp_2_ub(const std::vector<int>& task_time, int S) {
    int length = size(task_time);
    int up = 0;
    if (length ==1) return task_time[0]; //Only 1 element, C is its duration
    const int t_max = *std::max_element(task_time.begin(),task_time.end());
    const int t_sum = std::accumulate(task_time.begin(), task_time.end(), 0);
    if (length % 2 == 0) {
        up = (2 * t_sum + S-1)/ S;
    }
    else {
        up = (2 * t_sum + S-1) / (S +1);
    }
    return std::max(t_max, up);
}

void task_oriented_assignment(const ALBP& albp,ALBPSolution& solution) {
    std::vector<int> task_assignment(albp.N, -1);
    std::vector<int> ranking = solution.ranking;
    //std::vector<int> station_time(albp.N, 0);//Using N upper bound for station times
    std::vector station_cap(albp.N, albp.C);
    int task_index = 0;
    int n_stations = 0;
    while (!ranking.empty()) {
        const int current_task = ranking[task_index];
        bool can_assign = true;
        int left_station = 0;
        for (const int pred : albp.dir_pred[current_task]) {//Check if predecessors are all assigned
            if (task_assignment[pred] == -1) {
                can_assign = false;
                break;
            }
            if (task_assignment[pred] > left_station) {
                left_station = task_assignment[pred];
            }
        }
        if (can_assign) { //assign task to first available station
            for (int i = left_station; i < albp.N; ++i) {
                if (station_cap[i] - albp.task_time[current_task] >=0) {
                    task_assignment[current_task] = i;
                    station_cap[i] -= albp.task_time[current_task];
                    ranking.erase(ranking.begin()+task_index); //Potentially slow, will see
                    if (i + 1> n_stations) {
                        n_stations = i+1;
                    }
                    break;
                }
            }
        task_index = 0; //go back to check first task
        }
        else {
            task_index += 1;
        }

    }
    solution.n_stations = n_stations;
    solution.task_assignment = task_assignment;
    solution.task_to_station();
    solution.n_violations = count_violations(albp, solution.task_assignment);//real assignment violations
    solution.n_ranking_violations = count_violations(albp, solution.task_ranking); //ranking violations


    if (solution.n_violations > 0) {
        throw std::runtime_error("Deep task assignment failed to create feasible solution");
    }


}



std::vector<int> random_ranking(const ALBP&albp) {

    std::vector<int> ranking(albp.N);
    // Fill with 0 to N-1
    std::iota( ranking.begin(), ranking.end(), 0 );

    // Shuffle
    std::random_device rd;
    std::default_random_engine rng(rd());
    std::shuffle(ranking.begin(), ranking.end(), rng);
    return ranking;
}
std::vector<int> ascending_task_time_ranking(const ALBP&albp) {
    std::vector<int> ranking(albp.N);
    std::vector<int> indices(albp.N);
    std::iota( indices.begin(), indices.end(), 0 );
    // Sort indices based on the values in taskTimes
    std::sort(indices.begin(), indices.end(), [&](int i, int j) {
        return albp.task_time[i] < albp.task_time[j];
    });

    // Assign ranks

    for (int rank = 0; rank < albp.N; ++rank) {
        ranking[indices[rank]] = rank;
    }

    return ranking;
}
std::vector<int> descending_task_time_ranking(const ALBP&albp) {
    std::vector<int> ranking(albp.N);
    std::vector<int> indices(albp.N);
    std::iota( indices.begin(), indices.end(), 0 );
    // Sort indices based on the values in taskTimes
    std::sort(indices.begin(), indices.end(), [&](int i, int j) {
        return albp.task_time[i] > albp.task_time[j];
    });

    // Assign ranks

    for (int rank = 0; rank < albp.N; ++rank) {
        ranking[indices[rank]] = rank;
    }

    return ranking;
}


std::vector<int> inverse_task_number_ranking(const ALBP&albp) {
    std::vector<int> ranking(albp.N);
    // Fill with 0 to N-1
    std::iota( ranking.begin(), ranking.end(), 0 );
    std::reverse(ranking.begin(), ranking.end());
    return ranking;
}
std::vector<int> task_number_ranking(const ALBP&albp) {
    std::vector<int> ranking(albp.N);
    // Fill with 0 to N-1
    std::iota( ranking.begin(), ranking.end(), 0 );
    return ranking;
}


std::vector<int> n_suc_ranking(const ALBP&albp) {
    std::vector<int> ranking(albp.N);
    std::vector<int> indices(albp.N);
    std::iota( indices.begin(), indices.end(), 0 );
    // Sort indices based on the values in taskTimes
    std::sort(indices.begin(), indices.end(), [&](int i, int j) {
        return albp.suc[i].size() > albp.suc[j].size();
    });
    for (int rank = 0; rank < albp.N; ++rank) {
        ranking[indices[rank]] = rank;
    }

    return ranking;
}

std::vector<int> n_parents_ranking(const ALBP&albp) {
    std::vector<int> ranking(albp.N);
    std::vector<int> indices(albp.N);
    std::iota( indices.begin(), indices.end(), 0 );
    // Sort indices based on the values in taskTimes
    std::sort(indices.begin(), indices.end(), [&](int i, int j) {
        return albp.dir_pred[i].size() > albp.dir_pred[j].size();
    });
    for (int rank = 0; rank < albp.N; ++rank) {
        ranking[indices[rank]] = rank;
    }

    return ranking;
}

std::vector<int> n_children_ranking(const ALBP&albp) {
    std::vector<int> ranking(albp.N);
    std::vector<int> indices(albp.N);
    std::iota( indices.begin(), indices.end(), 0 );
    // Sort indices based on the values in taskTimes
    std::sort(indices.begin(), indices.end(), [&](int i, int j) {
        return albp.dir_suc[i].size() > albp.dir_suc[j].size();
    });
    for (int rank = 0; rank < albp.N; ++rank) {
        ranking[indices[rank]] = rank;
    }

    return ranking;
}

std::vector<int>  get_positional_weight(const ALBP &albp) {
    std::vector<int> weights(albp.N);
    for (int i = 0; i < albp.N; ++i) {
        weights[i] += albp.task_time[i]; // Own weight
        //Iterates through the transitive closure matrix and adds to weights
        for (const auto &suc : albp.suc[i]) {
            weights[i] += albp.task_time[suc];
        }
    }
    return weights;
}


std::vector<int>  get_alloc_station_ub(const ALBP &albp) {
    std::vector<int> weights(albp.N);
    for (int i = 0; i < albp.N; ++i) {
        int pw = albp.task_time[i]; // Own weight
        //Iterates through the transitive closure matrix and adds to weights
        for (const auto &suc : albp.suc[i]) {
            pw += albp.task_time[suc];
        }
        weights[i]  = albp.N + 1 - (albp.C + pw-1)/ albp.C; //using integer trick for rounding
    }
    return weights;
}

std::vector<int>  get_alloc_station_lb(const ALBP &albp) {
    std::vector<int> weights(albp.N);
    for (int i = 0; i < albp.N; ++i) {
        int pw = albp.task_time[i]; // Own weight
        //Iterates through the transitive closure matrix and adds to weights
        for (const auto &suc : albp.pred[i]) {
            pw += albp.task_time[suc];
        }
        weights[i]  = (albp.C + pw -1)/ albp.C; //using integer trick for rounding
    }
    return weights;
}


std::vector<int>  pw_ranking(const ALBP&albp) {
    std::vector<int> ranking(albp.N);

    //Gets the positional weight for each task
    std::vector<int> indices(albp.N);
    std::iota( indices.begin(), indices.end(), 0 );
    std::vector<int> weights = get_positional_weight(albp);

    //sorts

    // Sort indices based on corresponding weights
    std::sort(indices.begin(),indices.end(), [&](const int a, const int b) {
        return weights[a] > weights[b];  // Ascending order
    });
    for (int rank = 0; rank < indices.size(); ++rank) {
        ranking[indices[rank]] = rank;
    }
    return ranking;
}

std::vector<int> cumulative_pw_ranking(const ALBP&albp){
    //Gets the positional weight for each task
    std::vector<int> ranking(albp.N);
    std::vector<int> indices(albp.N);
    std::iota( indices.begin(), indices.end(), 0 );
    const std::vector<int> pw = get_positional_weight(albp);
    std::vector<int> cum_pw = albp.task_time;
    for (int task = 0; task < albp.N; ++task) {
        for (int const suc: albp.suc[task]) {
            cum_pw[task] += pw[suc];
        }
    }
    // Sort indices based on corresponding weights
    std::sort(indices.begin(),indices.end(), [&](const int a, const int b) {
        return cum_pw[a] > cum_pw[b];  // Ascending order
    });
    for (int rank = 0; rank < indices.size(); ++rank) {
        ranking[indices[rank]] = rank;
    }
    return ranking;
}

std::vector<int> average_pw_ranking(const ALBP&albp) {
    std::vector<int> ranking(albp.N);
    std::vector<int> indices(albp.N);
    std::iota( indices.begin(), indices.end(), 0 );
    std::vector<int> weights = get_positional_weight(albp);
    std::sort(indices.begin(),indices.end(), [&](const int a, const int b) {
    return weights[a]/ static_cast <double>(albp.suc[a].size()) > weights[b]/static_cast <double>(albp.suc[b].size());  // Ascending order
    });
    for (int rank = 0; rank < indices.size(); ++rank) {
        ranking[indices[rank]] = rank;
    }
    return ranking;
}

std::vector<int> as_up_avg_ranking(const ALBP&albp) {
    std::vector<int> ranking(albp.N);
    std::vector<int> indices(albp.N);
    std::iota( indices.begin(), indices.end(), 0 );
    std::vector<int> weights = get_alloc_station_ub(albp);
    std::sort(indices.begin(),indices.end(), [&](const int a, const int b) {
    return weights[a]/static_cast <double>( albp.suc[a].size() )< weights[b] /static_cast <double>(albp.suc[b].size());  // Ascending order
});
    for (int rank = 0; rank < indices.size(); ++rank) {
        ranking[indices[rank]] = rank;
    }
    return ranking;

}
std::vector<int> as_up_ranking(const ALBP&albp) {
    std::vector<int> ranking(albp.N);
    std::vector<int> indices(albp.N);
    std::iota( indices.begin(), indices.end(), 0 );
    std::vector<int> weights = get_alloc_station_ub(albp);
    std::sort(indices.begin(),indices.end(), [&](const int a, const int b) {
    return weights[a] < weights[b];  // Ascending order
});
    for (int rank = 0; rank < indices.size(); ++rank) {
        ranking[indices[rank]] = rank;
    }
    return ranking;
}

std::vector<int> t_over_as_up_ranking(const ALBP&albp) {
    std::vector<int> ranking(albp.N);
    std::vector<int> indices(albp.N);
    std::iota( indices.begin(), indices.end(), 0 );
    std::vector<int> weights = get_alloc_station_ub(albp);
    std::sort(indices.begin(),indices.end(), [&](const int a, const int b) {
    return static_cast<double>(albp.task_time[a]) / weights[a] < static_cast<double>(albp.task_time[b]) / weights[b];
});
    for (int rank = 0; rank < indices.size(); ++rank) {
        ranking[indices[rank]] = rank;
    }
    return ranking;
}
std::vector<int> as_lb_ranking(const ALBP&albp) {
    std::vector<int> ranking(albp.N);
    std::vector<int> indices(albp.N);
    std::iota( indices.begin(), indices.end(), 0 );
    std::vector<int> weights = get_alloc_station_lb(albp);


    std::sort(indices.begin(),indices.end(), [&](const int a, const int b) {
return  weights[a]  < weights[b];
    });
    for (int rank = 0; rank < indices.size(); ++rank) {
        ranking[indices[rank]] = rank;
    }
    return ranking;
}

std::vector<int> slack_ranking(const ALBP&albp) {
    std::vector<int> ranking(albp.N);
    std::vector<int> indices(albp.N);
    std::iota( indices.begin(), indices.end(), 0 );
    static std::vector<int> lb = get_alloc_station_lb(albp);
    static std::vector<int> ub = get_alloc_station_ub(albp);


    std::sort(indices.begin(),indices.end(), [&](const int a, const int b) {
return  ub[a]-lb[a]  < ub[b]-lb[b];
    });
    for (int rank = 0; rank < indices.size(); ++rank) {
        ranking[indices[rank]] = rank;
    }
    return ranking;
}

std::vector<int> suc_over_slack_ranking(const ALBP&albp) {
    std::vector<int> ranking(albp.N);
    std::vector<int> indices(albp.N);
    std::iota( indices.begin(), indices.end(), 0 );
    static std::vector<int> lb = get_alloc_station_lb(albp);
    static std::vector<int> ub = get_alloc_station_ub(albp);


    std::sort(indices.begin(),indices.end(), [&](const int a, const int b) {
return   static_cast<double>(albp.suc[a].size())/ (ub[a]-lb[a])  < static_cast<double>(albp.suc[b].size())/ (ub[b]-lb[b]);
    });
    for (int rank = 0; rank < indices.size(); ++rank) {
        ranking[indices[rank]] = rank;
    }
    return ranking;
}


int count_violations(const ALBP&albp, const std::vector<int>& task_assignment) {
    int violations = 0;
    for (const auto& rel : albp.precedence_relations) {
        int parent = rel.parent - 1; // Convert to 0-based index
        int child = rel.child - 1;   // Convert to 0-based index
        if (task_assignment[parent] > task_assignment[child]) {
            violations++;
        }
    }
    return violations;
}

void shallow_task_assignment( const ALBP&albp,  ALBPSolution& solution) {
    // Assign tasks to stations based on the ranking
    solution.n_stations = 0;
    solution.task_assignment.resize(albp.N, -1);
    solution.station_assignments.clear();
    solution.station_assignments.resize(albp.N); //using upper bound of N stations
    //creates a vector of size N with Cycle time
    std::vector<int> station_times(albp.N, albp.C);
    int max_station = 0;
    for (int i = 0; i < albp.N; ++i) {
        int task = solution.ranking[i];
        // Find the first station that can accommodate the task

        for (int j = 0; j <= albp.N; ++j) { //using upper bound of N stations
            if (station_times[j] >= albp.task_time[task]) {
                solution.task_assignment[task] = j;
                solution.station_assignments[j].push_back(task);
                station_times[j] -= albp.task_time[task];
                if (j > max_station) {
                    max_station = j;
                }
                break;
            }
        }
        solution.n_stations = max_station + 1;
        //calculates the number of violations
        solution.n_violations = count_violations(albp, solution.task_assignment);
        solution.n_ranking_violations = count_violations(albp, solution.task_ranking);
            
    }
}

std::vector<std::vector<int> > generate_rankings(const ALBP &albp, const int n_random) {

    //Generating task number ranking
    std::vector<std::function<std::vector<int>(const ALBP&)>> ranking_functions = {
         task_number_ranking,
         inverse_task_number_ranking,
         ascending_task_time_ranking,
         descending_task_time_ranking,
         n_parents_ranking,
         n_children_ranking,
        n_suc_ranking,
        pw_ranking,
        cumulative_pw_ranking,
        average_pw_ranking,
        as_up_ranking,
        as_lb_ranking,
        as_up_avg_ranking,
        t_over_as_up_ranking,
        slack_ranking,
        suc_over_slack_ranking,
    };
    std::vector<std::vector<int>> rankings;
    rankings.reserve(ranking_functions.size()+n_random);
    for (const auto& func : ranking_functions) {
        rankings.push_back(func(albp));
    }

    // // Generate more rankings
    for (int i = 0; i < n_random; ++i) {
        std::vector<int> ranking = random_ranking( albp);
        rankings.push_back(ranking);
    }

    return rankings;
}


std::vector<ALBPSolution> generate_solutions( const ALBP &albp, const int n_random) {
    // First generates a list of rankings for the tasks. It will be a vector of vectors, one for each ranking of size N
    // Then generates a solution for each ranking
    std::vector<ALBPSolution> solutions;
    std::cout << "generating initial solutions" << std::endl;
    for (const std::vector<std::vector<int>> rankings = generate_rankings(albp, n_random); const auto & ranking : rankings) {
        ALBPSolution solution(albp.N);
        solution.ranking = ranking;
        solution.ranking_to_task_ranking();
        shallow_task_assignment(albp, solution);
        solutions.push_back(solution);
    }
    //Add in hoffman solution
    solutions.push_back( mhh_solve_salbp1(albp));
    return solutions;
}

ALBPSolution process_init_solution( const ALBP &albp, const std::vector<int> &initial_solution) {
    /* Reads in an initial task assignment. Assumes indexing starts at zero and initial solution is feasible*/
    ALBPSolution solution(albp.N);
    solution.n_stations = *std::max_element(initial_solution.begin(), initial_solution.end()) + 1;
    solution.task_assignment = initial_solution;
    solution.task_to_station();
    solution.station_to_ranking();
    solution.ranking_to_task_ranking();
    solution.n_violations = 0;

    return solution;
}

ALBPSolution generate_approx_solution(const ALBP&albp, const int n_random, const std::vector<int> &initial_solution)
{
     std::vector<ALBPSolution>  candidates = generate_solutions(albp, n_random);
    //Add in optional initial solution to the pool of candidates
    if (!initial_solution.empty()) {
        candidates.push_back( process_init_solution(albp, initial_solution) );
    }


    // Select the best solution from the generated solutions
    std::cout <<"finding best initial solution" << std::endl;
    ALBPSolution best_solution = candidates[0];
    int n_stations = best_solution.n_stations;
    int n_violations = best_solution.n_violations;
    for (const auto& sol : candidates) {
        if (sol.n_stations < n_stations) { // selects solution with the least number of stations
            n_stations = sol.n_stations;
            n_violations = sol.n_violations;
            best_solution = sol;
        }
        else if (sol.n_stations == n_stations) { //tiebreaker-- smaller number of violations
            if (sol.n_violations < n_violations) {
                n_violations = sol.n_violations;
                best_solution = sol;
            }
        }
    }
    std::cout << "stations: " << n_stations << " violations: " << n_violations << std::endl;
    if (candidates.empty()) {
        throw std::runtime_error("No solutions generated");
    }

    return best_solution;
}
