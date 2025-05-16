
#include "albp_solution.h"
#include "ALBP.h"
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
}
void ALBPSolution::task_to_station(){
    // Convert task assignment to station assignment
    station_assignments.clear();
    station_assignments.resize(n_stations);
    for (int i = 0; i < task_assignment.size(); ++i) {

        station_assignments[task_assignment[i]].push_back(i);
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
    solution.task_assignment = task_assignment;
    solution.task_to_station();
    solution.n_violations = count_violations(albp, solution.task_assignment);
    solution.n_tasks = albp.N;
    solution.n_stations = n_stations;
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

std::vector<int>  pw_ranking(const ALBP&albp) {
    std::vector<int> ranking(albp.N);
    std::vector<int> weights(albp.N);
    //Gets the positional weight for each task
    std::vector<int> indices(albp.N);
    std::iota( indices.begin(), indices.end(), 0 );

    for (int i = 0; i < albp.N; ++i) {
        weights[i] += albp.task_time[i]; // Own weight
        //Iterates through the transitive closure matrix and adds to weights
        for (int j = 0; j < albp.N; ++j) {
            weights[i]  += albp.t_close_mat[i * albp.N + j] * albp.task_time[j];
        }
    }

    //sorts

    // Sort indices based on corresponding weights
    std::ranges::sort(indices, [&](const int a, const int b) {
        return weights[a] > weights[b];  // Ascending order
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
    solution.n_tasks = albp.N;
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
            
    }
}

std::vector<std::vector<int> > generate_rankings(const ALBP &albp, const int n_random) {
    //Generating task number ranking
    std::vector<std::vector<int>> rankings;

    // GEnerate PW ranking

    std::vector<int> pw_rank = pw_ranking( albp);
    rankings.push_back(pw_rank);
    // Generate more rankings
    for (int i = 0; i < n_random; ++i) {
        std::vector<int> ranking = random_ranking( albp);
        rankings.push_back(ranking);
    }
    std::vector<int> task_number_rank = task_number_ranking(albp);
    rankings.push_back(task_number_rank);
    std::vector<int> inverse_task_number_rank = inverse_task_number_ranking(albp);
    rankings.push_back(inverse_task_number_rank);
    return rankings;
}

void generate_solutions(const ALBP&albp, std::vector<ALBPSolution>& solutions, const int n_random=3) {
    // First generates a list of rankings for the tasks. It will be a vector of vectors, one for each ranking of size N
    std::vector<std::vector<int>>  rankings = generate_rankings(albp, n_random);
    // Then generates a solution for each ranking
    for (const auto & ranking : rankings) {
        ALBPSolution solution;
        solution.ranking = ranking;
        shallow_task_assignment(albp, solution);
        solutions.push_back(solution);
    }
}

ALBPSolution generate_approx_solution(const ALBP&albp)
{
    std::vector<ALBPSolution> best_solutions;
    generate_solutions(albp, best_solutions);
    // Select the best solution from the generated solutions

    ALBPSolution best_solution;
    int n_stations = INT_MAX;
    int n_violations = INT_MAX;
    for (const auto& solution : best_solutions) {
        if (solution.n_stations < n_stations) { // selects solution with the least number of stations
            n_stations = solution.n_stations;
            n_violations = solution.n_violations;
            best_solution = solution;
        }
        else if (solution.n_stations == n_stations) { //tiebreaker-- smaller number of violations
            if (solution.n_violations < n_violations) {
                n_violations = solution.n_violations;
                best_solution = solution;
            }
        }
    }
    return best_solution;
}
