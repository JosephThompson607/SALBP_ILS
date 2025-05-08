
#include "ils.h"
#include "ALBP.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <filesystem>
#include <numeric>
#include <algorithm>
#include <random> 

void ALBPSolution::print() {
    std::cout << "ALBP Solution: " << std::endl;
    for (int i = 0; i < n_stations; ++i) {
        std::cout << "Station " << i + 1 << ": ";
        for (int j = 0; j < station_assignments[i].size(); ++j) {
            std::cout << station_assignments[i][j] + 1 << " ";
        }
        std::cout << std::endl;
    }
}
void ALBPSolution::task_to_station(){
    // Convert task assignment to station assignment
    station_assignments.clear();
    station_assignments.resize(n_stations);
    for (int i = 0; i < task_assignment.size(); ++i) {
        int station = task_assignment[i];
        if (station >= 0 && station < n_stations) {
            station_assignments[station].push_back(i);
        }
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


int calc_lb_1(const std::vector<int>& task_time, int C) {
    int lb_1 = std::accumulate(task_time.begin(), task_time.end(), 0);
    lb_1 = (lb_1 + C - 1) / C; // ceil(lb_1 / C)
    return lb_1;
}

void local_search(ALBPSolution& solution, ALBP& albp, float op_probs) {
    // Implement the local search algorithm here
    // This function should modify the solution in place
    // For example, you can swap tasks between stations or move tasks around
    // to improve the solution
    std::cout << "Performing local search..." << std::endl;
}

ALBPSolution iterated_local_search(ALBP& albp, int max_iter, float op_probs) {

    int iter = 0;
    int lb_1 = calc_lb_1(albp.task_time, albp.C);
    // Initialize an initial (poetentially infeasible) solution
    ALBPSolution best_solution = generate_approx_solution(albp);
    //prints the best solution
    std::cout << "Best task assignment: ";
    for (int i = 0; i < best_solution.task_assignment.size(); ++i) {
        std::cout << best_solution.task_assignment[i] + 1 << " ";
    }
    std::cout << std::endl;
    std::cout << "Best solution found: " << best_solution.n_stations << " stations" << std::endl;
    std::cout << "Best solution found: " << best_solution.n_violations << " violations" << std::endl;
    //Improves the solution with local search
    local_search(best_solution, albp, op_probs);
    //improves solution until iteration are reached or lower bound is reached
    while (iter < max_iter && best_solution.n_stations > lb_1) {
        // Perform local search
        local_search(best_solution, albp, op_probs);
        iter++;
    }
    // Return the best solution found
    return best_solution;
}

void random_ranking(std::vector<int>& ranking, int N) {
    std::vector<int> vec(N);

    // Fill with 0 to N-1
    for (int i = 0; i < N; ++i) {
        vec[i] = i;
    }

    // Shuffle
    std::random_device rd;
    std::default_random_engine rng(rd());
    std::shuffle(vec.begin(), vec.end(), rng);
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

void shallow_task_assignment(ALBP&albp, const std::vector<int>& ranking, ALBPSolution& solution) {
    // Assign tasks to stations based on the ranking
    solution.n_stations = 0;
    solution.task_assignment.resize(albp.N, -1);
    solution.station_assignments.clear();
    solution.station_assignments.resize(albp.N); //using upper bound of N stations
    //creates a vector of size N with Cycle time
    std::vector<int> station_times(albp.N, albp.C);
    for (int i = 0; i < albp.N; ++i) {
        int task = ranking[i];
        // Find the first station that can accommodate the task
        int max_station = 0;
        for (int j = 0; j < solution.n_stations; ++j) {
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
        //caculates the number of violations
        solution.n_violations = count_violations(albp, solution.task_assignment);
            
    }
}
void generate_solutions(ALBP&albp, std::vector<ALBPSolution>& solutions) {
    // First generates a list of rankings for the tasks. It will be a vector of vectors, one for each ranking of size N
    std::vector<std::vector<int>> rankings;
    // Generate more rankings
    for (int i = 1; i < 10; ++i) {
        std::vector<int> ranking(albp.N);
        random_ranking(ranking, albp.N);
        rankings.push_back(ranking);
    }
    // Then generates a solution for each ranking
    for (size_t i = 0; i < rankings.size(); ++i) {
        ALBPSolution solution;
        shallow_task_assignment(albp, rankings[i], solution);
        solutions.push_back(solution);
    }}

ALBPSolution generate_approx_solution(ALBP&albp)
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