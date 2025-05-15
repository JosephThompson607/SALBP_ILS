
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
#include <climits>
#include <random>
#include <ctime>
#include <vector>

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



void exchange_op(std::vector<int>& ranking, const ALBP& _) {
    /*changes the ranking of an ALBPSolution by swapping two elements in list*/
    // Set up random number generation
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<> dist(0, static_cast<int>(ranking.size()) - 1);

    // Pick two distinct random indices
    const int index1 = dist(gen);
    int index2 = dist(gen);
    while (index1 == index2) {
        index2 = dist(gen);
    }



    // Swap the selected elements

    std::swap(ranking[index1], ranking[index2]);

}


void insertion(std::vector<int> &ranking, const int from, int to) {
    if(from < to ) {
        std::rotate(ranking.begin()+from, ranking.begin()+from+1,ranking.begin()+ to+1);
    } else if (from > to) {
        std::rotate(ranking.begin()+ to, ranking.begin()+from, ranking.begin()+from+1);
    }
}


void float_shift_op(std::vector<int>& ranking, const ALBP& albp) {
    /*changes the ranking of an ALBPSolution using float shift operator*/
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<> dist(0, static_cast<int>(ranking.size()) - 1);
    const int index1 = dist(gen);
    //Gets the latest predecessor and earliest successor
    int min_success = albp.N+1;
    int max_pred= 0;
    //goes through the precedence
    for (int i = index1; i < albp.N; ++i) {
        for (int j = index1; j < albp.N; ++j) {
            if (albp.prec_mat[i*albp.N + j] && ranking[j] < min_success) {
                min_success = ranking[j];
            }
            if (albp.prec_mat[j*albp.N + i] && ranking[i] > max_pred) {
                max_pred = ranking[i];
            }
        }
    }
    //Note: allowing the float shift op to leave the thing in place if range is restricted
    std::uniform_int_distribution<> source_dist(max_pred, min_success);
    insertion( ranking,  index1, source_dist(gen));

}



void insertion_op(std::vector<int>& ranking, const ALBP& albp, const int range_start =0, int range_end = 0) {
    /*changes the ranking of an ALBPSolution by reinserting it into list*/
    if (range_end==0) {
       range_end = static_cast<int>(ranking.size()) - 1;
    }
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<> dist(range_start, range_end);

    // Pick two distinct random indices
    const int from = dist(gen);
    int to = dist(gen);
    while (from == to) {
        to = dist(gen);
    }
    insertion(ranking, from, to);
}

void deep_task_assignment(ALBPSolution& solution, const ALBP& albp) {
    std::vector<int> task_assignment(albp.N);
    std::vector<int> ranking = solution.ranking;

}
void local_search(ALBPSolution& solution,const ALBP& albp,const float op_probs , const int n_tries=10) {
    // Implement the local search algorithm here
    // This function should modify the solution in place
    // For example, you can swap tasks between stations or move tasks around
    // to improve the solution
    std::cout << "Performing local search..." << std::endl;
    std::random_device rd;  // Uses hardware randomness if available
    std::mt19937 generator(rd());
    std::uniform_real_distribution<float> distribution(0.0f, 1.0f);
    int ni = 0;
    while (ni  <= n_tries){
        ALBPSolution new_solution =solution;
        if (solution.n_violations == 0) {
            float_shift_op(new_solution.ranking, albp);
        }
        else {
            std::cout << std::endl << " Before swap " <<std::endl;
            for (int num : new_solution.ranking) {
                std::cout << num << " " ;
            }
            if (op_probs < distribution(generator)) {
                exchange_op(new_solution.ranking, albp);

            }
            else {
                insertion_op(new_solution.ranking, albp);
            }
            std::cout <<std::endl << " after swap " <<std::endl;
            for (int num : new_solution.ranking) {
                std::cout << num << " " ;
            }
        }
        shallow_task_assignment(albp, new_solution);
        if (new_solution.n_violations < solution.n_violations && new_solution.n_stations <= solution.n_stations) {
            deep_task_assignment(new_solution, albp);
            solution = new_solution;
            ni = 0;
        }
        else {
            ni ++;
        }
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
        //caculates the number of violations
        solution.n_violations = count_violations(albp, solution.task_assignment);
            
    }
}

std::vector<std::vector<int> > generate_rankings(const ALBP &albp, const int n_random) {
    //Generating task number ranking
    std::vector<std::vector<int>> rankings;

    // GEnerate PW ranking

    // std::vector<int> pw_rank = pw_ranking( albp);
    // rankings.push_back(pw_rank);
    // // Generate more rankings
    // for (int i = 0; i < n_random; ++i) {
    //     std::vector<int> ranking = random_ranking( albp);
    //     rankings.push_back(ranking);
    // }
    // std::vector<int> task_number_rank = task_number_ranking(albp);
    // rankings.push_back(task_number_rank);
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

ALBPSolution iterated_local_search(const ALBP &albp, const int max_iter, float op_probs) {

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