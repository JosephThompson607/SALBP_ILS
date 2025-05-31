//
// Created by Joseph Thompson on 2025-05-16.
//

#include "ils.h"

#include <vector>
#include <numeric>
#include <algorithm>
#include <iostream>
#include <random>
#include <__ostream/basic_ostream.h>
#include "albp_solution.h"
#include "ALBP.h"

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

std::tuple<int, int> get_task_bracket(const ALBPSolution &sol, const ALBP &albp, const int task_index) {
    int min_suc = albp.N;
    int max_pred= 0;


    for (const int pred : albp.dir_suc[task_index]) {
        if ( sol.task_ranking[pred] < min_suc) {
            min_suc = sol.task_ranking[pred];
        }
    }
    for (const int pred : albp.dir_pred[task_index]) {
        if ( sol.task_ranking[pred] > max_pred) {
            max_pred = sol.task_ranking[pred];
        }
    }

    return std::make_tuple(max_pred, min_suc);
}

void insertion(std::vector<int> &ranking, const int from, const int to) {
    if (from == to) return;
    if(from < to ) {
        std::rotate(ranking.begin()+from, ranking.begin()+from+1,ranking.begin()+ to+1);
    } else if (from > to) {
        std::rotate(ranking.begin()+ to, ranking.begin()+from, ranking.begin()+from+1);
    }
}
void float_shift_op(ALBPSolution &sol, const ALBP& albp) {
    /*changes the ranking of an ALBPSolution using float shift operator*/
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<> dist(0, static_cast<int>(sol.ranking.size()) - 1);
    const int index1 = dist(gen);
    //Gets the latest predecessor and earliest successor
    auto [max_pred, min_suc] = get_task_bracket(sol, albp, index1);

    //Note: allowing the float shift op to leave the thing in place if range is restricted
    std::uniform_int_distribution<> source_dist(max_pred, min_suc);
    const int to = std::min(static_cast<int>(sol.ranking.size()) - 1, source_dist(gen));
    insertion( sol.ranking,  index1, to);
}


void inversion_op(std::vector<int>& ranking, const ALBP& albp) {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<> dist(0, albp.N - 1);

    // Pick two distinct random indices
    int from = dist(gen);
    int to = dist(gen);
    while (from == to) {
        to = dist(gen);
    }
    if (from > to) {
        int holder = from;
        from = to;
        to = holder;
    }
    std::reverse(ranking.begin()+from, ranking.begin()+to+1);

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

void local_search(ALBPSolution& solution,const ALBP& albp,const float op_probs , const int n_tries=10) {
    //solution.station_to_ranking();
    std::random_device rd;  // Uses hardware randomness if available
    std::mt19937 generator(rd());
    std::uniform_real_distribution<float> distribution(0.0f, 1.0f);
    int ni = 0;
    while (ni  <= n_tries){
        ALBPSolution new_solution =solution;
        if (new_solution.n_ranking_violations == 0) {
            float_shift_op(new_solution, albp);
        }
        else {
            if (op_probs < distribution(generator)) {
                exchange_op(new_solution.ranking, albp);

            }
            else {
                insertion_op(new_solution.ranking, albp);
            }
        }
        shallow_task_assignment(albp, new_solution);

        if (new_solution.n_violations <= solution.n_violations && new_solution.n_stations < solution.n_stations) {
            task_oriented_assignment( albp, new_solution);
            solution = new_solution;
            solution.station_to_ranking();
            ni = 0;
        }
        else {
            ni ++;
        }
    }
}


ALBPSolution iterated_local_search(const ALBP &albp, const int max_iter, float op_probs, bool verbose) {


    // Initialize an initial (potentially infeasible) solution

    ALBPSolution best_solution = generate_approx_solution(albp, 250);
    //prints the best solution
    ALBPSolution candidate = best_solution;
    //Improves the solution with local search
    int iter = 0;
    const int lb_1 = calc_lb_1(albp.task_time, albp.C);
    const int lb_2 = calc_lb_2(albp.task_time, albp.C);
    std::cout << "Performing local search. first iter " << std::endl;
    local_search(candidate, albp, op_probs, 50);
    // std::cout << "assigning tasks deep" << std::endl;
    // task_oriented_assignment(albp, candidate);

    //improves solution until iteration are reached or lower bound is reached
    // while (iter < max_iter && candidate.n_stations > lb_1) {
    while (iter < max_iter ) {
        if (verbose) {
            std::cout << "Performing local search. Iteration no: " << iter <<" n_violations: "<< candidate.n_violations << " n_stations " << candidate.n_stations << std::endl;
        }
        // Perform local search
        candidate.station_to_ranking();
        inversion_op(candidate.ranking, albp);
        shallow_task_assignment(albp, candidate);
        local_search(candidate, albp, op_probs);
        if ((candidate.n_stations < best_solution.n_stations && candidate.n_violations  <= best_solution.n_violations) || (candidate.n_stations <= best_solution.n_stations && candidate.n_violations < best_solution.n_violations)) {
            std::cout << "found good solution . Iteration no: " << iter << " n_violations: " <<  candidate.n_violations <<" n_stations: " << candidate.n_stations << std::endl;
            best_solution = candidate;
            if (candidate.n_stations == lb_1 && candidate.n_violations == 0) {
                std::cout << "Found optimal solution, proven by lb_1. terminating search at: " << iter << std::endl;
                break;
            }
            else if (candidate.n_stations == lb_2 && candidate.n_violations == 0) {
                std::cout << "Found optimal solution, proven by lb_1. terminating search at: " << iter << std::endl;
                break;
            }
        }
        iter++;
    }
    if (best_solution.n_violations > 0) {
        task_oriented_assignment(albp, best_solution );
        std::cout <<std::endl << "fixing broken solution" <<std::endl;

    }
    // Return the best solution found
    return best_solution;
}


ALBPSolution ils_solve_SALBP1(const int C,const int N, const std::vector<int> &task_times, const std::vector<std::vector<int> > &raw_precedence, int max_iter, float op_probs, bool verbose,  const std::vector<int> &initial_solution) {
    ALBP albp(C, N, task_times, raw_precedence);
    //TODO: Allow for passing initial solution
    ALBPSolution result =iterated_local_search(albp, max_iter, op_probs, verbose);
    return result;
}
