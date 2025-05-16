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

std::tuple<int, int> get_task_bracket(const std::vector<int> &ranking, const ALBP &albp, const int task_index) {
    int min_suc = albp.N+1;
    int max_pred= 0;
    //goes through the precedence
    // for (int i = index1; i < albp.N; ++i) {
    //     for (int j = index1; j < albp.N; ++j) {
    //         if (albp.prec_mat[i*albp.N + j] && ranking[j] < min_suc) {
    //             min_suc = ranking[j];
    //         }
    //         if (albp.prec_mat[j*albp.N + i] && ranking[i] > max_pred) {
    //             max_pred = ranking[i];
    //         }
    //     }
    // }
    for (const int pred : albp.dir_suc[task_index]) {
        if ( ranking[pred] < min_suc) {
            min_suc = ranking[pred];
        }
    }
    for (const int pred : albp.dir_pred[task_index]) {
        if ( ranking[pred] > max_pred) {
            max_pred = ranking[pred];
        }
    }

    return std::make_tuple(min_suc, max_pred);
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
    auto [max_pred, min_suc] = get_task_bracket(ranking, albp, index1);
    //Note: allowing the float shift op to leave the thing in place if range is restricted
    std::uniform_int_distribution<> source_dist(max_pred, min_suc);
    insertion( ranking,  index1, source_dist(gen));

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
            if (op_probs < distribution(generator)) {
                exchange_op(new_solution.ranking, albp);

            }
            else {
                insertion_op(new_solution.ranking, albp);
            }
        }
        shallow_task_assignment(albp, new_solution);
        if (new_solution.n_violations < solution.n_violations && new_solution.n_stations <= solution.n_stations) {
            task_oriented_assignment( albp, new_solution);
            solution = new_solution;
            ni = 0;
        }
        else {
            ni ++;
        }
    }
}


ALBPSolution iterated_local_search(const ALBP &albp, const int max_iter, float op_probs) {


    // Initialize an initial (potentially infeasible) solution
    ALBPSolution best_solution = generate_approx_solution(albp);
    //prints the best solution
    std::cout << "Best task assignment: ";
    for (int i = 0; i < best_solution.task_assignment.size(); ++i) {
        std::cout << best_solution.task_assignment[i] + 1 << " ";
    }
    ALBPSolution candidate = best_solution;
    //Improves the solution with local search
    int iter = 0;
    const int lb_1 = calc_lb_1(albp.task_time, albp.C);
    local_search(candidate, albp, op_probs);
    //improves solution until iteration are reached or lower bound is reached
    while (iter < max_iter && candidate.n_stations > lb_1) {
        // Perform local search
        inversion_op(candidate.ranking, albp);
        shallow_task_assignment(albp, candidate);
        local_search(candidate, albp, op_probs);
        if (candidate.n_stations < best_solution.n_stations && candidate.n_violations  == 0) {

            best_solution = candidate;
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