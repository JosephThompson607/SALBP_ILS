//
// Created by Joseph Thompson on 2025-09-11.
//

#ifndef SALBP_BASICS_H
#define SALBP_BASICS_H


int calc_lb_1(const std::vector<int>& task_time, int C);
int calc_lb_2(const std::vector<int>& task_time, int C);
int calc_salbp_2_lbs(const std::vector<int>& task_time, int S);
int calc_salbp_2_ub(const std::vector<int>& task_time, int S);
std::vector<int>  get_positional_weight(const ALBP &albp);
ALBPSolution generate_approx_solution( const ALBP&albp,  int n_random,  const std::vector<int> &initial_solution = std::vector<int>());
ALBPSolution process_init_solution( const ALBP &albp, const std::vector<int> &initial_solution);
std::vector<ALBPSolution> generate_priority_ranking_solutions(const ALBP&albp, int n_random);
void shallow_task_assignment( const ALBP&albp,  ALBPSolution& solution);

void task_oriented_assignment( const ALBP& albp,ALBPSolution& solution);
ALBPSolution station_oriented_assignment( const ALBP& albp,ALBPSolution& solution);
int count_violations(const ALBP&albp, const std::vector<int>& task_assignment);
std::vector<ALBPSolution>  priority_solve_salbp_1( int C, int N, const std::vector<int>& task_times, const std::vector<std::vector<int>>& raw_precedence,  int n_random =100) ;
#endif //SALBP_BASICS_H
