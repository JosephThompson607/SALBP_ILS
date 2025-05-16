#ifndef ALBP_SOLUTION_H
#define ALBP_SOLUTION_H
#include "ALBP.h"
#include <string>
#include <vector>


struct ALBPSolution{
    std::vector<int> task_assignment; //Solution
    std::vector<std::vector<int>> station_assignments; //Solution
    std::vector<int> ranking;
    int n_stations; 
    int n_tasks; //Number of tasks
    int n_violations; //Number of violations
    // int num_cycles; //Number of cycles
    // std::vector<int> cycle_times; //Cycle times
    

    void print() const;
    void task_to_station();
    void station_to_task();
};

int calc_lb_1(const std::vector<int>& task_time, int C);

ALBPSolution generate_approx_solution(ALBP&albp);
void shallow_task_assignment( const ALBP&albp,  ALBPSolution& solution);
void task_oriented_assignment( const ALBP& albp,ALBPSolution& solution);
ALBPSolution generate_approx_solution(const ALBP& albp);
int count_violations(const ALBP&albp, const std::vector<int>& task_assignment);
#endif