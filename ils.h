#ifndef ILS_H
#define ILS_H
#include "ALBP.h"
#include <string>
#include <vector>


struct ALBPSolution{
    std::vector<int> task_assignment; //Solution
    std::vector<std::vector<int>> station_assignments; //Solution
    int n_stations; 
    int n_tasks; //Number of tasks
    int n_violations; //Number of violations
    // int num_cycles; //Number of cycles
    // std::vector<int> cycle_times; //Cycle times
    

    void print();
    void task_to_station();
    void station_to_task();
};

int calc_lb_1(const std::vector<int>& task_time, int C);

ALBPSolution generate_approx_solution(ALBP&albp);
ALBPSolution iterated_local_search(ALBP& albp, int max_iter, float op_probs);

#endif // ILS_H