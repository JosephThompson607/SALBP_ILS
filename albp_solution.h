#ifndef ALBP_SOLUTION_H
#define ALBP_SOLUTION_H
#include "ALBP.h"
#include <string>
#include <vector>


struct ALBPSolution{
private:
    int n_tasks; //Number of tasks
public:
    std::vector<int> task_assignment; //Solution
    std::vector<std::vector<int>> station_assignments; //Solution
    std::vector<int> ranking; //Gives the tasks in order of ranking
    std::vector<int> task_ranking; //Gives the ranking for each task
    int n_stations{};

    int n_violations{}; //Number of violations
    int n_ranking_violations{}; // number of violations from ranking
    // int num_cycles; //Number of cycles
    // std::vector<int> cycle_times; //Cycle times

    //constructor
    explicit ALBPSolution(const int n_tasks) : n_tasks(n_tasks) {}

    [[nodiscard]] int get_n_tasks() const { return n_tasks; }  // Read-only access
    //utility functions
    void print() const;
    void task_to_station();
    void station_to_task();
    void station_to_ranking();
    void ranking_to_task_ranking();
};

int calc_lb_1(const std::vector<int>& task_time, int C);

ALBPSolution generate_approx_solution( const ALBP&albp,  int n_random);
void shallow_task_assignment( const ALBP&albp,  ALBPSolution& solution);
void task_oriented_assignment( const ALBP& albp,ALBPSolution& solution);
int count_violations(const ALBP&albp, const std::vector<int>& task_assignment);
#endif