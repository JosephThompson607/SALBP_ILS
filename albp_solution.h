#ifndef ALBP_SOLUTION_H
#define ALBP_SOLUTION_H
#include <chrono>
using std::chrono::milliseconds;
#include "ALBP.h"
#include <string>
#include <vector>


struct ALBPSolution{
private:
    int n_tasks; //Number of tasks
public:
    std::vector<int> task_assignment; //Solution
    std::vector<std::vector<int>> station_assignments; //Solution
    std::vector<int> loads;// how much time is used at each station (VDLS)
    std::vector<int> earliest; //Earliest station each task can be assigned to, based on predecessors/successors (VDLS)
    std::vector<int> latest; //Latest station each task can be assigned to, based on predecessors/successors (VDLS)
    milliseconds elapsed_ms; //Run time TODO implement for hoffman and ILS

    std::vector<int> ranking; //Gives the tasks in order of ranking (ILS)
    std::vector<int> task_ranking; //Gives the ranking for each task (ILS)
    int n_stations{};//For SALBP-1
    int cycle_time{};//For SALBP-2
    int n_violations{}; //Number of violations
    int n_ranking_violations{}; // number of violations from ranking
    int max_ovlo_;//Cycle time overloading
    int tot_ovlo2_;
    bool optimal = false;
    std::vector<int> overloads_;
    std::vector<int> sq_overloads_;// squared overloads at each station
    // int num_cycles; //Number of cycles
    // std::vector<int> cycle_times; //Cycle times
    std::string method = "unspecified";

    //constructor
    explicit ALBPSolution(const int n_tasks) : n_tasks(n_tasks), n_stations(n_tasks), task_assignment(n_tasks, -1), ranking(n_tasks, -1), task_ranking(n_tasks, -1), overloads_(0), sq_overloads_(0),max_ovlo_(0), tot_ovlo2_(0), elapsed_ms(0){
    }

    [[nodiscard]] int get_n_tasks() const { return n_tasks; }  // Read-only access
    //utility functions
    void print() const;
    void task_to_station();
    void task_to_station_and_load(const ALBP &albp);
    void station_to_task();
    void station_to_ranking();
    void station_to_load(const ALBP &albp);
    void find_all_earliest(const ALBP &albp);
    void find_all_latest(const ALBP &albp);
    void find_latest(const ALBP &albp, int i);
    void find_earliest(const ALBP &albp, int i);
    void update_pred_latest(const ALBP &albp, int i, int old_station);
    void update_suc_earliest(const ALBP &albp, int i, int old_station);
    void update_window(const ALBP &albp, int i, int old_station);
    void find_windows(const ALBP &albp);
    void reverse();
    void ranking_to_task_ranking();
    void task_ranking_to_ranking();

};


#endif