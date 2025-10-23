
#include "albp_solution.h"
#include "ALBP.h"
#include "Hoff.h"
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
#include <functional>
#include <stdexcept>  // For std::runtime_error

 void ALBPSolution::print() const {
    std::cout << "ALBP Solution with S: " << n_stations << "  C: " <<cycle_time << "found with: " << method<<std::endl;
    for (int i = 0; i < n_stations; ++i) {
        std::cout << "Station " << i + 1 <<" load "<< loads[i] << " assigned tasks : ";
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
    for (int i = 0; i < n_tasks; ++i) {

        station_assignments[task_assignment[i]].push_back(i);
    }
}

void ALBPSolution::station_to_load(const ALBP &albp) {
    loads.clear();
    loads.resize(n_stations);
     int max_load = 0;
    for (int i = 0; i < n_stations; ++i) {

        for ( const int task: station_assignments[i]) {
            loads[i] += (albp.task_time[task]);
        }
        if (loads[i] > max_load) {
            max_load = loads[i];
        }
    }
     cycle_time = max_load;
}
/* Finds earliest and latest stations for all task **/
void ALBPSolution::find_windows(const ALBP &albp) {
    find_all_earliest(albp);
    find_all_latest(albp);
}

/* Finds earliest stations for a task **/
void ALBPSolution::find_earliest(const ALBP &albp, int i) {
    for (int j: albp.dir_pred[i]){
        if (task_assignment[j] > earliest[i]){
            earliest[i] = task_assignment[j];
        }
    }
}

void ALBPSolution::find_all_earliest(const ALBP &albp) {
    earliest.resize(albp.N,0);
    for (int i =0; i < albp.N; i++){
        find_earliest(albp, i);

    }
}

/* Finds latest stations for a task **/
void ALBPSolution::find_latest(const ALBP &albp, const int i) {
    for (int j: albp.dir_suc[i]) {
        if (task_assignment[j] < latest[i]){
            latest[i] = task_assignment[j];
        }
    }
}
/* updates the latest station of the predecessors of a task, in the case where the task is before all other successors*/
void ALBPSolution::update_pred_latest(const ALBP &albp, const int i) {
    for (int j: albp.dir_pred[i]) {
        if (task_assignment[i] < latest[j]) {
            latest[j] = task_assignment[i];
        }
    }
}

void ALBPSolution::update_suc_earliest(const ALBP &albp, const int i) {
    for (int j: albp.dir_suc[i]) {
        if (task_assignment[i] > earliest[j]) {
            earliest[j] = task_assignment[i];
        }
    }
}

void ALBPSolution::update_window(const ALBP &albp, const int i) {
    update_pred_latest(albp,i);
    update_suc_earliest(albp,i);
}
 void ALBPSolution::find_all_latest(const ALBP &albp){
    latest.resize(albp.N,n_stations-1);
    for (int i =0; i < albp.N; i++){
        find_latest(albp, i);

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

void ALBPSolution::station_to_ranking() {
    // Convert task assignment to ranking
    ranking.clear();
    task_ranking.resize(n_tasks);
    int ranking_counter = 0;
    for (int i = 0; i < n_stations; ++i) {
        for (int j : station_assignments[i]) {
            task_ranking[j] = ranking_counter;
            ranking.push_back(j);
            ++ranking_counter;
        }
    }
}

void ALBPSolution::ranking_to_task_ranking() {
    //Converts list of tasks ordered by rank to vector of ranks for each task ordered sequentially
    task_ranking.clear();
    task_ranking.resize(n_tasks, 0);
    for (int i = 0; i < n_tasks; ++i) {
        task_ranking[ranking[i]] = i;
    }
}

