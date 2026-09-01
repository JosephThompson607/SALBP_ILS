
#include "albp_solution.h"
#include "ALBP.h"
#include "heuristics/Hoff.h"
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
#include <map>
#include <stdexcept>  // For std::runtime_error

 void ALBPSolution::print() const {
    std::cout << "ALBP Solution with S: " << n_stations << "  C: " <<cycle_time << " found with: " << method<<std::endl;
    for (int i = 0; i < n_stations; ++i) {
        std::cout << "Station " << i + 1 <<" load "<< loads[i] << " assigned tasks : ";
        for (int j : station_assignments[i]) {
            std::cout << j + 1 << " ";
        }
        std::cout << std::endl;
    }


}

void ALBPSolution::reverse()  {

     std::reverse(station_assignments.begin(), station_assignments.end());
     station_to_task();
     if (!loads.empty()) {
         std::reverse(loads.begin(), loads.end());
         auto max_it = std::max_element(loads.begin(),loads.end());
         max_station = std::distance(loads.begin(), max_it);
         cycle_time = *max_it;
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

void ALBPSolution::task_to_station_and_load(const ALBP &albp) {
     station_assignments.clear();
     station_assignments.resize(n_stations);
     loads.clear();
     loads.resize(n_stations);
     cycle_time= 0;
     for (int i = 0; i < n_tasks; ++i) {
         const int station = task_assignment[i];
         loads[station] += albp.task_time[i];
         station_assignments[station].push_back(i);
         if (loads[station] > cycle_time) {
             cycle_time = loads[station];
             max_station=station;
         }
     }
 }
void ALBPSolution::station_to_load(const ALBP &albp) {
    loads.clear();
    loads.resize(n_stations);
     cycle_time = 0;
    for (int i = 0; i < n_stations; ++i) {

        for ( const int task: station_assignments[i]) {
            loads[i] += (albp.task_time[task]);
        }
        if (loads[i] > cycle_time) {
            cycle_time = loads[i];
            max_station=i;
        }
    }
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
     earliest.clear();
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
void ALBPSolution::update_pred_latest(const ALBP &albp, const int i, int old_station) {
    for (int j: albp.dir_pred[i]) {
        if (task_assignment[i] < latest[j]) {
            latest[j] = task_assignment[i];
        }
        else if (task_assignment[i] > latest[j] && old_station == latest[j]) {
            //Search for latest station by taking the minimum of all successor stations
            int latest_st = task_assignment[albp.dir_suc[j][0]];
            for (int suc=1; suc<albp.dir_suc[j].size(); suc++) {
                if (task_assignment[albp.dir_suc[j][suc]] < latest_st) {
                    latest_st = task_assignment[albp.dir_suc[j][suc]];
                }
            }
            latest[j] = latest_st;
        }
    }
}

void ALBPSolution::update_suc_earliest(const ALBP &albp, const int i, int old_station) {
    for (int j: albp.dir_suc[i]) {
        if (task_assignment[i] > earliest[j]) {
            earliest[j] = task_assignment[i];
        }
        else if (task_assignment[i] < earliest[j] && old_station == earliest[j]) {
            //Search for earliest station by taking the maxiumum of all predecessor station assignments
            int earliest_station = task_assignment[albp.dir_pred[j][0]];
            for (int pred=1; pred<albp.dir_pred[j].size(); pred++) {
                if (task_assignment[albp.dir_pred[j][pred]] > earliest_station) {
                    earliest_station = task_assignment[albp.dir_pred[j][pred]];
                }
            }
            earliest[j] = earliest_station;
        }
    }
}

void ALBPSolution::update_window(const ALBP &albp, const int i, int old_station) {
    update_pred_latest(albp,i,old_station);
    update_suc_earliest(albp,i, old_station);
}
 void ALBPSolution::find_all_latest(const ALBP &albp){
     latest.clear();
     latest.resize(albp.N, n_stations - 1);
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

void ALBPSolution::station_to_ranking(bool sort_by_task) {
    // Convert task assignment to ranking
    ranking.clear();
    task_ranking.resize(n_tasks);
    int ranking_counter = 0;
    for (int i = 0; i < n_stations; ++i) {
        if (sort_by_task) {
            std::sort(station_assignments[i].begin(), station_assignments[i].end());
        }
        for (int j : station_assignments[i]) {
            task_ranking[j] = ranking_counter;
            ranking.push_back(j);
            ++ranking_counter;
        }
    }
}

void ALBPSolution::task_assignment_to_ranking() {
     ranking.clear();
     ranking.resize(n_tasks);

     std::iota(ranking.begin(), ranking.end(), 0);

     std::sort(ranking.begin(), ranking.end(),
         [&](int a, int b) {
             if (task_assignment[a] != task_assignment[b])
                 return task_assignment[a] < task_assignment[b];
             return a < b;
         });
     ranking_to_task_ranking();

 }

void ALBPSolution::ranking_to_task_ranking() {
    //Gets a vector with index of task and value of rank
    task_ranking.clear();
    task_ranking.resize(n_tasks, 0);
    for (int i = 0; i < n_tasks; ++i) {
        task_ranking[ranking[i]] = i;
    }
}
void ALBPSolution::task_ranking_to_ranking() {
     //Gets a vecor with the index of the rank and the value of the task
     ranking.clear();
     ranking.resize(n_tasks, 0);
     for (int i = 0; i < n_tasks; ++i) {
         ranking[task_ranking[i]] = i;
     }
 }

