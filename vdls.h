//
// Created by Joseph Thompson on 2025-07-02.
//

#ifndef VDLS_H
#define VDLS_H
#include "ALBP.h"
#include "albp_solution.h"


class VDLS {
public:
    //explicit ALBP(const std::string& filename) { loadFromFile(filename); }
    explicit VDLS(const ALBP& albp,  int max_attempts = 5000):albp_(albp), max_attempts_(max_attempts) {}
    ALBPSolution solve_type_1( );

private:
    //Variables
    const ALBP& albp_;
    ALBPSolution best_{albp_.N};
    int n_attempts_=0;
    int max_attempts_;
    int max_depth_=6; //number of shifts to consider
    int lb_;
    int n_perts_ = 6;
    std::vector<int>  pw_{};

    //Methods
    ALBPSolution vdls_heuristic(  int n_stations,   int lb);
    ALBPSolution hoff_search(int n_stations);
    bool local_search(ALBPSolution& local_best,ALBPSolution& prev_sol, int depth, int last_task=-1, bool improved = false);
    void perform_shift(ALBPSolution &sol, int task, int task_idx, int old_station, int new_station);
    void perturbation(ALBPSolution& incumbent_solution);
};

ALBPSolution vdls_solve_salbp1(const ALBP &albp);
ALBPSolution vdls_solve_salbp1( int C,int N, const std::vector<int>& task_times, const std::vector<std::vector<int>>& raw_precedence);
#endif //VDLS_H
