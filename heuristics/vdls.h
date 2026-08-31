//
// Created by Joseph Thompson on 2025-07-02.
//

#ifndef VDLS_H
#define VDLS_H
#include <chrono>
#include <utility>
#include <optional>
#include <random>

#include "../ALBP.h"
#include "../albp_solution.h"



class VDLS {
public:
    //explicit ALBP(const std::string& filename) { loadFromFile(filename); }
    explicit VDLS(const ALBP& albp,  int max_attempts, double time_limit, std::optional<unsigned> seed = std::nullopt):
    albp_(albp),
    max_attempts_(max_attempts),
    time_limit_(time_limit),
    rng_(seed.value_or(std::random_device{}()))
    {
        start_time_ = std::chrono::steady_clock::now();

    }
    ALBPSolution solve_type_1( );
    ALBPSolution solve_type_2( );
    bool time_exceeded() const;
    void add_init_solution(std::vector<int>init_solution);
private:
    //Variables
    const ALBP& albp_;
    ALBPSolution best_{albp_.N};
    int n_attempts_=0;
    std::mt19937 rng_;
    int seed_;
    int max_attempts_;

    int max_depth_=6; //number of shifts to consider
    int lb_;
    int n_perts_ = 6;

    std::chrono::steady_clock::time_point start_time_;
    std::chrono::duration<double> time_limit_;
    std::vector<int>  pw_{};

    //Methods
    ALBPSolution vdls_heuristic(  int n_stations,   int lb);
    ALBPSolution hoff_search(int n_stations) const;
    int random_selection(const std::vector<int>& int_vec);
    bool local_search(ALBPSolution& local_best, ALBPSolution& prev_sol, int depth, int last_task=-1);
    void perform_shift(ALBPSolution &sol, int task, int task_idx, int old_station, int new_station);
    void perturbation(ALBPSolution& incumbent_solution);
};

ALBPSolution vdls_solve_salbp1(const ALBP &albp, const std::vector<int> &initial_solution = std::vector<int>(), std::optional<int> max_attempts = std::nullopt, std::optional<double> time_limit = std::nullopt, std::optional<unsigned> seed = std::nullopt);

ALBPSolution vdls_solve_salbp1(int C, int N,
    const std::vector<int>& task_times,
    const std::vector<std::vector<int>>& raw_precedence,
    const std::vector<int> &initial_solution = std::vector<int>(),
    std::optional<int> max_attempts = std::nullopt, std::optional<double> time_limit = std::nullopt, std::optional<unsigned> seed = std::nullopt);
ALBPSolution vdls_solve_salbp2(const ALBP &albp, const std::vector<int> &initial_solution = std::vector<int>(), std::optional<int> max_attempts = std::nullopt, std
                               ::optional<double> time_limit = std::nullopt, std::optional<unsigned> seed = std::nullopt);
ALBPSolution vdls_solve_salbp2(int S, int N, const std::vector<int>& task_times, const std::vector<std::vector<int>>& raw_precedence, const std::vector<int> &initial_solution = std::vector<int>(), std::optional<int> max_attempts = std::nullopt,
    std::optional<double> time_limit = std::nullopt, std::optional<unsigned> seed = std::nullopt);

#endif //VDLS_H
