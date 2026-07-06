//
// Created by Joseph Thompson on 2026-06-30.
// Simple-Tabu search from Pape 2015, "Heuristics and lower bounds for the simple assembly line balancing problem type 1"
//

#ifndef TABU_H
#define TABU_H


#include <deque>
#include <random>


#include "salbp_basics.h"
#include "../ALBP.h"
#include "../albp_solution.h"



class TabuList {

public:
    explicit TabuList( int n, size_t max_size);
    bool contains(const std::pair<int,int> & move);
    void insert(const std::pair<int,int> & move);
    void reset();
private:
    [[nodiscard]] size_t index(const std::pair<int,int> & move) const;
    int n_;
    size_t max_size_;
    std::vector<bool> is_tabu_;
    std::deque<std::pair<int,int>> tabu_list_;
    bool del_move(const std::pair<int, int> &move);
};

struct TabuMove {
    ALBPSolution sol;
    int obj;
    std::pair<int,int> reversal;
    std::pair<int,int> move;
    int task_idx;
    std::optional<std::pair<int,int>> reversal2=std::nullopt;
    std::optional<std::pair<int,int>> move2=std::nullopt;
    std::optional<int> task2_idx=std::nullopt;
};

class Tabu {

    public:
        explicit Tabu(const ALBP& albp, double time_limit,std::optional<unsigned> seed = std::nullopt):
        albp_(albp),
        lb_(calc_salbp_1_bin_lbs(albp_.task_time, albp.C)),
        max_ovlo_station_(-1),
        tabu_(albp_.N, 10),
        time_limit_(time_limit),
        rng_(seed.value_or(std::random_device{}())){};
        ALBPSolution solve();
        [[nodiscard]] ALBPSolution best() const {return best_;}

    private:
        const ALBP& albp_;
        int lb_;
        TabuList tabu_;
        size_t max_ovlo_station_;
        ALBPSolution best_{albp_.N};
        ALBPSolution current_{albp_.N};
        std::chrono::steady_clock::time_point start_time_;
        std::chrono::duration<double> time_limit_;
        std::mt19937 rng_;

        bool is_optimal(ALBPSolution &sol) ;
        int roulette_select(ALBPSolution &sol);
        void add_init_solution(std::vector<int>init_solution);

        int try_shift(const ALBPSolution &sol, int task, int old_station, int new_station) const;

        int try_swap(const ALBPSolution &sol, int task_1, int task_2, int station_1, int station_2) const;

        void perform_shift(ALBPSolution &sol, int task, int task_idx, int old_station,int new_station);
        void perform_swap(ALBPSolution &sol, int task_1, int task_idx_1, int task_2, int task_idx_2, int station_1,int station_2);
        void elim_station(ALBPSolution &sol);

        TabuMove shift(int s, int &best_obj, const ALBPSolution &current);

        TabuMove swap(int s, int &local_obj, const ALBPSolution &current);

        void do_move(ALBPSolution &sol, const TabuMove &move, bool update_tabu);

        void shift_and_swap(ALBPSolution& sol, int  s);
        void quad_overloads(ALBPSolution &sol,  int C);
        bool time_exceeded() const;

    public:

};

ALBPSolution tabu_solve_salbp1(const ALBP &albp , std::optional<double> time_limit=std::nullopt);
ALBPSolution tabu_solve_salbp1(int C, int N,
             const std::vector<int> &task_times,
             const std::vector<std::vector<int> > &raw_precedence, std::optional<double> time_limit = std::nullopt);
#endif //TABU_H
