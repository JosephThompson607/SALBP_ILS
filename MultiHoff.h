//
// Created by Joseph Thompson on 2025-06-12.
//

#ifndef MULTIHOFF_H
#define MULTIHOFF_H


#include "albp_solution.h"
#include "ALBP.h"
#include <deque>
#include <map>

class MultiHoff {
    public:
        explicit MultiHoff(const ALBP& albp, int max_attempts = 5000);
        ALBPSolution solve( );

private:
    const ALBP& albp_;
    void add_new_available(std::vector<int> &eligible_tasks,  int task);
    void remove_new_available(std::vector<int> &eligible_tasks,  int task);
    void gen_load( int depth, int remaining_time,int start,  float cost, std::vector<int>eligible_tasks) ;
    std::vector<int> filter_eligible(std::vector<int>& elig);

        void calc_salbp_1_lbs(int &lb) const;

        void initialize_current_s_assignments() ;
    int one_packing_search(std::vector<int>& elig, int station);
    bool check_ub() const;

    void mark_task_assigned(int task, std::vector<int>& elig, bool remove_old=true);
    ALBPSolution mhh_sol_;
    std::vector<int> remaining_task_times_;
    std::vector<std::vector<int>> dir_pred_ ;
    std::vector<std::vector<int>> dir_suc_ ;
    std::vector<int> no_prec_tasks_{}; //tasks always available from assignment (forward)
    std::vector<int> no_suc_tasks_{};
    std::vector<int> s_task_assign_{}; //task assignments to a given station
    std::deque<std::vector<int>> s_forwards_{}; //task assignments in forward direction
    std::deque<std::vector<int>> s_backwards_{}; //task assignments in backwards direction
    std::vector<int> best_s_task_assign_{}; //task assignments to a given station
    std::vector<int> n_prec_{}; //number of predecessor unassigned. 0 if available, 1 if not
    std::vector<int> n_suc_{};
    std::vector<int> n_prec_orig_{}; //number of predecessor unassigned. 0 if available, 1 if not
    std::vector<int> n_suc_orig_{};
    std::map<int,  int> reverse_s_assignment_; //{task:(front, back)} note that back is starting counting from right.That won't cause a bug, won't it?
    bool back_pass_=false; //After mf forward station assigments, filling in rest of stations from back
    bool reverse_ = false;  //After filling in stations starting from forward, fill in stations from back (phase 2)
    int n_attempts_;
    int max_attempts_;
    float  min_cost_=0;
    int ub_;
    int lb_;
    int forward_station_;
    int backward_station_;

    std::vector<int>  pw_ranking_{};
    std::vector<int>  rpw_ranking_{};
    //ALBPSolution solve_pass( );
    void reverse_solve_order();
    ALBPSolution solve_one_pass();

};
inline void swap_and_pop(int item, std::vector<int> &vec);
inline void remove_tasks_unordered(std::vector<int>& vec, const std::vector<int>& to_remove);
ALBPSolution mhh_solve_salbp1(int C, int N, const std::vector<int>& task_times, const std::vector<std::vector<int>>& raw_precedence);
ALBPSolution mhh_solve_salbp1(const ALBP &albp );

#endif //MULTIHOFF_H
