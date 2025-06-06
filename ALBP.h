#ifndef ALBP_H
#define ALBP_H

#include <string>
#include <vector>

struct PrecedenceRelation {
    int parent, child;
};




struct ALBP{
    std::string name; //Name of the ALBPß
    int C; //Cycle time
    int N; //Number of tasks
    int S; //Number of stations (not used for SALBP-1)
    std::vector<int> task_time; //Task times
    std::vector<int> prec_mat; //Precedence matrix
    std::vector<int> t_close_mat;
    std::vector<std::vector<int>> dir_suc;
    std::vector<std::vector<int>> dir_pred;
    std::vector<std::vector<int>> suc;
    std::vector<std::vector<int>> pred;
    std::vector<PrecedenceRelation> precedence_relations; //Precedence relations
    std::vector<int> task_assignment; //Original task assignment (if applicable)


    ALBP() = default;
    explicit ALBP(const std::string& filename) { loadFromFile(filename); }
    ALBP(int C, int N, const std::vector<int>& task_times, const std::vector<std::vector<int>>& raw_precedence);
    void print(bool print_prec_mat);

    void calc_trans_closure();

    bool loadFromFile(const std::string& filename);
};

#endif // ALBP_H