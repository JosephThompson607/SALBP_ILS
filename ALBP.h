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
    std::vector<int> task_time; //Task times
    std::vector<int> prec_mat; //Precedence matrix
    std::vector<int> t_close_mat;
    std::vector<std::vector<int>> dir_suc;
    std::vector<std::vector<int>> dir_pred;
    std::vector<PrecedenceRelation> precedence_relations; //Precedence relations
    std::vector<int> task_assignment; //Task assignment
    void print(bool print_prec_mat);
    bool loadFromFile(const std::string& filename);
};

#endif // ALBP_H