#ifndef ALBP_H
#define ALBP_H

#include <string>
#include <vector>

struct PrecedenceRelation {
    int parent, child;
};




struct ALBP{
    std::string name; //Name of the ALBPß
    int C; //Cycle time (for SALBP-1)
    int N; //Number of tasks
    int S; //Number of stations (For SALBP-2)
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
    static ALBP type_1(int C, int N, const std::vector<int>& task_times, const std::vector<std::vector<int>>& raw_precedence, bool reverse=false);
    static ALBP type_2(int S, int N, const std::vector<int>& task_times_, const std::vector<std::vector<int>>& raw_precedence, bool reverse=false);
    [[nodiscard]] ALBP reverse() const;
    void print(bool print_prec_mat);

    void calc_trans_closure();

    bool loadFromFile(const std::string& filename);
    private:
        ALBP(int C_, int S_, int N_,
            const std::vector<int>& task_times_,
            const std::vector<std::vector<int>>& raw_precedence,
            bool reverse);

        ALBP(int C_, int S_, int N_,
             const std::vector<int>& task_times_,
             const std::vector<PrecedenceRelation>& raw_precedence,
             bool reverse);
        void initialize_precedence(int C_, int S_, int N_,
                           const std::vector<int>& task_times_,
                           bool reverse);

        void add_relation(int u, int v, bool reverse);

};

#endif // ALBP_H