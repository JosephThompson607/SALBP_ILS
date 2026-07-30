#ifndef ALBP_H
#define ALBP_H

#include <string>
#include <vector>

struct PrecedenceRelation {
    int parent, child;
};

struct CritPaths {
    std::vector<int> depth;
    std::vector<int> height;
    std::vector<int> depth_pred;
    std::vector<int> height_suc;
};

struct PathStats {
    int total;
    int total_sq;
    int min;
    int max;
    int mean;
    int variance;
};


struct ALBP{
    std::string name; //Name of the ALBPß
    int C; //Cycle time (for SALBP-1)
    int N; //Number of tasks
    int S; //Number of stations (For SALBP-2)
    int total_time;
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
    static ALBP type_1(int C, int N, const std::vector<int>& task_times, const std::vector<std::vector<int>>& raw_precedence, bool reverse=false, bool light=false, bool is_topological=false);
    static ALBP type_2(int S, int N, const std::vector<int>& task_times_, const std::vector<std::vector<int>>& raw_precedence, bool reverse=false, bool light=false, bool is_topological=false);
    [[nodiscard]] ALBP reverse() const;
    void print(bool print_prec_mat);
    void add_precedence_relation(std::vector<int> prec);

    CritPaths get_critical_paths() const;
    PathStats get_path_stats(const std::vector<int>& nodes) const;
    void calc_trans_closure();
    void calc_fast_trans_closure(bool is_topological=false);

    bool loadFromFile(const std::string& filename);
    private:
        ALBP(int C_, int S_, int N_,
            const std::vector<int>& task_times_,
            const std::vector<std::vector<int>>& raw_precedence,
            bool reverse, bool light, bool is_topological);

        ALBP(int C_, int S_, int N_,
             const std::vector<int>& task_times_,
             const std::vector<PrecedenceRelation>& raw_precedence,
             bool reverse, bool light, bool is_topological);
        void initialize_precedence(int C_, int S_, int N_,
                           const std::vector<int>& task_times_,
                           bool reverse);
        void update_prec_and_suc(std::vector<int> new_pred, std::vector<int> new_suc);


        void add_relation(int u, int v, bool reverse);

};
 std::vector<int>fast_transitive_closure(const std::vector<std::vector<int>>& dir_preds, const std::vector<std::vector<int>>& dir_sucs, const std::vector<int>& prec_mat,  bool alreadyTopo = true);


std::vector<int>get_topo_sort(const std::vector<std::vector<int>>& dir_preds, const std::vector<std::vector<int>>& dir_sucs);
#endif // ALBP_H