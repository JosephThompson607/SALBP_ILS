
#include "ALBP.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <filesystem>
#include <numeric>


//function that returns a matrix that is the transitiive closure of the precedence matrix
std::vector<int> transitive_closure(const std::vector<int>& prec_mat, int N) {
    std::vector<int> t_close_mat(N * N, 0);
    t_close_mat = prec_mat;

    for (int k = 0; k < N; ++k) {
        for (int i = 0; i < N; ++i) {
            for (int j = 0; j < N; ++j) {
                t_close_mat[i * N + j] = t_close_mat[i * N + j] || (t_close_mat[i * N + k] && t_close_mat[k * N + j]);
            }
        }
    }
    return t_close_mat;
}


std::vector<std::vector<int>> all_successors(const std::vector<int>& t_close_mat, const int N) {
    std::vector<std::vector<int>> all_suc(N);
    for (int i = 0; i <N; ++i) {
        for (int j = 0; j < N; ++j) {
            if (t_close_mat[i * N + j]) {
                all_suc[i ].push_back(j);
            }
        }

    }
    return all_suc;
}
std::vector<std::vector<int>> all_predecessors(const std::vector<int>& t_close_mat, const int N) {
    std::vector<std::vector<int>> all_pred(N);
    for (int j = 0; j < N; ++j) {
        for (int i = 0; i <N; ++i) {
            if (t_close_mat[i * N + j]) {
                all_pred[j ].push_back(i);
            }
        }

    }
    return all_pred;
}

void ALBP::add_relation(int u, int v, bool reverse) {
    if (reverse) std::swap(u, v);

    if (u < 1 || u > N || v < 1 || v > N) {
        std::cerr << "Invalid precedence pair (" << u << ", " << v
                  << "). Assuming 1-indexed.\n";
        return;
    }

    precedence_relations.push_back({u, v});
    prec_mat[(u - 1) * N + (v - 1)] = 1;
    dir_suc[u - 1].push_back(v - 1);
    dir_pred[v - 1].push_back(u - 1);
}

void ALBP::initialize_precedence(int C_, int S_, int N_,
                                 const std::vector<int>& task_times_,
                                 bool reverse) {
    if (task_times_.size() != static_cast<size_t>(N_)) {
        throw std::invalid_argument("task_times size does not match N");
    }
    C = C_;
    N = N_;
    S = S_;
    task_time = task_times_;
    name = "constructed_from_data";

    prec_mat.assign(N * N, 0);
    dir_suc.assign(N, {});
    dir_pred.assign(N, {});
    suc.assign(N, {});
    pred.assign(N, {});
}

ALBP::ALBP(const int C_, const int S_, const int N_,
           const std::vector<int>& task_times_,
           const std::vector<std::vector<int>>& raw_precedence,
           bool reverse) {
    initialize_precedence(C_, S_, N_, task_times_, reverse);

    for (const auto& pair : raw_precedence) {
        if (pair.size() < 2) continue;
        add_relation(pair[0], pair[1], reverse);
    }

    calc_trans_closure();
}
ALBP::ALBP(int C_, int S_, int N_,
           const std::vector<int>& task_times_,
           const std::vector<PrecedenceRelation>& raw_precedence,
           bool reverse) {
    initialize_precedence(C_, S_, N_, task_times_, reverse);

    for (const auto& rel : raw_precedence) {
        add_relation(rel.parent, rel.child, reverse);
    }

    calc_trans_closure();
}

ALBP ALBP::type_2(int S_, int N_, const std::vector<int>& task_times_, const std::vector<std::vector<int>>& raw_precedence, const bool reverse) {
    int C_ub = std::accumulate(task_times_.begin(), task_times_.end(), 0);
   return ALBP(C_ub, S_, N_, task_times_, raw_precedence, reverse );
}

ALBP ALBP::type_1(int C_, int N_, const std::vector<int>& task_times_, const std::vector<std::vector<int>>& raw_precedence, const bool reverse) {
    int S_ub = N_;
    return ALBP(C_, S_ub, N_, task_times_, raw_precedence, reverse);
}



//print function
void ALBP::print(bool print_prec_mat = false) {
    std::cout << "ALBP Name: " << name << std::endl;
    std::cout << "Cycle time (SALBP-1 only: " << C << std::endl;
    std::cout << "Number of stations (SALBP-2 only): " << S << std::endl;
    std::cout << "Number of tasks: " << N << std::endl;
    std::cout << "Task times: ";
    for (int i = 0; i < N; ++i) {
        std::cout << task_time[i] << " ";
    }
    std::cout << std::endl;
    std::cout << "Number of precedence relations: " << precedence_relations.size() << std::endl;
    std::cout << "Precedence relations: ";
    for (const auto& rel : precedence_relations) {
        std::cout << "(" << rel.parent << ", " << rel.child << ") ";
    }
    if (print_prec_mat) {
        std::cout << std::endl << "Precedence matrix: " << std::endl;
        for (int i = 0; i < N; ++i) {
            for (int j = 0; j < N; ++j) {
                std::cout << prec_mat[i * N + j] << " ";
            }
            std::cout << std::endl;
        }
        std::cout << std::endl << "Transitive closure matrix: " << std::endl;
        for (int i = 0; i < N; ++i) {
            for (int j = 0; j < N; ++j) {
                std::cout << t_close_mat[i * N + j] << " ";
            }
            std::cout << std::endl;
        }
    }
    std::cout << std::endl;
}


void ALBP::calc_trans_closure() {
    t_close_mat = transitive_closure(prec_mat, N);
    suc = all_successors(t_close_mat, N);
    pred = all_predecessors(t_close_mat, N);
}


// Load from .alb file; returns true on success
bool ALBP::loadFromFile(const std::string& filename) {
    std::ifstream infile(filename);
    if (!infile.is_open()) {
        std::cerr << "Error: cannot open file " << filename << std::endl;
        return false;
    }
    // 
    std::filesystem::path p(filename);
    name = p.filename().string();
    std::string line;
    enum Section { None, NumTasks, OrderStrength, CycleTime, TaskTimes, Precedences } section = None;

    while (std::getline(infile, line)) {
        // Trim whitespace including \r and \n
        line.erase(line.find_last_not_of(" \t\r\n") + 1);
        line.erase(0, line.find_first_not_of(" \t\r\n"));
        if (line.empty()) continue;
        if (line == "<number of tasks>") {
            section = NumTasks;
            continue;
        }
        if (line == "<order strength>") {
            // Skip order strength for now as it's not in the struct
            section = OrderStrength;
            continue;
        }
        if (line == "<cycle time>") {
            section = CycleTime;
            continue;
        }
        if (line == "<task times>") {
            section = TaskTimes;
            // prepare container
            task_time.assign(N, 0);
            continue;
        }
        if (line == "<precedence relations>") {
            section = Precedences;
            // initialize matrix
            prec_mat.assign(N * N, 0);
            dir_pred.resize(N);
            dir_suc.resize(N);
            continue;
        }
        if (line == "<end>") {
            break;
        }

        std::istringstream iss(line);
        if (section == NumTasks) {
            iss >> N;
        }
        else if (section == CycleTime) {
            iss >> C;
        }
        else if (section == TaskTimes) {
            int id, t;
            if ((iss >> id >> t)) {
                if (id >= 1 && id <= N)
                    task_time[id - 1] = t;
            }
        }
        else if (section == Precedences) {
            // format: u,v
            int u, v;
            char comma;
            if ((iss >> u >> comma >> v) && comma == ',') {
                if (u >= 1 && u <= N && v >= 1 && v <= N) {
                    precedence_relations.push_back({u, v});
                    prec_mat[(u - 1) * N + (v - 1)] = 1;
                    dir_suc[u - 1].push_back(v-1);
                    dir_pred[v-1].push_back(u-1);
                }
            }
        }
    }

    infile.close();
    // Calculate transitive closure
    calc_trans_closure();
    return true;
}

ALBP ALBP::reverse() const {
    ALBP new_salbp = ALBP (C, S, N, task_time, precedence_relations, true);


    return new_salbp;
}


