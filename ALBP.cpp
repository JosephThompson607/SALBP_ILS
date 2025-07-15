
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

ALBP::ALBP(int C_, int S_, int N_, const std::vector<int>& task_times_, const std::vector<std::vector<int>>& raw_precedence) {
    if (task_times_.size() != static_cast<size_t>(N_)) {
        throw std::invalid_argument("task_times size does not match N");
    }
    C = C_;
    N = N_;
    S= S_; //Not used for SALBP_1
    task_time = task_times_;
    name = "constructed_from_data";

    // Initialize empty vectors for size
    prec_mat.resize(N * N, 0);
    dir_suc.resize(N);
    dir_pred.resize(N);
    suc.resize(N);
    pred.resize(N);
    for (const auto& pair : raw_precedence) {
        const int u = pair[0];
        const int v = pair[1];
        if (u < 1 || u > N || v < 1 || v > N) {
            std::cerr << "Invalid precedence pair. Assuming 1 indexed based pairs: (" << u << ", " << v << ")\n";
            continue;
        }
        precedence_relations.push_back({u, v});
        prec_mat[(u - 1) * N + (v - 1)] = 1;
        dir_suc[u - 1].push_back(v-1);
        dir_pred[v-1].push_back(u-1);
    }
    calc_trans_closure();
}

ALBP ALBP::type_2(int S_, int N_, const std::vector<int>& task_times_, const std::vector<std::vector<int>>& raw_precedence) {
    int C_ub = std::accumulate(task_times_.begin(), task_times_.end(), 0);
   return ALBP(C_ub, S_, N_, task_times_, raw_precedence);
}

ALBP ALBP::type_1(int C_, int N_, const std::vector<int>& task_times_, const std::vector<std::vector<int>>& raw_precedence) {
    int S_ub = N_;
    return ALBP(C_, S_ub, N_, task_times_, raw_precedence);
}



//print function
void ALBP::print(bool print_prec_mat = false) {
    std::cout << "ALBP Name: " << name << std::endl;
    std::cout << "Cycle time: " << C << std::endl;
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


