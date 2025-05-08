
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
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            t_close_mat[i * N + j] = prec_mat[i * N + j];
        }
    }

    for (int k = 0; k < N; ++k) {
        for (int i = 0; i < N; ++i) {
            for (int j = 0; j < N; ++j) {
                t_close_mat[i * N + j] = t_close_mat[i * N + j] || (t_close_mat[i * N + k] && t_close_mat[k * N + j]);
            }
        }
    }
    return t_close_mat;
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
    enum Section { None, NumTasks, CycleTime, TaskTimes, Precedences } section = None;

    while (std::getline(infile, line)) {
        if (line.empty()) continue;
        if (line == "<number of tasks>") {
            section = NumTasks;
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
                }
            }
        }
    }

    infile.close();
    // Calculate transitive closure
    t_close_mat = transitive_closure(prec_mat, N);
    return true;
}


