#include "ALBP.h"
#include "ils.h"
#include <iostream>
#include <filesystem>

int main() {
    ALBP problem;
    namespace fs = std::filesystem;

    //
    // std::string path = "../"; // Current directory
    //
    // for (const auto& entry : fs::directory_iterator(path)) {
    //     std::cout << entry.path() << std::endl;
    // }
    if (!problem.loadFromFile("../test_12_task.alb")) {
        return 1;
    }
    std::cout << "Loaded " << problem.N << " tasks, cycle time " << problem.C << "\n";
    
    //runs the local search algorithm
    iterated_local_search(problem, 1, 1.0);
    return 0;
}