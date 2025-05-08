#include "ALBP.h"
#include "ils.h"
#include <iostream>

int main() {
    ALBP problem;
    if (!problem.loadFromFile("../test.alb")) {
        return 1;
    }
    std::cout << "Loaded " << problem.N << " tasks, cycle time " << problem.C << "\n";
    
    //runs the local search algorithm
    iterated_local_search(problem, 1, 0.5);
    return 0;
}