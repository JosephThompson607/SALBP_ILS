#include "ALBP.h"
#include "albp_solution.h"
#include "ils.h"
#include <iostream>
#include <filesystem>

int default_run() {
    ALBP problem;
    namespace fs = std::filesystem;

    if (!problem.loadFromFile("/Users/letshopethisworks2/Documents/phd_paper_material/MMABPWW/SALBP_benchmark/medium data set_n=50/instance_n=50_108.alb")) {
        return 1;
    }
    std::cout << "Loaded " << problem.N << " tasks, cycle time " << problem.C << "\n";

    //runs the local search algorithm
    ALBPSolution result =iterated_local_search(problem, 50000, 0.5);
    std::cout << "Here is the result" << std::endl;
    result.print();
    return 0;
}

int python_constructor_test() {
    int C = 10;
    int N = 5;
    std::vector<int> task_times = {1, 2, 3, 4, 5};

    // Precedence constraints: each pair is (pred, succ), using 1-based indexing
    std::vector<std::vector<int>> precedence = {
        {1, 2},
        {1, 3},
        {2, 4},
        {3, 5}
    };
    std::vector<int> test_assignments = {0,1,2,3,4};
    ALBP albp(C, N, task_times, precedence);
    ALBPSolution result =  ils_solve_SALBP1(C, N, task_times, precedence, 10, 0.5, true, test_assignments);
    std::cout << "Here is the result" << std::endl;
    result.print();
    // std::cout << "Name: " << albp.name << std::endl;
    // std::cout << "Cycle time: " << albp.C << std::endl;
    // std::cout << "Number of tasks: " << albp.N << std::endl;
    //
    // std::cout << "Precedence matrix:" << std::endl;
    // for (int i = 0; i < N; ++i) {
    //     for (int j = 0; j < N; ++j) {
    //         std::cout << albp.prec_mat[i * N + j] << " ";
    //     }
    //     std::cout << std::endl;
    // }
    //
    // std::cout << "Precedence relations:" << std::endl;
    // for (const auto& rel : albp.precedence_relations) {
    //     std::cout << rel.parent << " -> " << rel.child << std::endl;
    // }
    //
    return 0;
}

int main(int argc, char* argv[]) {
    // Check command line arguments
    if (argc != 2) {
        std::cerr << "Usage: " << argv[0] << " <path_to_alb_file>" << std::endl;
        std::cerr << "Example: " << argv[0] << " problem.alb" << std::endl;
        std::cerr << "Performing default run to test system" << std::endl;
        //default_run();
        python_constructor_test();
        return 1;
    }

    std::string filepath = argv[1];

    // Check if file exists
    if (!std::filesystem::exists(filepath)) {
        std::cerr << "Error: File '" << filepath << "' does not exist." << std::endl;
        return 1;
    }

    // Check file extension (optional)
    std::filesystem::path path(filepath);
    if (path.extension() != ".alb") {
        std::cout << "Warning: File does not have .alb extension" << std::endl;
    }

    // Create ALBP instance and load file
    ALBP problem;

    std::cout << "Loading ALBP file: " << filepath << std::endl;

    if (!problem.loadFromFile(filepath)) {
        std::cerr << "Error: Failed to load file '" << filepath << "'" << std::endl;
        return 1;
    }

    // Print success message and basic info
    std::cout << "Successfully loaded ALBP problem!" << std::endl;

    // Print the problem details (assuming you have a print method)
    std::cout << "\n--- Problem Details ---" << std::endl;
    problem.print(false);  // print with precedence matrix

    ALBPSolution result =iterated_local_search(problem, 1000, 0.5);
    std::cout << "Here is the result" << std::endl;
    result.print();
    return 0;

}