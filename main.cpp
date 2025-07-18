#include "ALBP.h"
#include "albp_solution.h"
#include "ils.h"
#include "mhh.h"
#include "vdls.h"
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
    ALBP albp = ALBP::type_1(C, N, task_times, precedence);
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
    // } n_attempts_++;
    //
    return 0;
}
int mhh_test() {
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
    ALBP albp = ALBP::type_1(C, N, task_times, precedence);

    ALBPSolution result =  mhh_solve_salbp1(C, N, task_times, precedence);
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

int vdls_salbp_1_test() {
    int C = 20;
    int N = 8;
    std::vector<int> task_times = {11, 17, 9, 5, 8, 12, 10, 3};

    // Precedence constraints: each pair is (pred, succ), using 1-based indexing
    std::vector<std::vector<int>> precedence = {
        {1, 2},
        {2, 3},
        {2, 4},
        {3, 5},
        {3, 6},
        {4, 6},
        {5, 7},
        {6, 8}
    };
    //std::vector<int> test_assignments = {0,1,2,3,4};
    ALBP albp = ALBP::type_1(C, N, task_times, precedence);
    ALBPSolution result =  vdls_solve_salbp1(C, N, task_times, precedence);
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


int vdls_salbp_2_test() {
    int S = 10;
    int N = 20;
    std::vector<int> task_times = {
        132, 120, 514, 190, 209, 457, 163, 491, 503, 138,
        138, 247, 230, 169, 29, 120, 247, 104, 286, 154
    };


    // Precedence constraints: each pair is (pred, succ), using 1-based indexing
    std::vector<std::vector<int>> precedence = {
        {1, 9}, {2, 9}, {3, 9},
        {4, 5}, {4, 6}, {4, 7}, {4, 8},
        {5, 9}, {6, 9},
        {7, 10}, {7, 11}, {7, 16},
        {9, 12}, {9, 13}, {9, 14}, {9, 15},
        {10, 17}, {10, 18},
        {11, 18},
        {12, 16},
        {13, 16}, {13, 17},
        {16, 18},
        {17, 19}, {17, 20}
    };
    //std::vector<int> test_assignments = {0,1,2,3,4};
    ALBP albp = ALBP::type_2(S, N, task_times, precedence);
    ALBPSolution result =  vdls_solve_salbp2(S, N, task_times, precedence, {}, 1000, 2000);
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
    if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << " <path_to_alb_file>" << std::endl;
        std::cerr << "Example: " << argv[0] << " problem.alb" << std::endl;
        std::cerr << "Performing default run to test system" << std::endl;
        //default_run();
        vdls_salbp_2_test();
        return 1;
    }
    bool salbp2 =false;
    std::optional<int> time_limit = std::nullopt;
    std::optional<int> max_attempts = std::nullopt;
    int n_stations;
    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];

        if (arg == "--n_stations" && i + 1 < argc) {
            n_stations = std::stoi(argv[++i]);
            salbp2 = true;

        } else if (arg == "--time_limit" && i + 1 < argc) {
            time_limit = std::stoi(argv[++i]);
        } else if (arg == "--max_attempts") {
            max_attempts = std::stoi(argv[++i]);
        } else {
            std::cerr << "Unknown argument: " << arg << "\n";
        }
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

    // Create ALBP instance and loads file
    ALBP problem;

    std::cout << "Loading ALBP file: " << filepath << std::endl;

    if (!problem.loadFromFile(filepath)) {
        std::cerr << "Error: Failed to loads file '" << filepath << "'" << std::endl;
        return 1;
    }

    // Print success message and basic info
    std::cout << "Successfully loaded ALBP problem!" << std::endl;

    // Optional: parse second argument if provided
    if (salbp2) {
        try {
            problem.S = n_stations;
            std::cout<< "Stations " << problem.S << " detected, performing SALBP-2" << std::endl;
            // Print the problem details (assuming you have a print method)
            std::cout << "\n--- Problem Details ---" << std::endl;
            problem.print(false);  // print without precedence matrix
            std::cout << "\n--- Solving problem with vdls  " << "---"<< std::endl;
            ALBPSolution result = vdls_solve_salbp2(problem, {},max_attempts, time_limit);
            std::cout << "Here is the result" << std::endl;
            result.print();
        } catch (const std::exception& e) {
            std::cerr << "Error: Invalid number of stations: " << argv[2] << std::endl;
            return 1;
        }
    }

    else{ ALBPSolution result =vdls_solve_salbp1(problem);
        std::cout << "Here is the result" << std::endl;
        result.print();
        return 0;}



}