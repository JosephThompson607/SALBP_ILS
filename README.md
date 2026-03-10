# SALBP1_heuristics — Build & CLI Reference
## Building
mkdir build && cd build
cmake ..
cmake --build .

The executable will be at build/ILS_ALBP_exe.
To build only the executable without the Python module:

cmake .. -DBUILD_PYTHON_MODULE=OFF
cmake --build .

## Basic Usage
./build/ILS_ALBP_exe <path_to_file.alb> [options]
Running with no arguments triggers a default internal test run.

## Options
- --heuristic <name>        Solver to use (default: priority)
Choices: priority, MHH, hoffman, VDLS, ILS
- --n_stations <int>        Number of stations — switches to SALBP-2 mode
- --time_limit <int>        Time limit in seconds
- --max_attempts <int>      Max solver iterations
- --random_seed <int>       Random seed for reproducibility
- --priority_n_random <int> Number of random solutions for priority heuristic (default: 100)

## Examples

SALBP-1 with default priority heuristic:
./build/ILS_ALBP_exe problem.alb

SALBP-1 with ILS, 30 second limit:
./build/ILS_ALBP_exe problem.alb --heuristic ILS --time_limit 30

SALBP-2 with 5 stations:
./build/ILS_ALBP_exe problem.alb --n_stations 5

SALBP-1 with reproducible priority run:
./build/ILS_ALBP_exe problem.alb --heuristic priority --priority_n_random 200 --random_seed 42


## Python binding

This project also has a python binding that requires Pybind 11. Sample usage can be seen in test_module.py


## References

The heuristics implemented in this project are based on the following works:

- **Hoffman heuristic** — Hoffman, T.R. (1963). *Assembly Line Balancing with a Precedence Matrix*. Management Science, 9(4), 551–562.

- **Multi-Hoffman heuristic (MHH)** — Sternatz, J. (2014). *Enhanced multi-Hoffmann heuristic for efficiently solving real-world assembly line balancing problems in automotive industry*. European Journal of Operational Research, 235(3), 740–754.

- **Variable-Depth Local Search (VDLS)** — Álvarez-Miranda, E., Pereira, J., Vargas, C., & Vilà, M. (2023). *Variable-depth local search heuristic for assembly line balancing problems*. International Journal of Production Research, 61(9), 3103–3121.

- **Iterated Local Search (ILS)** — Ghandi, S., & Masehian, E. (2025). *An efficient solution to the simple assembly line balancing problem type 1 using iterated local search*. Engineering Applications of Artificial Intelligence, 144, 110162.

- **Priority ranking methods** — Scholl, A. (1999). *Balancing and Sequencing of Assembly Lines*. Physica-Verlag HD.