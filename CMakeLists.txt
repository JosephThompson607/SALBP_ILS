cmake_minimum_required(VERSION 3.10)
project(ILS_ALBP)

set(CMAKE_CXX_STANDARD 26)  # Or whatever C++ standard you're using

# Add executable target with source files
add_executable(ILS_ALBP
        main.cpp
        ALBP.cpp
        albp_solution.cpp
        ils.cpp
        ils.h)