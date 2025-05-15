cmake_minimum_required(VERSION 3.10)
project(YourProjectName)

set(CMAKE_CXX_STANDARD 26)  # Or whatever C++ standard you're using

# Add executable target with source files
add_executable(YourProjectName
        main.cpp
        ALBP.cpp
        ils.cpp)