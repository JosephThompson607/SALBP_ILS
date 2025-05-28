cmake_minimum_required(VERSION 3.31)
project(ILS_ALBP)
#
set(CMAKE_PREFIX_PATH "/Users/letshopethisworks2/opt/anaconda3/envs/torch_env" CACHE PATH "Prefix for finding Anaconda packages")

find_package(Python COMPONENTS Interpreter Development REQUIRED)


set(CMAKE_CXX_STANDARD 26)  # Or whatever C++ standard you're using
# Add pybind11 from submodule
add_subdirectory(extern/pybind11/)

message(STATUS "Found Python interpreter: ${Python_EXECUTABLE}")
message(STATUS "Found Python include directories: ${Python_INCLUDE_DIRS}")
message(STATUS "Found Python libraries: ${Python_LIBRARIES}")
#
#
## Add executable target with source files
#add_executable(ILS_ALBP
#        main.cpp
#        ALBP.cpp
#        albp_solution.cpp
#        ils.cpp
#        ils.h
#       bindings.cpp
#)



## Add sources
set(SOURCES
        bindings.cpp
        ALBP.cpp
)

## Create the Python module
pybind11_add_module(ILS_ALBP ${SOURCES})