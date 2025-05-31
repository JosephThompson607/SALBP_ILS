cmake_minimum_required(VERSION 3.31)
project(ILS_ALBP)

# Build option
option(BUILD_PYTHON_MODULE "Build Python module" ON)
option(BUILD_EXECUTABLE "Build standalone executable" OFF)

set(CMAKE_CXX_STANDARD 23)

# Common sources
set(COMMON_SOURCES
        ALBP.cpp
        albp_solution.cpp
        ils.cpp
)

if(BUILD_PYTHON_MODULE)
    # Find Python and pybind11
    set(Python3_ROOT_DIR "/Users/letshopethisworks2/miniforge3/envs/ils_env")
    set(Python3_INCLUDE_DIR "/Users/letshopethisworks2/miniforge3/envs/ils_env/include/python3.11")
    set(Python3_LIBRARY "/Users/letshopethisworks2/miniforge3/envs/ils_env/lib/libpython3.11.dylib")

    find_package(Python COMPONENTS Interpreter Development REQUIRED)
    find_package(pybind11 CONFIG REQUIRED)

    # Python module sources
    set(PYTHON_SOURCES
            ${COMMON_SOURCES}
            bindings.cpp
    )

    # Create Python module
    pybind11_add_module(ILS_ALBP ${PYTHON_SOURCES})

    message(STATUS "Building Python module")
endif()

if(BUILD_EXECUTABLE)
    # Executable sources
    set(EXECUTABLE_SOURCES
            ${COMMON_SOURCES}
            main.cpp
    )

    # Create executable
    add_executable(ILS_ALBP_exe ${EXECUTABLE_SOURCES})

    message(STATUS "Building standalone executable")
endif()