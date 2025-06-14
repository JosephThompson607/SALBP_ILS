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
        mhh.cpp
)

if(BUILD_PYTHON_MODULE)
    # Find Python and pybind11
    set(PYBIND11_FINDPYTHON ON)
    find_package(pybind11 CONFIG REQUIRED)
    # Python module sources
    set(PYTHON_SOURCES
            ${COMMON_SOURCES}
            bindings.cpp
    )

    # Create Python modulea
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