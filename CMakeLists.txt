cmake_minimum_required(VERSION 3.27)
project(ILS_ALBP)

# Build option
option(BUILD_PYTHON_MODULE "Build Python module" ON)
option(BUILD_EXECUTABLE "Build standalone executable" ON)

set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
# Common sources
set(COMMON_SOURCES
        ALBP.cpp
        albp_solution.cpp
        ils.cpp
        mhh.cpp
        vdls.cpp
        vdls.h
)

if(BUILD_PYTHON_MODULE)
    # Find Python and pybind11
    set(PYBIND11_FINDPYTHON ON)
include(FetchContent)
    FetchContent_Declare(
        pybind11
        GIT_REPOSITORY https://github.com/pybind/pybind11.git
        GIT_TAG        v2.13.6  # Choose the version you prefer
    )
    FetchContent_MakeAvailable(pybind11)
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