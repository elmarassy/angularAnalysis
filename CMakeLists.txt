cmake_minimum_required(VERSION 3.28)
project(CERN)

set(CMAKE_CXX_STANDARD 17)

find_package(ROOT REQUIRED COMPONENTS RIO Net Hist Tree RooFit RooFitCore)

include(${ROOT_USE_FILE})

include(FetchContent)
FetchContent_Declare(
        pybind11
        GIT_REPOSITORY https://github.com/pybind/pybind11.git
        GIT_TAG        v2.12.0
)
FetchContent_MakeAvailable(pybind11)

add_executable(cern src/main.cpp
        include/fit.h
        include/timeDependentFit.h
        src/fit.cpp
        src/timeDependentFit.cpp
        include/timeDependentFit.h
        src/newTimeDependentFit.cpp
        include/newTimeDependentFit.h
        include/parameterValues.h)

target_link_libraries(cern ${ROOT_LIBRARIES})

target_include_directories(cern PRIVATE ${ROOT_INCLUDE_DIRS})
