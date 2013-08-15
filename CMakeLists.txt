cmake_minimum_required(VERSION 2.8)

project(baoab-constrained CXX)
enable_language(CXX C)

set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)

include_directories(${PROJECT_SOURCE_DIR}/contrib/cubature)

find_package(OpenMP)

find_package(Gnuplot)
add_definitions(-DGNUPLOT_EXECUTABLE="${GNUPLOT_EXECUTABLE}" -D_GNU_SOURCE)

find_package(Armadillo)
include_directories(${ARMADILLO_INCLUDE_DIRS})

set(CXX_WARN_FLAGS "-W -Wall -Wextra -Wconversion -Wshadow -Wpointer-arith -Wcast-qual -Wcast-align -Wwrite-strings -Wfloat-equal")
set(CMAKE_CXX_FLAGS_DEBUG "-O0")
set(CMAKE_CXX_FLAGS_RELEASE "-DNDEBUG -O3 -funroll-loops -flto")
set(CMAKE_CXX_FLAGS "-g ${CXX_WARN_FLAGS} ${OpenMP_CXX_FLAGS} -std=c++11")

add_library(baoab BAOAB.h BAOAB.cpp Plotter.h Plotter.cpp linspace.h Experiment.h Experiment.cpp Average.h Average.cpp Force.h Force.cpp)

add_executable(simul simul.cpp Options.cpp)
target_link_libraries(simul baoab ${ARMADILLO_LIBRARIES})
