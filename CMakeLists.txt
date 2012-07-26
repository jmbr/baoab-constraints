cmake_minimum_required(VERSION 2.8)

project(baoab-constrained CXX)
enable_language(CXX C)

set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)

find_package(OpenMP)

find_package(GSL REQUIRED)
include_directories(${GSL_INCLUDE_DIRS})

find_package(Armadillo REQUIRED)
include_directories(${ARMADILLO_INCLUDE_DIRS})

find_package(Gnuplot)
add_definitions(-DGNUPLOT_EXECUTABLE="${GNUPLOT_EXECUTABLE}")

set(CXX_WARN_FLAGS "-W -Wall -Wextra -Wconversion -Wshadow -Wpointer-arith -Wcast-qual -Wcast-align -Wwrite-strings -Wfloat-equal")
set(CMAKE_CXX_FLAGS_DEBUG "-O2")
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -DNDEBUG -funroll-loops -march=native")
set(CMAKE_CXX_FLAGS "-ggdb3 ${CXX_WARN_FLAGS} ${OpenMP_CXX_FLAGS}")

add_library(baoab SHARED
  BAOAB.h BAOAB.cpp
  Plotter.h Plotter.cpp
  linspace.h
  Experiment.h Experiment.cpp
  Average.h Average.cpp)

add_executable(simul simul.cpp Options.cpp)
target_link_libraries(simul baoab ${GSL_LIBRARIES} ${ARMADILLO_LIBRARIES})
