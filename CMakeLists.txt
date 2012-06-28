cmake_minimum_required(VERSION 2.8)
project(baoab-constrained CXX)
enable_language(CXX C)

set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)

find_package(GSL REQUIRED)
include_directories(${GSL_INCLUDE_DIRS})

find_package(Gnuplot)
add_definitions(-DGNUPLOT_EXECUTABLE="${GNUPLOT_EXECUTABLE}")

find_package(OpenMP)

#find_package(Threads)
#find_package(GTest)

find_package(Armadillo)
include_directories(${ARMADILLO_INCLUDE_DIRS})

# add_definitions(-DHAVE_SSE)

set(CXX_WARN_FLAGS "-W -Wall -Wextra -Wconversion -Wshadow -Wpointer-arith -Wcast-qual -Wcast-align -Wwrite-strings -Wfloat-equal")
set(CMAKE_CXX_FLAGS_DEBUG "-O2")
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -funroll-loops -march=native")
set(CMAKE_CXX_FLAGS "-ggdb3 ${CXX_WARN_FLAGS} ${OpenMP_CXX_FLAGS}")
# set(CMAKE_CXX_FLAGS "-ggdb3 -std=c++0x ${CXX_WARN_FLAGS} ${OpenMP_CXX_FLAGS}")

set(LINK_SEARCH_END_STATIC "ON")

add_library(baoab
  BAOAB.h BAOAB.cpp
  Plotter.h Plotter.cpp
  Histogram.h Histogram.cpp
  linspace.h
  Experiment.h Experiment.cpp)

add_executable(simul simul.cpp Options.cpp)
target_link_libraries(simul baoab
 ${GSL_LIBRARIES} ${ARMADILLO_LIBRARIES})

# enable_testing()
# include_directories(${GTEST_INCLUDES})
# add_executable(test-suite test-suite.cpp)
# target_link_libraries(test-suite
#   ${GSL_LIBRARIES} ${ARMADILLO_LIBRARIES}
#   ${GTEST_BOTH_LIBRARIES} ${CMAKE_THREAD_LIBS_INIT})
# add_test(test-suite test-suite)
