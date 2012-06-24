cmake_minimum_required(VERSION 2.8)
project(baoab-constrained CXX)
enable_language(CXX C)

set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)

set(Boost_USE_STATIC_LIBS ON)
set(Boost_USE_MULTITHREADED OFF)
set(Boost_USE_STATIC_RUNTIME ON)
find_package(Boost REQUIRED program_options)
include_directories(${Boost_INCLUDE_DIRS})

find_package(GSL REQUIRED)
include_directories(${GSL_INCLUDE_DIRS})

find_package(Gnuplot)
add_definitions(-DGNUPLOT_EXECUTABLE="${GNUPLOT_EXECUTABLE}")

find_package(OpenMP)

#find_package(Threads)
#find_package(GTest)
  
add_definitions(-DHAVE_SSE)

set(CXX_WARN_FLAGS "-W -Wall -Wextra -Wconversion -Wshadow -Wpointer-arith -Wcast-qual -Wcast-align -Wwrite-strings -Wfloat-equal")
set(CMAKE_CXX_FLAGS_DEBUG "-O2")
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -DNDEBUG -march=native")
# set(CMAKE_CXX_FLAGS_RELEASE "-O3 -march=native")
set(CMAKE_CXX_FLAGS "-ggdb3 ${CXX_WARN_FLAGS} ${OpenMP_CXX_FLAGS}")

set(LINK_SEARCH_END_STATIC "ON")

add_library(baoab
  BAOAB.h BAOAB.cpp
  BAOAB_with_Rotation.cpp
  BAOAB_with_RATTLE.cpp
  BAOAB_with_DoPri.cpp
  Plotter.h Plotter.cpp
  Histogram.h Histogram.cpp
  linspace.h double2.h
  dop853.h dop853.c
  dopri5.h dopri5.c
  Experiment.h Experiment.cpp)

add_executable(simul simul.cpp Options.cpp)
target_link_libraries(simul baoab ${Boost_LIBRARIES} ${GSL_LIBRARIES})

# enable_testing()
# include_directories(${GTEST_INCLUDES})
# add_executable(test-suite test-suite.cpp)
# target_link_libraries(test-suite simul
#   ${Boost_LIBRARIES} ${GSL_LIBRARIES}
#   ${GTEST_BOTH_LIBRARIES} ${CMAKE_THREAD_LIBS_INIT})
# add_test(test-suite test-suite)
