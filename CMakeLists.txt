cmake_minimum_required(VERSION 2.8)
project(baoab-constrained CXX)
enable_language(CXX C)

set(Boost_USE_STATIC_LIBS    ON)
set(Boost_USE_MULTITHREADED  OFF)
set(Boost_USE_STATIC_RUNTIME ON)
find_package(Boost)
# find_package(Threads)
# find_package(OpenMP)
find_package(Gnuplot)
find_package(Armadillo)

include_directories(${ARMADILLO_INCLUDE_DIRS})

add_definitions(-DGNUPLOT_EXECUTABLE="${GNUPLOT_EXECUTABLE}")

set(CXX_WARN_FLAGS "-W -Wall -Wextra -Wconversion -Wshadow -Wpointer-arith -Wcast-qual -Wcast-align -Wwrite-strings -Wfloat-equal")
set(CMAKE_CXX_FLAGS_DEBUG "-O0")
set(CMAKE_CXX_FLAGS_RELEASE "-DNDEBUG -Ofast -march=native -mtune=native -msse4.2 -mssse3")
set(CMAKE_CXX_FLAGS "-ggdb3 -std=c++0x ${CXX_WARN_FLAGS} ${OpenMP_CXX_FLAGS}")

set(LINK_SEARCH_END_STATIC "ON")

add_executable(baoab baoab.cpp Plotter.cpp dop853.c)
target_link_libraries(baoab ${Boost_LIBRARIES})

# add_executable(baoab-backup baoab-backup.cpp Plotter.cpp)
# target_link_libraries(baoab-backup ${Boost_LIBRARIES})

# add_executable(rattle rattle.cpp Plotter.cpp)
# target_link_libraries(rattle ${Boost_LIBRARIES})

# add_executable(ode-test ode-test.cpp dop853.c)
# target_link_libraries(ode-test ${Boost_LIBRARIES})

# add_executable(O O.cpp)
# target_link_libraries(O ${Boost_LIBRARIES} ${ARMADILLO_LIBRARIES})
