cmake_minimum_required(VERSION 3.22 FATAL_ERROR)

add_executable(triag main.cpp)

find_package(GTest REQUIRED)
include_directories(${GTEST_INCLUDE_DIRS})

target_link_libraries(triag ${GTEST_LIBRARIES} pthread)