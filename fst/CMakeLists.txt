cmake_minimum_required (VERSION 2.6)
project (fst)

message(STATUS "Configuring..." ${CMAKE_PROJECT_NAME})


set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g -O9 -Werror -mpopcnt -pthread -std=c++11")

# TODO: remove those compile flask before benchmarking (those are for optimizing code - not good for debugging)
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -O0")
set(CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG} -O0")

include_directories("${CMAKE_CURRENT_SOURCE_DIR}/include")

add_subdirectory(src)
add_subdirectory(test)
