cmake_minimum_required (VERSION 2.6)
project (art)

message(STATUS "Configuring..." ${CMAKE_PROJECT_NAME})

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g -O9 -Werror -mpopcnt -pthread -std=c++11")

include_directories("${CMAKE_CURRENT_SOURCE_DIR}/include")

add_subdirectory(src)

