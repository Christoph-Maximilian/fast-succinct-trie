cmake_minimum_required(VERSION 2.6)
project(fst)

message(STATUS "Configuring..." ${CMAKE_PROJECT_NAME})


set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g2   -mpopcnt -pthread -std=c++14") #-Werror -O9

# TODO: remove those compile flask before benchmarking (those are for optimizing code - not good for debugging)
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} ")
set(CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG} ")

include_directories("${CMAKE_CURRENT_SOURCE_DIR}/include")

add_subdirectory(src)
add_subdirectory(test)

#add_library(fstlib SHARED
#        src/FST.cpp
#        src/bitmap-rank.cc
#        src/bitmap-rankF.cc
#        src/bitmap-select.cc
#        include/bitmap-rank.h
#        include/bitmap-rankF.h
#        include/bitmap-select.h
#        include/common.h
#        include/FST.hpp
#        include/popcount.h
#        include/shared.h
#        )

#target_include_directories(fstlib PUBLIC ..)
