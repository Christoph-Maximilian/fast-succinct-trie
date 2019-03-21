# TODO: remove those compile flask before benchmarking (those are for optimizing code - not good for debugging)
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG}")
set(CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG}")

add_library(FST SHARED FST.cpp bitmap-rank.cc bitmap-rankF.cc bitmap-select.cc)
