# TODO: remove those compile flask before benchmarking (those are for optimizing code - not good for debugging)

if (NOT TARGET sdsl)
    add_subdirectory(../../../sdsl-lite sdsl)
endif()

set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -
O0")
set(CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG} -O0")

add_library(FST SHARED FST.cpp bitmap-rank.cc bitmap-rankF.cc bitmap-select.cc)
target_link_libraries(FST sdsl)
