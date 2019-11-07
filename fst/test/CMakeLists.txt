include_directories(${PROJECT_SOURCE_DIR}/src)
find_package(GTest REQUIRED)
include_directories(${GTEST_INCLUDE_DIR})

if (NOT TARGET s2)
  # add FindGFlags.cmake and FindGlog.cmake to cmake module path
  list(APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/../../s2geometry/third_party/cmake")
  add_subdirectory(../../../s2geometry s2geometry)
endif()

function (add_fst_test file_name )
  add_executable(${file_name} ${file_name}.cpp)
  target_link_libraries(${file_name} FST gtest s2)
  add_test(NAME ${file_name} COMMAND ${CMAKE_CURRENT_BINARY_DIR}/${file_name})
endfunction()

add_fst_test(TestFST)

configure_file(bulkload_sort test/bulkload_sort COPYONLY)
configure_file(bulkload_example test/bulkload_example COPYONLY)
configure_file(polygon_ids test/polygon_ids COPYONLY)
configure_file(simple_point_test test/simple_point_test COPYONLY)
configure_file(polygon_ids_simplified test/polygon_ids_simplified COPYONLY)
