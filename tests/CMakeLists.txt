cmake_minimum_required(VERSION 3.18)
include(FetchContent)
FetchContent_Declare(googletest
  URL https://github.com/google/googletest/archive/refs/tags/v1.14.0.tar.gz
  DOWNLOAD_EXTRACT_TIMESTAMP TRUE)
FetchContent_MakeAvailable(googletest)
include(GoogleTest)

foreach(t sequence matching align io refdiff pairing)
  add_executable(test_${t} test_${t}.cpp)
  target_link_libraries(test_${t} PRIVATE g2g_lib GTest::gtest_main)
  gtest_discover_tests(test_${t})
endforeach()

add_executable(test_decompose test_decompose.cpp)
target_link_libraries(test_decompose PRIVATE g2g_lib GTest::gtest_main)
gtest_discover_tests(test_decompose)

if(G2G_HAVE_OPENMP)
  add_executable(test_parallel test_parallel.cpp)
  target_link_libraries(test_parallel PRIVATE g2g_lib GTest::gtest_main)
  gtest_discover_tests(test_parallel)
endif()

# -- Phase 5 correctness tests ----------------------------------
add_executable(test_phase5_exp1 test_phase5_exp1_synthetic_gt.cpp)
target_link_libraries(test_phase5_exp1 PRIVATE g2g_lib GTest::gtest_main)
gtest_discover_tests(test_phase5_exp1)

add_executable(test_phase5_exp2 test_phase5_exp2_minimap2.cpp)
target_link_libraries(test_phase5_exp2 PRIVATE g2g_lib GTest::gtest_main)
gtest_discover_tests(test_phase5_exp2)
