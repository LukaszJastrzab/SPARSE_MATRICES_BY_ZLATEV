cmake_minimum_required( VERSION 3.19 )

project( SparseMatricesByZlatev LANGUAGES CXX )

if( NOT EXISTS "${CMAKE_SOURCE_DIR}/build" )
    file( MAKE_DIRECTORY "${CMAKE_SOURCE_DIR}/build" )
endif()

set( CMAKE_CXX_STANDARD 17)

include( FetchContent )

FetchContent_Declare(
  googletest
  GIT_REPOSITORY https://github.com/google/googletest.git
  GIT_TAG v1.14.0
)

set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)

FetchContent_MakeAvailable( googletest )

enable_testing()

add_library(sparse_matrices_lib
    src/sparse_matrices.cpp
    src/sparse_matrices.hpp
)

add_executable( ${PROJECT_NAME}
    tests/main.cpp
)

set( EXE_NAME ${PROJECT_NAME}.exe )

target_link_libraries( ${PROJECT_NAME}
    sparse_matrices_lib
    GTest::gtest
    GTest::gtest_main
)

include( GoogleTest )
gtest_discover_tests( ${PROJECT_NAME} )
