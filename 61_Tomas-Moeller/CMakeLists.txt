cmake_minimum_required(VERSION 3.4.3)
# Project name and programming language
project(raycaster CXX)

if(CMAKE_CXX_COMPILER_ID MATCHES "Clang")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fPIC -std=c++17 -O0 -Wall -Wextra")
  set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -O0 -g")
elseif(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
  # -fPIC; Create shared library as a Position Independent Code. (be able to
  # implement regardless of memory position)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fPIC -std=c++17 -O0 -Wall -Wextra")
  set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -O0 -g")
endif()

set(INCLUDE_DIR "${PROJECT_SOURCE_DIR}/inc")

add_executable(main ${PROJECT_SOURCE_DIR}/src/TomasMoeller.cpp)
target_include_directories(main PRIVATE ${INCLUDE_DIR})
