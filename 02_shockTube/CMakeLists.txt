cmake_minimum_required(VERSION 3.0.0)
project(shockTube CXX)

if(CMAKE_CXX_COMPILER_ID MATCHES "Clang")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fPIC -std=c++17 -O0 -Wall -Wextra")
  set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -O0 -g")
elseif(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fPIC -std=c++17 -O0 -Wall -Wextra")
  set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -O0 -g")
endif()

list(APPEND FLOW_SOURCE_FILES
     ${CMAKE_CURRENT_SOURCE_DIR}/src/ShockTubeSolver.cpp)

list(APPEND FLOW_HEADER_FILES
     ${CMAKE_CURRENT_SOURCE_DIR}/include/ShockTubeSolver.hpp)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

set(FLOW_INCLUDE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/include")

add_executable(shockTube ${FLOW_HEADER_FILES} ${FLOW_SOURCE_FILES})
target_include_directories(shockTube PRIVATE ${FLOW_INCLUDE_DIR})
