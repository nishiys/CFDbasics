cmake_minimum_required(VERSION 3.0.0)
project(su2meshparser CXX)

if(CMAKE_CXX_COMPILER_ID MATCHES "Clang")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fPIC -std=c++17 -O0 -Wall -Wextra")
  set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -O0 -g")
elseif(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fPIC -std=c++17 -O0 -Wall -Wextra")
  set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -O0 -g")
endif()

list(
  APPEND
  SOURCE_FILES
  ${CMAKE_CURRENT_SOURCE_DIR}/src/SU2meshparser.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/src/Node2d.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/src/Face2d.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/src/CellQuad4.cpp)

list(
  APPEND
  HEADER_FILES
  ${CMAKE_CURRENT_SOURCE_DIR}/inc/SU2meshparser.hpp
  ${CMAKE_CURRENT_SOURCE_DIR}/inc/Node2d.hpp
  ${CMAKE_CURRENT_SOURCE_DIR}/inc/Face2d.hpp
  ${CMAKE_CURRENT_SOURCE_DIR}/inc/CellQuad4.hpp)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

set(INCLUDE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/inc")

add_library(su2meshparser SHARED)
target_sources(su2meshparser PRIVATE ${SOURCE_FILES} ${HEADER_FILES})
target_include_directories(su2meshparser PUBLIC ${INCLUDE_DIR})

add_executable(main ${CMAKE_CURRENT_SOURCE_DIR}/driver/main.cpp)
target_include_directories(main PRIVATE ${INCLUDE_DIR})
target_link_libraries(main PRIVATE su2meshparser)
