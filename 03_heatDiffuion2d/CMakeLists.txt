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
  MESH_SOURCE_FILES
  ${CMAKE_CURRENT_SOURCE_DIR}/mesh/src/SU2meshparser.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/mesh/src/Node2d.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/mesh/src/Face2d.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/mesh/src/CellQuad4.cpp)

list(
  APPEND
  MESH_HEADER_FILES
  ${CMAKE_CURRENT_SOURCE_DIR}/mesh/inc/SU2meshparser.hpp
  ${CMAKE_CURRENT_SOURCE_DIR}/mesh/inc/Node2d.hpp
  ${CMAKE_CURRENT_SOURCE_DIR}/mesh/inc/Face2d.hpp
  ${CMAKE_CURRENT_SOURCE_DIR}/mesh/inc/CellQuad4.hpp)

list(APPEND HEAT_SOURCE_FILES
     ${CMAKE_CURRENT_SOURCE_DIR}/heatdiffusion2d/src/HeatDiffusion2d.cpp)

list(APPEND HEAT_HEADER_FILES
     ${CMAKE_CURRENT_SOURCE_DIR}/heatdiffusion2d/inc/HeatDiffusion2d.hpp)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

set(MESH_INCLUDE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/mesh/inc")
set(HEAT_INCLUDE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/heatdiffusion2d/inc")

add_library(heatdiffusion2d SHARED)
target_sources(
  heatdiffusion2d PRIVATE ${MESH_SOURCE_FILES} ${MESH_HEADER_FILES}
                          ${HEAT_SOURCE_FILES} ${HEAT_HEADER_FILES})
target_include_directories(heatdiffusion2d PUBLIC ${MESH_INCLUDE_DIR}
                                                  ${HEAT_INCLUDE_DIR})

add_executable(main ${CMAKE_CURRENT_SOURCE_DIR}/driver/main.cpp)
target_link_libraries(main PRIVATE heatdiffusion2d)
