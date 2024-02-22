cmake_minimum_required(VERSION 3.0 FATAL_ERROR)
project(femtoMerge)

find_package(ROOT REQUIRED COMPONENTS Core RIO Hist)
include_directories(${CMAKE_CURRENT_SOURCE_DIR})

add_library(FemtoMerge SHARED femtoMerge.cxx)
target_link_libraries(FemtoMerge PUBLIC ROOT::Core ROOT::RIO ROOT::Hist)

add_executable(femtoMerge femtoMerge.cxx)
target_link_libraries(femtoMerge FemtoMerge)