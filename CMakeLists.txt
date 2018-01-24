cmake_minimum_required(VERSION 3.9)
project(bstat C)
set(GCC_COVERAGE_LINK_FLAGS "-lm")
set(CMAKE_C_STANDARD 11)
include_directories(include/)
include_directories(src/)


#add_subdirectory(include)
#add_subdirectory(src)
set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${GCC_COVERAGE_LINK_FLAGS}")
add_executable(bstat main.c include/bstat.h include/hmm.h src/hmm.c src/util.c)