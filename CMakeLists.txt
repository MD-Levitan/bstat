cmake_minimum_required(VERSION 3.9)
project(bstat C)

option(BUILD_SHARED_LIBS "Build shared libraries." ON)
option(BUILD_TESTS "Build tests." ON)


if(CMAKE_COMPLILER_IS_GNUCC)
    set(GCC_COVERAGE_LINK_FLAGS "-lm")
    set(CMAKE_C_STANDARD 99)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wall")
endif(CMAKE_COMPLILER_IS_GNUCC)

set(dir bin/)
set(${PROJECT_BINARY_DIR} ${dir})
set(${CMAKE_BINARY_DIR} ${dir})
set(${CMAKE_RUNTIME_OUTPUT_DIRECTORY} ${dir})
set(${CMAKE__OUTPUT_DIRECTORY} ${dir})

include_directories(include/)
include_directories(src/)

add_subdirectory(src)

add_library(bstat STATIC
        $<TARGET_OBJECTS:bstat_src>)

if(BUILD_TESTS)
    enable_testing()
    add_subdirectory(tests)
endif()

set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${GCC_COVERAGE_LINK_FLAGS}")
#add_executable(bstat main.c include/bstat.h include/hmm.h include/statistic.h include/cmm.h
#                            src/cmm.c src/hmm.c src/util.c src/statistic.c include/cmms.h src/cmms.c include/dcmm.h src/dcmm.c)
#add_executable(bstat tests/cmm_test.c)