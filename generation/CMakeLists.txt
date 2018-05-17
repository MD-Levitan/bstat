cmake_minimum_required(VERSION 3.8)
project(belt_test C)

include_directories(include)

set(GCC_COVERAGE_LINK_FLAGS "-lm")
set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${GCC_COVERAGE_LINK_FLAGS}")
set(CMAKE_C_STANDARD 99)

set(SOURCE_FILES main.c generation.c)
add_executable(generation ${SOURCE_FILES})


target_link_libraries(generation bstat)

install(TARGETS generation
        DESTINATION ${PROJECT_BINARY_DIR}
        PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE GROUP_READ GROUP_EXECUTE WORLD_READ WORLD_EXECUTE)