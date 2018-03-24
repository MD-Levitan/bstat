include_directories(include)
#add_subdirectory(./src)
set(GCC_COVERAGE_LINK_FLAGS "-lm")
set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${GCC_COVERAGE_LINK_FLAGS}")

add_executable(cmm_test cmm_test.c)
target_link_libraries(cmm_test bstat)
add_test(cmm_test ${PROJECT_BINARY_DIR}/tests/cmm_test)

add_executable(cmms_test cmms_test.c)
target_link_libraries(cmms_test bstat)
add_test(cmms_test ${PROJECT_BINARY_DIR}/tests/cmms_test)


add_executable(hmm_test hmm_test.c)
target_link_libraries(hmm_test bstat)
add_test(hmm_test ${PROJECT_BINARY_DIR}/tests/hmm_test)

add_executable(dcmm_test dcmm_test.c)
target_link_libraries(dcmm_test bstat)
add_test(dcmm_test ${PROJECT_BINARY_DIR}/tests/dcmm_test)

#add_executable(cmm_test cmm_test.c)

