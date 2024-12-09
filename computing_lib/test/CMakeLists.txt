find_package(GTest REQUIRED)
if (GTEST_FOUND)
    message("GTest Found")
endif ()

find_package(Eigen3 REQUIRED)

message("Tests building has begun")

file(GLOB_RECURSE TESTS test_*.cpp)

foreach (TEST IN LISTS TESTS)
    get_filename_component(FILE_NAME ${TEST} NAME_WE)
    set(EXE_NAME run_${FILE_NAME})
    add_executable(${EXE_NAME} ${TEST}
            mesh/test_INMOST_understanding.cpp)
    message("TEST: ${EXE_NAME}")
    target_include_directories(${EXE_NAME} PUBLIC ../src)
    target_link_libraries(${EXE_NAME} ${PROJECT_NAME} GTest::gtest_main computing_lib Eigen3::Eigen)
    target_compile_options(${EXE_NAME} PRIVATE -Werror)
    add_test(${EXE_NAME} ${EXE_NAME})
endforeach ()

