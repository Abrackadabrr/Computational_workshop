add_executable(integration Lab1/main.cpp)
add_executable(fvm Lab2/fvm.cpp)
add_executable(fem Lab2/fem.cpp
        Lab3/dynamic_problem.cpp)
add_executable(dynamic_fem Lab3/dynamic_problem.cpp)

find_package(OpenMP)
if (OPENMP_FOUND)
    message("OpenMP Found")
endif ()

target_link_libraries(integration PUBLIC ${PROJECT_NAME})
target_include_directories(integration PUBLIC ../src/)
target_link_libraries(fvm PUBLIC ${PROJECT_NAME} OpenMP::OpenMP_CXX)
target_include_directories(fvm PUBLIC ../src/)
target_link_libraries(fem PUBLIC ${PROJECT_NAME} OpenMP::OpenMP_CXX)
target_include_directories(fem PUBLIC ../src/)
target_link_libraries(dynamic_fem PUBLIC ${PROJECT_NAME} OpenMP::OpenMP_CXX)
target_include_directories(dynamic_fem PUBLIC ../src/)
