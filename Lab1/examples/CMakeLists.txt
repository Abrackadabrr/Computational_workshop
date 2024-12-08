add_executable(integration main.cpp)
add_executable(fvm Lab2/fvm.cpp
        Lab2/TestFucntions.hpp)

find_package(OpenMP)
if (OPENMP_FOUND)
    message("OpenMP Found")
endif ()

target_link_libraries(integration PUBLIC ${PROJECT_NAME})
target_include_directories(integration PUBLIC ../src/)
target_link_libraries(fvm PUBLIC ${PROJECT_NAME} OpenMP::OpenMP_CXX)
target_include_directories(fvm PUBLIC ../src/)