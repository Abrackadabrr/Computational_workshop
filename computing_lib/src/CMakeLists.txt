cmake_minimum_required(VERSION 3.22)
project(computing_lib)

set(CMAKE_LIBRARY_PATH "${CMAKE_LIBRARY_PATH} /usr/local/")

find_package(inmost REQUIRED)
find_package(Eigen3 REQUIRED)

file(GLOB_RECURSE SOURCES *.cpp)

message(STATUS "DIRS: ${INMOST_INCLUDE_DIRS}")
message(STATUS "DEFS: ${INMOST_DEFENITIONS}")
message(STATUS "LIBS: ${INMOST_LIBRARIES}")
message(STATUS "LIB_DIR: ${INMOST_LIBRARY_DIR}")

add_library(${PROJECT_NAME} ${SOURCES}
        fem/assembly/Utils.hpp)

target_compile_definitions(inmost PUBLIC ${INMOST_DEFINITIONS})
target_link_directories(inmost INTERFACE ${INMOST_LIBRARY_DIRS})
target_link_libraries(inmost INTERFACE ${INMOST_LIBRARIES})
target_include_directories(inmost INTERFACE ${INMOST_INCLUDE_DIRS})

target_include_directories(${PROJECT_NAME} PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(${PROJECT_NAME} PUBLIC inmost Eigen3::Eigen)
