cmake_minimum_required(VERSION 3.22)
project(computing_lib)

add_subdirectory(src)

enable_testing()
add_subdirectory(test)

add_subdirectory(examples)