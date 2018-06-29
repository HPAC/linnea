cmake_minimum_required(VERSION 3.0)
project({0})

find_package(LibTests REQUIRED)

add_executable({0} runner.cpp)
target_link_libraries({0} PRIVATE libtests)
set_target_properties({0} PROPERTIES CXX_STANDARD 14)