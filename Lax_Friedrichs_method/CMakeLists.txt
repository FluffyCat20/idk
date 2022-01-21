cmake_minimum_required (VERSION 3.12)
project (lax_friedrichs)
add_executable(1dimension main_1d.cpp data_1d.h data_1d.cpp exact_solution.cpp exact_solution.h riemann_v2.c)
add_executable(2dimensions main_2d.cpp data_2d.h data_2d.cpp)

set (nlohmann_json_DIR json_parser/)
find_package(nlohmann_json 3.2.0 REQUIRED)
target_link_libraries(1dimension PRIVATE nlohmann_json::nlohmann_json)
target_link_libraries(2dimensions PRIVATE nlohmann_json::nlohmann_json)
