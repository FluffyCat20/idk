cmake_minimum_required (VERSION 3.12)
project (lax_friedrichs)
add_executable(lax_friedrichs main.cpp data.h data.cpp exact_solution.cpp exact_solution.h riemann_v2.c)

set (nlohmann_json_DIR json_parser/)
find_package(nlohmann_json 3.2.0 REQUIRED)
target_link_libraries(lax_friedrichs PRIVATE nlohmann_json::nlohmann_json)
