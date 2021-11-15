# CMake generated Testfile for 
# Source directory: D:/json-develop/test/cmake_add_subdirectory
# Build directory: D:/gas_dyn/Lax_Friedrichs_method/json_parser/test/cmake_add_subdirectory
# 
# This file includes the relevant testing commands required for 
# testing this directory and lists subdirectories to be tested as well.
add_test(cmake_add_subdirectory_configure "C:/Qt/Tools/CMake_64/bin/cmake.exe" "-G" "Ninja" "-DCMAKE_CXX_COMPILER=C:/Program Files (x86)/Microsoft Visual Studio/2019/Community/VC/Tools/MSVC/14.28.29333/bin/HostX86/x64/cl.exe" "-Dnlohmann_json_source=D:/json-develop" "D:/json-develop/test/cmake_add_subdirectory/project")
set_tests_properties(cmake_add_subdirectory_configure PROPERTIES  FIXTURES_SETUP "cmake_add_subdirectory" LABELS "not_reproducible" _BACKTRACE_TRIPLES "D:/json-develop/test/cmake_add_subdirectory/CMakeLists.txt;1;add_test;D:/json-develop/test/cmake_add_subdirectory/CMakeLists.txt;0;")
add_test(cmake_add_subdirectory_build "C:/Qt/Tools/CMake_64/bin/cmake.exe" "--build" ".")
set_tests_properties(cmake_add_subdirectory_build PROPERTIES  FIXTURES_REQUIRED "cmake_add_subdirectory" LABELS "not_reproducible" _BACKTRACE_TRIPLES "D:/json-develop/test/cmake_add_subdirectory/CMakeLists.txt;8;add_test;D:/json-develop/test/cmake_add_subdirectory/CMakeLists.txt;0;")
