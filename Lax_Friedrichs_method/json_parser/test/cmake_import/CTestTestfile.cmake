# CMake generated Testfile for 
# Source directory: D:/json-develop/test/cmake_import
# Build directory: D:/gas_dyn/Lax_Friedrichs_method/json_parser/test/cmake_import
# 
# This file includes the relevant testing commands required for 
# testing this directory and lists subdirectories to be tested as well.
add_test(cmake_import_configure "C:/Qt/Tools/CMake_64/bin/cmake.exe" "-G" "Ninja" "-A" "" "-DCMAKE_CXX_COMPILER=C:/Program Files (x86)/Microsoft Visual Studio/2019/Community/VC/Tools/MSVC/14.28.29333/bin/HostX86/x64/cl.exe" "-Dnlohmann_json_DIR=D:/gas_dyn/Lax_Friedrichs_method/json_parser" "D:/json-develop/test/cmake_import/project")
set_tests_properties(cmake_import_configure PROPERTIES  FIXTURES_SETUP "cmake_import" LABELS "not_reproducible" _BACKTRACE_TRIPLES "D:/json-develop/test/cmake_import/CMakeLists.txt;1;add_test;D:/json-develop/test/cmake_import/CMakeLists.txt;0;")
add_test(cmake_import_build "C:/Qt/Tools/CMake_64/bin/cmake.exe" "--build" ".")
set_tests_properties(cmake_import_build PROPERTIES  FIXTURES_REQUIRED "cmake_import" LABELS "not_reproducible" _BACKTRACE_TRIPLES "D:/json-develop/test/cmake_import/CMakeLists.txt;9;add_test;D:/json-develop/test/cmake_import/CMakeLists.txt;0;")
