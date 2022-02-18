cmake_minimum_required (VERSION 3.4)
project (bcfdelta CXX)

list(APPEND CMAKE_PREFIX_PATH "${CMAKE_CURRENT_SOURCE_DIR}/submodules/b.i.o./build_system")

find_package (bio REQUIRED)

add_executable (bcfdelta main.cpp)
target_link_libraries (bcfdelta bio::bio seqan3::seqan3)
