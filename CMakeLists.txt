cmake_minimum_required (VERSION 3.4)
project (bcfdelta CXX)

list(APPEND CMAKE_PREFIX_PATH "${CMAKE_CURRENT_SOURCE_DIR}/submodules/b.i.o./build_system")

find_package (bio REQUIRED)

add_executable (bcfdelta main.cpp)
target_link_libraries (bcfdelta bio::bio seqan3::seqan3)

########################################################
## clang-format
########################################################

find_program(CLANG_FORMAT "clang-format-14")

if (CLANG_FORMAT STREQUAL "CLANG_FORMAT-NOTFOUND")
    find_program(CLANG_FORMAT "clang-format")
endif()

if (NOT CLANG_FORMAT STREQUAL "CLANG_FORMAT-NOTFOUND")
    add_custom_target (check_format "find" "${CMAKE_CURRENT_SOURCE_DIR}/" "-maxdepth" "1" "-name" "'*.[ch]pp'" "-exec"
                       ${CLANG_FORMAT} "-style=file" "-n" "-Werror" "{}" "+"
                       COMMENT "Checking the app with clang-format.")

    add_custom_target (do_format "find" "${CMAKE_CURRENT_SOURCE_DIR}/" "-maxdepth" "1" "-name" "'*.[ch]pp'" "-exec"
                       ${CLANG_FORMAT} "-style=file" "-i" "{}" "+"
                       COMMENT "Format the app with clang-format.")
else ()
    message(FATAL_ERROR "clang-format not found. No targets made.")
endif()
