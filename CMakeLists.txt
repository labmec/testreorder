cmake_minimum_required(VERSION 3.14)

project (testreorder LANGUAGES C CXX)


set_property(GLOBAL PROPERTY USE_FOLDERS ON)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
# this option is useful for static analyzers
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)


# Finding the pz package
find_package(NeoPZ REQUIRED HINTS ../NeoPZ_install ../neopz_install)

# creating a common library for all the targets in this repo
add_library(libreorder SHARED "")
target_link_libraries(libreorder PUBLIC NeoPZ::pz)
target_include_directories(libreorder PUBLIC ${PZ_INCLUDE_DIRS})

add_subdirectory(sources)
get_target_property(lib_files libreorder SOURCES)

source_group(TREE ${CMAKE_CURRENT_SOURCE_DIR} FILES ${lib_files})  

# Making pz "seeable" in Xcode
file(GLOB_RECURSE PZ_SRC
     "../neopz/*.h"
     "../neopz/*.cpp"
     "../NeoPzDevelop/*.h"
     "../NeoPzDevelop/*.cpp"
)
add_custom_target( GlobalFiles SOURCES ${PZ_SRC})


# creating the projects
add_subdirectory(targets)
