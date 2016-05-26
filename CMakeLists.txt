cmake_minimum_required (VERSION 2.6)
project (NARROW_BAND)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
set(EXECUTABLE_OUTPUT_PATH /home/dede/bin/)

find_package(PCL 1.6 REQUIRED)
include_directories(${PCL_INCLUDE_DIRS})
link_directories(${PCL_LIBRARY_DIRS})
add_definitions(${PCL_DEFINITIONS})

find_package(VTK REQUIRED)
include(${VTK_USE_FILE})

add_executable(narrow_band main.cpp NarrowBand.cpp)
target_link_libraries(narrow_band ${PCL_LIBRARIES} ${VTK_LIBRARIES})
