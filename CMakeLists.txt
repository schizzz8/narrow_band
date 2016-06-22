cmake_minimum_required (VERSION 2.6)
project (narrow_band)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
set(EXECUTABLE_OUTPUT_PATH /home/dede/level_set_ws/bin/)

find_package(PCL 1.6 REQUIRED)
include_directories(${PCL_INCLUDE_DIRS})
link_directories(${PCL_LIBRARY_DIRS})
add_definitions(${PCL_DEFINITIONS})

find_package(VTK REQUIRED)
include(${VTK_USE_FILE})

aux_source_directory(. SRC_LIST)
add_executable(${PROJECT_NAME} ${SRC_LIST})
target_link_libraries(${PROJECT_NAME} ${PCL_LIBRARIES} ${VTK_LIBRARIES})
