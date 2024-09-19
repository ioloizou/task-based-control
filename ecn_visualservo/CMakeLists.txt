cmake_minimum_required(VERSION 2.8.3)
project(ecn_visualservo)

find_package(VISP REQUIRED)
find_package(log2plot REQUIRED)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_BUILD_TYPE Debug)

add_definitions(-DBASE_PATH="${CMAKE_CURRENT_SOURCE_DIR}/")

include_directories(
  include
  ${VISP_INCLUDE_DIRS}
  ${LOG2PLOT_INCLUDE_DIRS}
)

add_library(${PROJECT_NAME} SHARED include/simulator.h src/simulator.cpp
  include/feature_stack.h src/feature_stack.cpp)
target_link_libraries(${PROJECT_NAME} ${VISP_LIBRARIES} ${LOG2PLOT_LIBRARIES} )


add_executable(compare_features main.cpp)
target_link_libraries(compare_features ${PROJECT_NAME})

