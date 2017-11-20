cmake_minimum_required(VERSION 2.8)
project(green_cage_deformer)

set(PROJECT_PATH ${CMAKE_CURRENT_BINARY_DIR}/${PROJECT})
set(CMAKE_INSTALL_PREFIX ${PROJECT_PATH})
set(CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

file(GLOB SOURCES src/*.cpp src/*.cc src/*.h src/*.hh)
include_directories(libs)
add_library(${PROJECT_NAME} SHARED ${SOURCES})

find_package(Maya REQUIRED)
find_package(OpenGL)

SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${MAYA_CXX_FLAGS} ${COMPILE_FLAGS}")

include_directories(${MAYA_INCLUDE_DIR})
link_directories(${MAYA_LIBRARY_DIR})
target_link_libraries(${PROJECT_NAME} ${MAYA_LIBRARIES} ${OPENGL_LIBRARIES})

MAYA_PLUGIN(${PROJECT_NAME})
install(TARGETS ${PROJECT_NAME} ${MAYA_TARGET_TYPE} DESTINATION plug-ins)
