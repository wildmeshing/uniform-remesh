cmake_minimum_required(VERSION 3.16)
project(example)

list(PREPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

# Libigl
option(LIBIGL_WITH_OPENGL            "Use OpenGL"         OFF)
option(LIBIGL_WITH_OPENGL_GLFW       "Use GLFW"           OFF)
option(LIBIGL_WITH_CGAL       "Use CGAL"           ON)
include(libigl)

# Add your project files
file(GLOB SRC_FILES cgal_iso_remesh.cpp)
add_executable(remesh ${SRC_FILES})
target_link_libraries(remesh PUBLIC igl::core igl::cgal)

file(GLOB SRC_FILES cgal_decimate.cpp)
add_executable(sec ${SRC_FILES})
target_link_libraries(sec PUBLIC igl::core igl::cgal)