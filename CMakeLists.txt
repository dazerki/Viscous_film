# Minimum version
cmake_minimum_required(VERSION 3.12.2)

# Project name
project(Project CXX)

# Exec name
set(EXEC_NAME project)
set(CMAKE_CXX_FLAGS "-O3")

# CMake module
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/")


# Python
find_package(PythonLibs REQUIRED)
include_directories(${PYTHON_INCLUDE_DIRS})

# OpenMP
find_package(OpenMP REQUIRED)
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
include_directories(${OPENMP_CXX_INCLUDE_DIRS})


# Deps
include_directories(${PROJECT_SOURCE_DIR}/deps)

# Project files
include_directories(${PROJECT_SOURCE_DIR}/src ${PROJECT_SOURCE_DIR}/inc)
file(GLOB C_SRCS
    "${PROJECT_SOURCE_DIR}/src/*.c"
    "${PROJECT_SOURCE_DIR}/src/*.cpp"
)

# Use G++ for compiling/linking the source files
set_source_files_properties(${C_SRCS} PROPERTIES LANGUAGE CXX)
add_executable(${EXEC_NAME} ${C_SRCS})
target_link_libraries(${EXEC_NAME} ${PYTHON_LIBRARIES} ${OPENMP_CXX_LIBRARIES} pthread m)
