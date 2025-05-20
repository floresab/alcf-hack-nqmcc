cmake_minimum_required(VERSION 3.20)
include(FindPackageHandleStandardArgs)
# ---------------------------------------------------------------------
project(alcf-hack-2025 VERSION 0.9
        DESCRIPTION ""
        LANGUAGES Fortran)
# ---------------------------------------------------------------------
enable_language(Fortran)
# ---------------------------------------------------------------------
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY bin)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY lib)
set(CMAKE_INSTALL_PREFIX ".") 
set(CMAKE_Fortran_MODULE_DIRECTORY ${PROJECT_BINARY_DIR}/mod)
# ---------------------------------------------------------------------
add_library(alcf_hack  init.f90)
add_executable(qmc nqmcc.f90)
target_link_libraries(qmc PUBLIC alcf_hack)
# ---------------------------------------------------------------------
