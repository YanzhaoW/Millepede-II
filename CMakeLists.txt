cmake_minimum_required(VERSION 3.25)

project(Millepede2 LANGUAGES Fortran CXX C)

option(SUPPORT_READ_C "Support read CFiles" ON)
option(SUPPORT_C_RFIO "Support read CFiles" OFF)
option(SUPPORT_ZLIB "Support zlib" ON)
option(SUPPORT_OPENMP "Support OpenMP" ON)
option(SUPPORT_LAPACK64 "Support OpenMP" OFF)

set(CMAKE_Fortran_FLAGS
    "${CMAKE_Fortran_FLAGS} -Wall -fautomatic -fno-backslash -O3")
set(CMAKE_C_FLAGS
    "${CMAKE_C_FLAGS} -Wall -Wno-unused-function -O3 -Df2cFortran")
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})

add_subdirectory(src)
