cmake_minimum_required(VERSION 3.25)

project(Millepede2 LANGUAGES Fortran CXX C)

include(GNUInstallDirs)

set(LIBRARY_OUTPUT_PATH "${CMAKE_BINARY_DIR}/lib")
set(EXECUTABLE_OUTPUT_PATH "${CMAKE_BINARY_DIR}/bin")

option(SUPPORT_READ_C "Support read CFiles" ON)
option(SUPPORT_ZLIB "Support zlib" ON)
option(SUPPORT_OPENMP "Support OpenMP" ON)
option(SUPPORT_LAPACK64 "Support LAPACK64" OFF)

set(CMAKE_Fortran_FLAGS
    "${CMAKE_Fortran_FLAGS} -Wall -fautomatic -fno-backslash -O3 -cpp")
set(CMAKE_C_FLAGS
    "${CMAKE_C_FLAGS} -Wall -Wno-unused-function -O3 -Df2cFortran")
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

if(SUPPORT_ZLIB)
  message(STATUS "Zlib is supported")
  find_package(ZLIB REQUIRED)
endif()

if(SUPPORT_OPENMP)
  message(STATUS "OPENMP is supported")
  find_package(OpenMP REQUIRED)
endif()

if(SUPPORT_LAPACK64)
  set(BLA_VENDOR OPENBLAS)
  message(STATUS "LAPACK64 is supported using vendor: ${BLA_VENDOR}")
  find_package(BLAS REQUIRED)
endif()

add_subdirectory(src)
