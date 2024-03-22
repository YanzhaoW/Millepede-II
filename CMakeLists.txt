cmake_minimum_required(VERSION 3.25)

project(Millepede2 LANGUAGES Fortran CXX C)
enable_language(Fortran)

set(CMAKE_Fortran_FLAGS
    "${CMAKE_Fortran_FLAGS} -Wall -fautomatic -fno-backslash -O3")
set(CMAKE_C_FLAGS
    "${CMAKE_C_FLAGS} -Wall -Wno-unused-function -O3 -Df2cFortran")

set(MILLEPEDE_modules
    Dbandmatrix.f90
    linesrch.f90
    mille.f90
    minresDataModule.f90
    minresModule.f90
    minresqlpBlasModule.f90
    minresqlpDataModule.f90
    minresqlpModule.f90
    mpbits.f90
    mpdalc.f90
    mpdef.f90
    mphistab.f90
    mpmod.f90
    mpmon.f90
    mpnum.f90
    mpqldec.f90
    mptest1.f90
    mptest2.f90
    mptext.f90
    randoms.f90
    vertpr.f90)
add_library(millepede SHARED ${MILLEPEDE_modules})
target_include_directories(millepede PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

set(MILLESRC Mille.cc)
add_library(mille SHARED ${MILLESRC})
target_include_directories(mille PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

add_executable(pede pede.f90)
target_link_libraries(pede PUBLIC millepede)
