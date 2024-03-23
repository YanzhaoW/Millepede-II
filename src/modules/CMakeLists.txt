set(MILLEPEDE_modules
    Dbandmatrix.f90
    linesrch.f90
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

add_library(millepede STATIC ${MILLEPEDE_modules})
target_include_directories(millepede PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
set_target_properties(millepede PROPERTIES Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
