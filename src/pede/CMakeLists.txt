add_executable(pede pede.f90)
target_link_libraries(pede PUBLIC millepede mille)
target_compile_options(pede PUBLIC -DREAD_C_FILES)
