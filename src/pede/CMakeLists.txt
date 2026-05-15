add_executable(pede pede.f90)
target_link_libraries(pede PRIVATE millepede_core mille_c)
if(SUPPORT_READ_C)
    message(STATUS "READ_C is supported")
    target_compile_definitions(pede PUBLIC -DREAD_C_FILES)
endif()

install(TARGETS pede EXPORT millepede)
