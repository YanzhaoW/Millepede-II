add_library(mille SHARED)
target_sources(
    mille
    PUBLIC
        FILE_SET public_headers
            TYPE HEADERS
            BASE_DIRS ${CMAKE_CURRENT_SOURCE_DIR}
            FILES
                MilleEntry.hpp
                MilleWriter.hpp
                ParResultReader.hpp
                SteerWriter.hpp
    PRIVATE MilleWriter.cpp ParResultReader.cpp SteerWriter.cpp
)

add_library(mille_c STATIC)
target_sources(
    mille_c
    PUBLIC
        FILE_SET internal_headers
            TYPE HEADERS
            BASE_DIRS ${CMAKE_CURRENT_SOURCE_DIR}
    PRIVATE Mille.cc readc.c mille.f90
)
target_link_libraries(mille_c PUBLIC millepede PRIVATE ZLIB::ZLIB)
if(ZLIB_FOUND)
    target_compile_definitions(mille_c PUBLIC -DUSE_ZLIB)
endif()

add_library(millepede::mille ALIAS mille)

install(TARGETS mille EXPORT millepede FILE_SET public_headers)
