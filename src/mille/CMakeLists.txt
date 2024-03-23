set(MILLESRC Mille.cc readc.c mille.f90)
add_library(mille STATIC ${MILLESRC})
target_link_libraries(
  mille
  PUBLIC millepede
  PRIVATE ZLIB::ZLIB)
target_include_directories(mille PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

if(ZLIB_FOUND)
  target_compile_definitions(mille PUBLIC -DUSE_ZLIB)
endif()
