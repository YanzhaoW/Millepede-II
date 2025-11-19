configure_file(${CMAKE_CURRENT_SOURCE_DIR}/readMilleBinary.py
               ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/readMilleBinary.py COPYONLY USE_SOURCE_PERMISSIONS)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/compareResults.py
               ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/compareResults.py COPYONLY USE_SOURCE_PERMISSIONS)
