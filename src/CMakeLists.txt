add_executable(minimc minimc.cpp)
add_dependencies(minimc xerces)
target_compile_options(minimc PRIVATE -Wall)
target_link_libraries(minimc xerceslib)
set_target_properties(minimc PROPERTIES
  RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}
  )
