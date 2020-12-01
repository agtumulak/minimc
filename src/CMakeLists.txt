# Add new source files here
set(MINIMC_SOURCES
  hello_world.cpp
  )

add_library(minimc ${MINIMC_SOURCES})
add_dependencies(minimc xerces_download)
target_include_directories(minimc PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_compile_options(minimc PRIVATE -Wall)
target_link_libraries(minimc xerces)
set_target_properties(minimc PROPERTIES
  RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}
  )
