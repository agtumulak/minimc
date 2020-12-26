# Add new source files here
set(MINIMC_SOURCES
  BasicTypes.hpp
  CSGSurface.cpp
  CSGSurface.hpp
  Point.cpp
  Point.hpp
  XMLDocument.cpp
  XMLDocument.hpp
  )

add_library(minimc ${MINIMC_SOURCES})
add_dependencies(minimc xerces_download pugixml_download)
target_include_directories(minimc PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_compile_options(minimc PRIVATE -Wall
  $<$<AND:$<CONFIG:Debug>,$<BOOL:${APPLE}>>:-fprofile-instr-generate -fcoverage-mapping>)
target_link_options(minimc PUBLIC
  $<$<AND:$<CONFIG:Debug>,$<BOOL:${APPLE}>>:-fprofile-instr-generate>)
target_link_libraries(minimc xerces pugixml)
set_target_properties(minimc PROPERTIES
  RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}
  )
