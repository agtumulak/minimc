find_package(Doxygen)
if(Doxygen_FOUND)
  set(DOXYGEN_SOURCES
    mainpage.dox
    )
  file(GLOB_RECURSE MD_FILES "${PROJECT_SOURCE_DIR}/*.md")
  set(DOXYGEN_USE_MATHJAX YES)
  set(DOXYGEN_QUIET YES)
  doxygen_add_docs(doc ${DOXYGEN_SOURCES} ${CMAKE_SOURCE_DIR}/src)
endif()
