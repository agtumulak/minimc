find_package(Doxygen)
if(Doxygen_FOUND)
  set(DOXYGEN_SOURCES
    mainpage.dox
    developer_guide.dox
    style_guide.dox
    )
  file(GLOB_RECURSE MD_FILES "${PROJECT_SOURCE_DIR}/*.md")
  set(DOXYGEN_USE_MATHJAX YES)
  set(DOXYGEN_MATHJAX_EXTENSIONS "TeX/AMSmath")
  set(DOXYGEN_QUIET YES)
  set(DOXYGEN_EXCLUDE_PATTERNS "*.cpp")
  doxygen_add_docs(doc ${DOXYGEN_SOURCES} ${CMAKE_SOURCE_DIR}/src)
endif()
