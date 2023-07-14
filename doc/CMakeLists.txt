find_package(Doxygen)
if(Doxygen_FOUND)
  set(DOXYGEN_SOURCES
    mainpage.dox
    developer_guide.dox
    style_guide.dox
    stream_delegates.dox
    estimators.dox
    )
  file(GLOB_RECURSE MD_FILES "${PROJECT_SOURCE_DIR}/*.md")
  set(DOXYGEN_USE_MATHJAX YES)
  set(DOXYGEN_MATHJAX_EXTENSIONS "TeX/AMSmath")
  set(DOXYGEN_MATHJAX_CODEFILE "mycode.js")
  set(DOXYGEN_QUIET YES)
  set(DOXYGEN_EXCLUDE_PATTERNS "*.cpp")
  set(DOXYGEN_GENERATE_TREEVIEW YES)
  set(DOXYGEN_HTML_HEADER "${CMAKE_CURRENT_SOURCE_DIR}/header.html")
  set(DOXYGEN_HTML_EXTRA_FILES
    ${CMAKE_CURRENT_SOURCE_DIR}/doxygen-awesome-css/doxygen-awesome-darkmode-toggle.js
    ${CMAKE_CURRENT_SOURCE_DIR}/doxygen-awesome-css/doxygen-awesome-fragment-copy-button.js
    ${CMAKE_CURRENT_SOURCE_DIR}/doxygen-awesome-css/doxygen-awesome-paragraph-link.js
    ${CMAKE_CURRENT_SOURCE_DIR}/doxygen-awesome-css/doxygen-awesome-interactive-toc.js)
  set(DOXYGEN_HTML_EXTRA_STYLESHEET
    ${CMAKE_CURRENT_SOURCE_DIR}/doxygen-awesome-css/doxygen-awesome.css
    ${CMAKE_CURRENT_SOURCE_DIR}/doxygen-awesome-css/doxygen-awesome-sidebar-only.css
    ${CMAKE_CURRENT_SOURCE_DIR}/doxygen-awesome-css/doxygen-awesome-sidebar-only-darkmode-toggle.css)
  set(DOXYGEN_IMAGE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/resources)
  doxygen_add_docs(minimcdoc ${DOXYGEN_SOURCES} ${CMAKE_SOURCE_DIR}/src)
endif()
