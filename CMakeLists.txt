cmake_minimum_required(VERSION 3.18)
project(minimc VERSION 0.1)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
set(CMAKE_CXX_STANDARD 17)

include(ExternalProject)
include(CTest)

# add xerces as an external project and create a static library
set(XERCES_PREFIX ${CMAKE_BINARY_DIR}/xerces)
set(XERCES_BINARY_DIR ${XERCES_PREFIX}/build)
set(XERCES_INSTALL_DIR ${XERCES_PREFIX}/install)
# logic to determine where static library will be found
if (WIN32)
  set(XERCES_SUFFIX "_3")
  if (CMAKE_BUILD_TYPE STREQUAL "Debug")
    set(XERCES_SUFFIX "${XERCES_SUFFIX}D")
  endif()
endif()
set(XERCES_LIB ${XERCES_INSTALL_DIR}/lib/${CMAKE_STATIC_LIBRARY_PREFIX}xerces-c${XERCES_SUFFIX}${CMAKE_STATIC_LIBRARY_SUFFIX})
if (UNIX)
  # since xerces is an external library, we do not care about compiler warnings
  set(UNIX_XERCES_CMAKE_ARGS -DCMAKE_CXX_FLAGS="-w")
endif()
ExternalProject_Add(
  xerces_download
  PREFIX ${XERCES_PREFIX}
  SOURCE_DIR ${XERCES_PREFIX}/src
  TMP_DIR ${XERCES_BINARY_DIR}
  STAMP_DIR ${XERCES_BINARY_DIR}
  DOWNLOAD_DIR ${XERCES_BINARY_DIR}
  BINARY_DIR ${XERCES_BINARY_DIR}
  INSTALL_DIR ${XERCES_INSTALL_DIR}
  URL https://apache.claz.org//xerces/c/3/sources/xerces-c-3.2.3.tar.xz
  URL_HASH SHA256=12fc99a9fc1d1a79bd0e927b8b5637a576d6656f45b0d5e70ee3694d379cc149
  CMAKE_ARGS
    -DBUILD_SHARED_LIBS:BOOL=OFF
    -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
    -DCMAKE_INSTALL_PREFIX=${XERCES_INSTALL_DIR}
    ${UNIX_XERCES_CMAKE_ARGS}
  BUILD_BYPRODUCTS ${XERCES_LIB}
  )
add_library(xerces STATIC IMPORTED)
set(XERCES_INCLUDE_DIR ${XERCES_INSTALL_DIR}/include)
file(MAKE_DIRECTORY ${XERCES_INCLUDE_DIR})
set_target_properties(xerces PROPERTIES
  IMPORTED_LOCATION ${XERCES_LIB}
  INTERFACE_INCLUDE_DIRECTORIES ${XERCES_INCLUDE_DIR})
if (APPLE)
  target_link_options(xerces INTERFACE -framework CoreServices -lcurl)
endif()

add_subdirectory(src)
add_subdirectory(test)
add_subdirectory(doc)
