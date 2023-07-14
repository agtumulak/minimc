cmake_minimum_required(VERSION 3.18)
project(minimc VERSION 0.5.1)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
set(CMAKE_CXX_STANDARD 17)

# Set a default build type if none was specified
# https://blog.kitware.com/cmake-and-the-default-build-type/
set(DEFAULT_BUILD_TYPE "Release")
if(EXISTS "${CMAKE_SOURCE_DIR}/.git")
  set(DEFAULT_BUILD_TYPE "Debug")
endif()
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  message(STATUS "No build type specified. Defaulting to ${DEFAULT_BUILD_TYPE}")
  set(CMAKE_BUILD_TYPE "${DEFAULT_BUILD_TYPE}" CACHE
    STRING "Choose the type of build." FORCE)
  # Set the possible values of build type for cmake-gui
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release")
endif()

# platform-dependent setup
if (MSVC)
  add_compile_definitions(_USE_MATH_DEFINES) # https://stackoverflow.com/a/63810121/5101335
endif()

# add third party libraries
# https://github.com/catchorg/Catch2/blob/devel/docs/cmake-integration.md#cmake-targets
add_subdirectory(third-party/Catch2)
add_subdirectory(third-party/eigen)
add_subdirectory(third-party/pugixml)

include(ExternalProject)

# add xerces as an external project and create a static library
set(XERCES_PREFIX ${CMAKE_BINARY_DIR}/xerces)
set(XERCES_BINARY_DIR ${XERCES_PREFIX}/build)
set(XERCES_INSTALL_DIR ${XERCES_PREFIX}/install)
# logic to determine where static library will be found
if (MSVC)
  if (CMAKE_BUILD_TYPE MATCHES "Debug")
    set(XERCES_SUFFIX "_3D")
  else()
    set(XERCES_SUFFIX "_3")
  endif()
endif()
set(XERCES_LIB ${XERCES_INSTALL_DIR}/lib/${CMAKE_STATIC_LIBRARY_PREFIX}xerces-c${XERCES_SUFFIX}${CMAKE_STATIC_LIBRARY_SUFFIX})
if (UNIX)
  # since xerces is an external library, we do not care about compiler warnings
  set(UNIX_XERCES_CMAKE_ARGS -DCMAKE_CXX_FLAGS="-w")
  if (NOT APPLE)
    set(UNIX_XERCES_CMAKE_ARGS ${UNIX_XERCES_CMAKE_ARGS} -Dtranscoder=gnuiconv)
  endif()
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
  URL https://dlcdn.apache.org//xerces/c/3/sources/xerces-c-3.2.4.tar.xz
  URL_HASH SHA256=075bc57940da0f9be6dd183c550c8ce0b9833e4550dc382048377a1a5e3b2bd9
  CMAKE_ARGS
    -DBUILD_SHARED_LIBS:BOOL=OFF
    -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
    -DCMAKE_INSTALL_PREFIX=${XERCES_INSTALL_DIR}
    ${UNIX_XERCES_CMAKE_ARGS}
  BUILD_BYPRODUCTS ${XERCES_LIB}
  )
add_library(xerces STATIC IMPORTED)
add_dependencies(xerces xerces_download)
set(XERCES_INCLUDE_DIR ${XERCES_INSTALL_DIR}/include)
file(MAKE_DIRECTORY ${XERCES_INCLUDE_DIR})
set_target_properties(xerces PROPERTIES
  IMPORTED_LOCATION ${XERCES_LIB}
  INTERFACE_INCLUDE_DIRECTORIES ${XERCES_INCLUDE_DIR})
if (APPLE)
  target_link_options(xerces INTERFACE -framework CoreServices -lcurl)
endif()

# add hdf5 as an external project and create a static library
set(HDF5_PREFIX ${CMAKE_BINARY_DIR}/hdf5)
set(HDF5_BINARY_DIR ${HDF5_PREFIX}/build)
set(HDF5_INSTALL_DIR ${HDF5_PREFIX}/install)
# logic to determine where static libraries will be found
if (MSVC)
  set(HDF5_PREFIX "lib")
  if (CMAKE_BUILD_TYPE MATCHES "Debug")
    set(HDF5_SUFFIX _D)
  endif()
else()
  set(HDF5_PREFIX ${CMAKE_STATIC_LIBRARY_PREFIX})
  if (CMAKE_BUILD_TYPE MATCHES "Debug")
    set(HDF5_SUFFIX _debug)
  endif()
endif()
set(HDF5_SUFFIX ${HDF5_SUFFIX}${CMAKE_STATIC_LIBRARY_SUFFIX})
set(HDF5_LIB ${HDF5_INSTALL_DIR}/lib/${HDF5_PREFIX}hdf5${HDF5_SUFFIX})
set(HDF5_CPP_LIB ${HDF5_INSTALL_DIR}/lib/${HDF5_PREFIX}hdf5_cpp${HDF5_SUFFIX})
if (UNIX)
  # since hddf5 is an external library, we do not care about compiler warnings
  set(UNIX_HDF5_CMAKE_ARGS -DCMAKE_CXX_FLAGS="-w")
endif()
ExternalProject_Add(
  hdf5_download
  PREFIX ${HDF5_PREFIX}
  SOURCE_DIR ${HDF5_PREFIX}/src
  TMP_DIR ${HDF5_BINARY_DIR}
  STAMP_DIR ${HDF5_BINARY_DIR}
  DOWNLOAD_DIR ${HDF5_BINARY_DIR}
  BINARY_DIR ${HDF5_BINARY_DIR}
  GIT_REPOSITORY https://github.com/HDFGroup/hdf5
  GIT_TAG hdf5-1_8_22
  UPDATE_COMMAND ""
  CMAKE_ARGS
    -DBUILD_SHARED_LIBS:BOOL=OFF
    -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
    -DCMAKE_INSTALL_PREFIX=${HDF5_INSTALL_DIR} # https://stackoverflow.com/a/32287517/5101335
    ${UNIX_XERCES_CMAKE_ARGS}
  BUILD_BYPRODUCTS ${HDF5_LIB} ${HDF5_CPP_LIB}
  )
set(HDF5_INCLUDE_DIR ${HDF5_INSTALL_DIR}/include)
file(MAKE_DIRECTORY ${HDF5_INCLUDE_DIR})

add_library(hdf5 STATIC IMPORTED)
add_dependencies(hdf5 hdf5_download)
set_target_properties(hdf5 PROPERTIES
  IMPORTED_LOCATION ${HDF5_LIB}
  INTERFACE_INCLUDE_DIRECTORIES ${HDF5_INCLUDE_DIR})
add_library(hdf5_cpp STATIC IMPORTED)
add_dependencies(hdf5_cpp hdf5_download)
set_target_properties(hdf5_cpp PROPERTIES
  IMPORTED_LOCATION ${HDF5_CPP_LIB}
  INTERFACE_INCLUDE_DIRECTORIES ${HDF5_INCLUDE_DIR})
target_link_libraries(hdf5_cpp INTERFACE hdf5 ${CMAKE_DL_LIBS}) # https://stackoverflow.com/a/56842079/5101335

add_subdirectory(src)
add_subdirectory(test)
add_subdirectory(doc)
