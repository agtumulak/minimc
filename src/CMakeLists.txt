# https://stackoverflow.com/a/29871891/5101335
find_package(Threads REQUIRED)

# Add new source files here
set(MINIMC_SOURCES
  BasicTypes.hpp
  Bank.cpp
  Bank.hpp
  Bins.cpp
  Bins.hpp
  CSGSurface.cpp
  CSGSurface.hpp
  Cell.cpp
  Cell.hpp
  Driver.cpp
  Driver.hpp
  Continuous.cpp
  Continuous.hpp
  ContinuousMap.hpp
  ContinuousEvaluation.cpp
  ContinuousEvaluation.hpp
  ContinuousReaction.cpp
  ContinuousReaction.hpp
  Estimator.cpp
  Estimator.hpp
  Constants.hpp
  FixedSource.cpp
  FixedSource.hpp
  HDF5DataSet.hpp
  Interaction.cpp
  Interaction.hpp
  KEigenvalue.cpp
  KEigenvalue.hpp
  Material.cpp
  Material.hpp
  Multigroup.cpp
  Multigroup.hpp
  Nuclide.cpp
  Nuclide.hpp
  Particle.cpp
  Particle.hpp
  Perturbation.cpp
  Perturbation.hpp
  Point.cpp
  Point.hpp
  Reaction.cpp
  Reaction.hpp
  Sensitivity.cpp
  Sensitivity.hpp
  ScalarField.cpp
  ScalarField.hpp
  Scorable.cpp
  Scorable.hpp
  Source.cpp
  Source.hpp
  ThermalScattering.cpp
  ThermalScattering.hpp
  TransportMethod.cpp
  TransportMethod.hpp
  World.cpp
  World.hpp
  XMLDocument.cpp
  XMLDocument.hpp
  )

add_library(minimc ${MINIMC_SOURCES})
target_include_directories(minimc PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_compile_options(minimc PUBLIC
  $<$<AND:$<CONFIG:Debug>,$<BOOL:${APPLE}>>:-fprofile-instr-generate -fcoverage-mapping -O0 -Wall -Wextra -Wdocumentation>
  $<$<AND:$<CONFIG:Release>,$<CXX_COMPILER_ID:Clang,GNU>>:-O3>)
target_link_options(minimc PUBLIC
  $<$<AND:$<CONFIG:Debug>,$<BOOL:${APPLE}>>:-fprofile-instr-generate>)
target_link_libraries(minimc PUBLIC xerces pugixml hdf5_cpp Threads::Threads)
set_target_properties(minimc PROPERTIES
  RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}
  )
add_executable(runminimc minimc.cpp)
target_link_libraries(runminimc minimc)
