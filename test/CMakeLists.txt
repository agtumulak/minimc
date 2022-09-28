# https://github.com/catchorg/Catch2/blob/devel/docs/cmake-integration.md#cmake-targets
Include(FetchContent)
FetchContent_Declare(
  catch_download
  GIT_REPOSITORY https://github.com/catchorg/Catch2.git
  GIT_TAG v3.0.1
  )
FetchContent_MakeAvailable(catch_download)

# add new unit test source files here
set(UNIT_TEST_SOURCES
  test_CSGSurface.cpp
  test_Cell.cpp
  test_FixedSource.cpp
  test_HDF5.cpp
  test_Material.cpp
  test_Nuclide.cpp
  test_Particle.cpp
  test_Point.cpp
  test_ThermalScattering.cpp
  test_XMLDocument.cpp
  test_World.cpp
  )
add_executable(all_tests ${UNIT_TEST_SOURCES})
target_link_libraries(all_tests PRIVATE Catch2::Catch2WithMain minimc)

# Add unit test input files here
# Prefixes are stripped so duplicate filenames must not appear
set(UNIT_TEST_INPUTS
  data/endfb8/n-H1_0-total.hdf5
  data/endfb8/n-H1_0-gamma.hdf5
  data/endfb8/n-H1_0-elastic.hdf5
  data/endfb8/n-O16_0-total.hdf5
  data/endfb8/n-O16_0-gamma.hdf5
  data/endfb8/n-O16_0-elastic.hdf5
  data/endfb8/n-U235_0-total.hdf5
  data/endfb8/n-U235_0-gamma.hdf5
  data/endfb8/n-U235_0-elastic.hdf5
  data/endfb8/n-U235_0-fission.hdf5
  data/endfb8/n-U235_0-nubar.hdf5
  data/tsl/endfb8/majorant_endfb8_safety_factor_1.05.hdf5
  data/tsl/endfb8/scatter_xs_endfb8_T_coeffs.hdf5
  data/tsl/endfb8/scatter_xs_endfb8_S_coeffs.hdf5
  data/tsl/endfb8/scatter_xs_endfb8_E_coeffs.hdf5
  data/tsl/endfb8/beta_endfb8_0_CDF_coeffs.hdf5
  data/tsl/endfb8/beta_endfb8_1_CDF_coeffs.hdf5
  data/tsl/endfb8/beta_endfb8_2_CDF_coeffs.hdf5
  data/tsl/endfb8/beta_endfb8_3_CDF_coeffs.hdf5
  data/tsl/endfb8/beta_endfb8_0_S_coeffs.hdf5
  data/tsl/endfb8/beta_endfb8_1_S_coeffs.hdf5
  data/tsl/endfb8/beta_endfb8_2_S_coeffs.hdf5
  data/tsl/endfb8/beta_endfb8_3_S_coeffs.hdf5
  data/tsl/endfb8/beta_endfb8_0_E_T_coeffs.hdf5
  data/tsl/endfb8/beta_endfb8_1_E_T_coeffs.hdf5
  data/tsl/endfb8/beta_endfb8_2_E_T_coeffs.hdf5
  data/tsl/endfb8/beta_endfb8_3_E_T_coeffs.hdf5
  data/tsl/endfb8/alpha_endfb8_0_CDF_coeffs.hdf5
  data/tsl/endfb8/alpha_endfb8_1_CDF_coeffs.hdf5
  data/tsl/endfb8/alpha_endfb8_2_CDF_coeffs.hdf5
  data/tsl/endfb8/alpha_endfb8_3_CDF_coeffs.hdf5
  data/tsl/endfb8/alpha_endfb8_0_S_coeffs.hdf5
  data/tsl/endfb8/alpha_endfb8_1_S_coeffs.hdf5
  data/tsl/endfb8/alpha_endfb8_2_S_coeffs.hdf5
  data/tsl/endfb8/alpha_endfb8_3_S_coeffs.hdf5
  data/tsl/endfb8/alpha_endfb8_0_beta_T_coeffs.hdf5
  data/tsl/endfb8/alpha_endfb8_1_beta_T_coeffs.hdf5
  data/tsl/endfb8/alpha_endfb8_2_beta_T_coeffs.hdf5
  data/tsl/endfb8/alpha_endfb8_3_beta_T_coeffs.hdf5
  test/continuous.xml
  test/continuous_invalid_nuclide_data_path.xml
  test/minimc.xsd
  test/minimc_basictypes.xsd
  test/multigroup.xml
  test/multigroup_malformed_capture.xml
  test/multigroup_malformed_scatter.xml
  test/multigroup_missing_particle_data.xml
  test/multigroup_nonexistent_nuclide.xml
  test/point_source_leakage.xml
  test/test_HDF5_valid.hdf5
  test/test_HDF5_invalid.hdf5
  test/test_ThermalScattering.xml
  test/test_XMLDocument.xsd
  test/test_XMLDocument_valid.xml
  test/test_XMLDocument_invalid.xml
  )

list(TRANSFORM UNIT_TEST_INPUTS PREPEND ${CMAKE_SOURCE_DIR}/)
add_custom_target(
  copy_xml_inputs ALL
  COMMAND ${CMAKE_COMMAND} -E copy_if_different ${UNIT_TEST_INPUTS} ${CMAKE_CURRENT_BINARY_DIR}
  )

# replace CTest with Catch2
add_custom_target(test
  COMMAND all_tests
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
  COMMENT "Run Catch2 instead of CTest")
add_dependencies(test all_tests copy_xml_inputs)

# code coverage with llvm
if (CMAKE_BUILD_TYPE MATCHES "Debug")
  add_custom_target(coverage
    COMMAND make test
    COMMAND xcrun llvm-profdata merge -sparse default.profraw -o default.profdata
    COMMAND xcrun llvm-cov show -instr-profile default.profdata all_tests -format=html -o report
    COMMAND open report/index.html
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
    COMMENT "Create coverage report on macOS")
  add_dependencies(coverage all_tests)
endif()
