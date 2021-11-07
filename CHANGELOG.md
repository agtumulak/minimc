# Changelog
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]
- Add isotropic flux source direction
- Fixed bug in Bins error parsing

## [0.0.4] - 2021-10-21
- Implemented free gas scattering
- Implemented global temperature specification
- Implemented estimators with user-specified binning structures
- Fixed bug in calculation of scattering cosine in thermal scattering

## [0.0.3] - 2021-09-09
- Moved RNG to be member of Particle
- Implemented HDF5 data support
- Implemented thermal scattering physics
- Add scripts for parsing, processing, and generating S(a,b) data using
  functional expansions in temperature
- Added improved Doxygen stylesheet

## [0.0.2] - 2021-05-13
- Implemented k-eigenvalue criticality solver
- Implemented fission physics and data
- Added continuous energy physics support

## [0.0.1] - 2021-01-16
### Added
- Basic World construction: cells, surfaces, materials, nuclides, and nuclear
  data are parsed from XML input files
- XML input files are validated with Xerces-C++ and loaded into memory with
  pugixml
- Implemented multigroup fixed-source runs with multithreading and simple
  estimators
- Unit tests are executed with Catch2
- Code is documented with Doxygen
- Continuous integration is performed with GitHub Actions

[Unreleased]: https://github.com/agtumulak/minimc/compare/v0.0.1...develop
[0.0.1]: https://github.com/agtumulak/minimc/releases/tag/v0.0.1
[0.0.2]: https://github.com/agtumulak/minimc/releases/tag/v0.0.2
[0.0.3]: https://github.com/agtumulak/minimc/releases/tag/v0.0.3
[0.0.4]: https://github.com/agtumulak/minimc/releases/tag/v0.0.4
