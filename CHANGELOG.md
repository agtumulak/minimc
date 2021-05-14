# Changelog
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

## [0.0.2] - 2021-05-13
- Implemented k-eigenvalue criticality solver
- Implemented fission physics and data
- Add continuous energy physics support

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
