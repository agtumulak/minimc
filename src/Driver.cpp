#include "Driver.hpp"

#include "FixedSource.hpp"
#include "KEigenvalue.hpp"
#include "Particle.hpp"
#include "TransportMethod.hpp"
#include "XMLDocument.hpp"
#include "pugixml.hpp"

#include <cassert>
#include <iosfwd>
#include <string>

// Driver

//// public

std::unique_ptr<Driver>
Driver::Create(const std::filesystem::path& xml_filepath) {
  // limit the lifetime of the XMLDocument to input parsing time
  auto doc = std::make_unique<XMLDocument>(xml_filepath);
  const std::string problem_type =
      doc->root.child("problemtype").first_child().name();
  if (problem_type == "fixedsource") {
    return std::make_unique<FixedSource>(doc->root);
  }
  else if (problem_type == "keigenvalue") {
    return std::make_unique<KEigenvalue>(doc->root);
  }
  else {
    assert(false); // this should have been caught by the validator
    return {};
  }
}

Driver::Driver(const pugi::xml_node& root)
    : world{root}, perturbations{root.child("perturbations"), world},
      batchsize(
          std::stoi(root.child("general").child("histories").child_value())),
      init_estimator_set{
          root.child("estimators"), world, perturbations,
          static_cast<Real>(batchsize)},
      threads{std::stoul(root.child("general").child("threads").child_value())},
      seed(std::stoi(
          root.child("general").child("seed")
              ? root.child("general").child("seed").child_value()
              : "1")) {
  // All Particle objects will use the selected TransportMethod
  Particle::transport_method = TransportMethod::Create(root, world);
}

Driver::~Driver() noexcept {}
