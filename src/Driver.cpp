#include "Driver.hpp"

#include "FixedSource.hpp"
#include "XMLDocument.hpp"

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
    assert(false); // to be implemented
  }
  else {
    assert(false); // this should have been caught by the validator
  }
}

Driver::Driver(const pugi::xml_node& root)
    : world{root}, batchsize(std::stoi(
                       root.child("general").child("histories").child_value())),
      threads{std::stoul(root.child("general").child("threads").child_value())},
      chunksize{
          std::stoul(root.child("general").child("chunksize").child_value())} {}

Driver::~Driver() noexcept {}