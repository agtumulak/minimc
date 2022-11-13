#include "Driver.hpp"

#include "Cell.hpp"
#include "FixedSource.hpp"
#include "Material.hpp"
#include "Nuclide.hpp"
#include "Particle.hpp"
#include "Reaction.hpp"
#include "Perturbation.hpp"
#include "StreamDelegate.hpp"
#include "XMLDocument.hpp"
#include "pugixml.hpp"

#include <cassert>
#include <iosfwd>
#include <string>

// Driver

//// public

std::unique_ptr<Driver>
Driver::Create(const std::filesystem::path& xml_filepath) {
  // Bind the lifetime of the XMLDocument (which may be large) to the scope of
  // this factory function
  auto doc = std::make_unique<XMLDocument>(xml_filepath);
  const std::string problem_type =
      doc->root.child("problemtype").first_child().name();
  if (problem_type == "fixedsource") {
    return std::make_unique<FixedSource>(doc->root);
  }
  else {
    assert(false); // this should have been caught by the validator
    return {};
  }
}

Driver::Driver(const pugi::xml_node& root)
    : world{root}, perturbations{[this, &root]() {
        std::vector<std::unique_ptr<const Perturbation>> result;
        for (const auto& perturbation_node : root.child("perturbations")) {
          result.push_back(Perturbation::Create(perturbation_node, world));
        }
        return result;
      }()},
      stream_delegate{StreamDelegate::Create(root, world)},
      batchsize(
          std::stoi(root.child("general").child("histories").child_value())),
      seed(std::stoi(
          root.child("general").child("seed")
              ? root.child("general").child("seed").child_value()
              : "1")),
      init_estimator_set{
          root.child("estimators"), world, perturbations,
          static_cast<Real>(batchsize)},
      threads{
          std::stoul(root.child("general").child("threads").child_value())} {}

Driver::~Driver() noexcept {}

void Driver::Transport(
    Particle& p, std::vector<EstimatorProxy>& estimators) noexcept {
  while (p.IsAlive()) {
    // sample the next position
    stream_delegate->StreamToNextCollision(p, estimators, world);
    if (p.reaction == Reaction::leak) {
      break;
    }
    // sample the next Nuclide, currently no need to delegate this
    const auto& sampled_nuclide = [&p]() {
      const MicroscopicCrossSection threshold =
          p.cell->material->GetMicroscopicTotal(p) * p.Sample();
      MicroscopicCrossSection accumulated = 0;
      for (const auto& [nuclide_ptr, afrac] : p.cell->material->afracs) {
        accumulated += afrac * nuclide_ptr->GetTotal(p);
        if (accumulated > threshold) {
          return nuclide_ptr;
        }
      }
      // Material total cross section is computed by adding Nuclide total cross
      // sections so this should never be reached
      assert(false);
    }();
    // interact with the sampled nuclide
    sampled_nuclide->Interact(p, estimators);
  }
}
