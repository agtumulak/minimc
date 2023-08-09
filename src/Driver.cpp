#include "Driver.hpp"

#include "Cell.hpp"
#include "Constants.hpp"
#include "Estimator/Estimator.hpp"
#include "Estimator/Proxy.hpp"
#include "Estimator/Visitor.hpp"
#include "FixedSource.hpp"
#include "Material.hpp"
#include "Nuclide.hpp"
#include "Particle.hpp"
#include "Perturbation/Perturbation.hpp"
#include "Point.hpp"
#include "Reaction.hpp"
#include "XMLDocument.hpp"
#include "pugixml.hpp"

#include <cassert>
#include <map>
#include <random>
#include <string>

class CSGSurface;

// Driver

//// public

std::unique_ptr<Driver> Driver::Create(
    const std::filesystem::path& xml_filepath,
    const std::filesystem::path& output_filepath) {
  // Bind the lifetime of the XMLDocument (which may be large) to the scope of
  // this factory function
  auto doc = std::make_unique<XMLDocument>(xml_filepath);
  const std::string problem_type =
      doc->root.child("problemtype").first_child().name();
  if (problem_type == "fixedsource") {
    return std::make_unique<FixedSource>(doc->root, output_filepath);
  }
  else {
    assert(false); // this should have been caught by the validator
    return {};
  }
}

Driver::Driver(
    const pugi::xml_node& root, const std::filesystem::path& output_filepath)
    : world{root}, perturbations{[this, root]() {
        std::vector<std::unique_ptr<const Perturbation::Interface>> result;
        for (const auto& perturbation_node : root.child("perturbations")) {
          result.push_back(
              Perturbation::Interface::Create(perturbation_node, world));
        }
        return result;
      }()},
      estimators{[this, root]() {
        std::vector<std::unique_ptr<Estimator::Interface>> result;
        for (const auto& estimator_node : root.child("estimators")) {
          result.push_back(Estimator::Interface::Create(
              estimator_node, world, perturbations));
        }
        return result;
      }()},
      output_filepath{output_filepath},
      total_weight{
          std::stoull(root.child("general").child("histories").child_value())},
      threads{std::stoul(root.child("general").child("threads").child_value())},
      seed(std::stoi(
          root.child("general").child("seed")
              ? root.child("general").child("seed").child_value()
              : "1")) {}

Driver::~Driver() noexcept {}

void Driver::Transport(
    Particle& p,
    std::vector<Estimator::Proxy>& estimator_proxies) const noexcept {
  // Associate each sensitivity proxy with the corresponding indirect effect
  for (auto& estimator_proxy : estimator_proxies) {
    estimator_proxy.SetIndirectEffects(p);
  }
  // Reassign the indirect effects of each
  // Perturbation::Sensitivity::Proxy::Interface at the beginning of Transport
  // rather than search it each time an Estimator is scored
  while (p.IsAlive()) {
    // sample the next collision point
    Stream(p, estimator_proxies);
    if (p.reaction == Reaction::leak) {
      break;
    }
    // sample the next Nuclide
    const auto& sampled_nuclide = [&p]() {
      const MicroscopicCrossSection threshold =
          p.cell->material->GetCellMajorant(p) * p.Sample();
      MicroscopicCrossSection accumulated = 0;
      for (const auto& [nuclide_ptr, afrac] : p.cell->material->afracs) {
        accumulated += afrac * nuclide_ptr->GetCellMajorant(p);
        if (accumulated > threshold) {
          return nuclide_ptr;
        }
      }
      // Material total cross section is computed by adding Nuclide total cross
      // sections so this should never be reached
      assert(false);
    }();
    // interact with the sampled nuclide
    sampled_nuclide->Interact(p, estimator_proxies);
  }
}

//// private

void Driver::Stream(
    Particle& p,
    std::vector<Estimator::Proxy>& estimator_proxies) const noexcept {
  while (true) {
    const auto distance_to_collision = std::exponential_distribution{
        p.GetCell().material->number_density *
        p.GetCell().material->GetCellMajorant(p)}(p.rng);
    const auto [nearest_surface, distance_to_surface_crossing] =
        p.GetCell().NearestSurface(p.GetPosition(), p.GetDirection());
    // check if collision within Cell has occured
    if (distance_to_collision < distance_to_surface_crossing) {
      // move the Particle to new position
      p.SetPosition(p.GetPosition() + p.GetDirection() * distance_to_collision);
      // this is where I would score track length estimators...if I had one
      // caller handles collision
      return;
    }
    else {
      // collision did not occur so Particle has streamed to adjacent Cell
      p.SetPosition(
          p.GetPosition() +
          p.GetDirection() * (distance_to_collision + constants::nudge));
      const auto& new_cell = world.FindCellContaining(p.GetPosition());
      p.SetCell(new_cell);
      // update Estimator objects
      for (auto& estimator_proxy : estimator_proxies) {
        class Visitor : public Estimator::Visitor {
        public:
          Visitor(const Particle& p, const CSGSurface& s)
              : Estimator::Visitor{p}, s{s} {}
          Score Visit(const Estimator::Current& current_estimator)
              const noexcept final {
            return current_estimator.surface.get() == &s ? 1 : 0;
          }

        private:
          const CSGSurface& s;
        };
        estimator_proxy.Visit(Visitor(p, *nearest_surface));
      }
      // if Particle leaked, update Reaction
      if (new_cell.IsVoid()) {
        p.reaction = Reaction::leak;
        break;
      }
    }
  }
}
