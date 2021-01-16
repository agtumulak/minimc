#include "Material.hpp"

#include <algorithm>
#include <numeric>
#include <random>

// Cell

//// public

pugi::xml_node Material::FindNode(
    const pugi::xml_node& root, const std::string& material_name) {
  const auto& materials_node{root.child("materials")};
  const auto& material_node{
      materials_node.find_child_by_attribute("name", material_name.c_str())};
  if (!material_node) {
    throw std::runtime_error(
        "Material node \"" + material_name + "\" not found. Must be one of: [" +
        std::accumulate(
            materials_node.begin(), materials_node.end(), std::string{},
            [](const auto& accumulated, const auto& material_node) noexcept {
              return accumulated + "\"" +
                     material_node.attribute("name").as_string() + "\", ";
            }) +
        "]");
  }
  return material_node;
}

Material::Material(
    const pugi::xml_node& root, const std::string& name,
    const std::vector<std::shared_ptr<const Nuclide>>& all_nuclides)
    : name{name},
      number_density{FindNode(root, name).attribute("aden").as_double()},
      afracs{AssignNuclides(root, name, all_nuclides)} {}

Material::CrossSection
Material::GetMicroscopicTotal(const Particle& p) const noexcept {
  // TODO: Memoize this call for performance
  return std::accumulate(
      afracs.cbegin(), afracs.cend(), CrossSection{0},
      [&p](const auto& accumulated, const auto& nuclide_afrac_pair) {
        const auto& [nuclide, afrac] = nuclide_afrac_pair;
        return accumulated + afrac * nuclide->GetTotal(p);
      });
}

Real Material::SampleCollisionDistance(
    std::minstd_rand& rng, const Particle& p) const noexcept {
  return std::exponential_distribution{number_density * GetMicroscopicTotal(p)}(
      rng);
}

const Nuclide& Material::SampleNuclide(
    std::minstd_rand& rng, const Particle& p) const noexcept {
  const CrossSection threshold{
      std::uniform_real_distribution{}(rng)*GetMicroscopicTotal(p)};
  CrossSection accumulated{0};
  for (const auto& [nuclide_ptr, afrac] : afracs) {
    accumulated += afrac * nuclide_ptr->GetTotal(p);
    if (accumulated > threshold) {
      return *nuclide_ptr;
    }
  }
  // This should never be reached since Material total cross section is
  // computed from constituent Nuclide total cross sections
  assert(false);
}

//// private

Material::AtomFractionMap Material::AssignNuclides(
    const pugi::xml_node& root, const std::string& name,
    const std::vector<std::shared_ptr<const Nuclide>>& all_nuclides) {
  AtomFractionMap material_atom_fractions;
  const auto& material_node = Material::FindNode(root, name);
  for (const auto& nuclide_node : material_node) {
    const auto nuclide_name = nuclide_node.attribute("name").as_string();
    const auto nuclide_it = std::find_if(
        all_nuclides.cbegin(), all_nuclides.cend(),
        [&nuclide_name](const std::shared_ptr<const Nuclide>& nuclide) {
          return nuclide->name == nuclide_name;
        });
    // World::CreateNuclides should have created all Nuclide objects needed
    assert(nuclide_it != all_nuclides.cend());
    material_atom_fractions.emplace(
        *nuclide_it, nuclide_node.attribute("afrac").as_double());
  }
  // Normalize atom fractions to unity because saying unity makes you sound smart
  const auto sum_afrac = std::accumulate(
      material_atom_fractions.cbegin(), material_atom_fractions.cend(), Real{0},
      [](const auto& accumulated, const auto& nuclide_afrac_pair) {
        return accumulated + nuclide_afrac_pair.second;
      });
  std::for_each(
      material_atom_fractions.begin(), material_atom_fractions.end(),
      [&sum_afrac](auto& nuclide_afrac_pair) {
        nuclide_afrac_pair.second = nuclide_afrac_pair.second / sum_afrac;
      });
  return material_atom_fractions;
}
