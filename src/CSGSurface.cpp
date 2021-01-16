#include "CSGSurface.hpp"

#include <cmath>
#include <limits>
#include <numeric>
#include <stdexcept>

// CSGSurface

//// public

std::unique_ptr<const CSGSurface>
CSGSurface::Create(const pugi::xml_node& root, const std::string& name) {
  const auto& surfaces_node{root.child("surfaces")};
  const auto& surface_node{
      surfaces_node.find_child_by_attribute("name", name.c_str())};
  if (!surface_node) {
    throw std::runtime_error(
        "Surface node \"" + name + "\" not found. Must be one of: [" +
        std::accumulate(
            surfaces_node.begin(), surfaces_node.end(), std::string{},
            [](const auto& accumulated, const auto& surface_node) noexcept {
              return accumulated + "\"" +
                     surface_node.attribute("name").as_string() + "\", ";
            }) +
        "]");
  }
  // construct appropriate Surface object
  std::unique_ptr<const CSGSurface> surface{}; // hope this triggers NRVO
  switch (ToSurfaceType(surface_node.name())) {
  case CSGSurface::SurfaceType::Sphere:
    surface = std::make_unique<const Sphere>(surface_node);
    break;
  }
  return surface;
}

CSGSurface::~CSGSurface() noexcept {}

//// protected

CSGSurface::CSGSurface(const pugi::xml_node& surface_node) noexcept
    : name{surface_node.attribute("name").as_string()} {}

//// private

CSGSurface::SurfaceType
CSGSurface::ToSurfaceType(const std::string& surface_name) noexcept {
  if (surface_name == "sphere") {
    return CSGSurface::SurfaceType::Sphere;
  }
  assert(false); // this should have been caught by the validator
}

Real CSGSurface::SolveQuadratic(Real a, Real b, Real c) noexcept {
  const Real discriminant{b * b - 4 * a * c};
  if (discriminant <= 0) {
    // discriminant < 0: no intersection
    // discriminant == 0: grazing, considered no intersection
    return std::numeric_limits<Real>::infinity();
  }
  const Real lesser{(-b - std::sqrtf(discriminant)) / (2 * a)};
  const Real greater{(-b + std::sqrtf(discriminant)) / (2 * a)};
  if (lesser > 0) {
    // Outside sphere; headed towards sphere
    return lesser;
  }
  else if (greater > 0) {
    // On sphere or inside sphere; leaving sphere
    return greater;
  }
  else {
    // On sphere or outside sphere; headed away from sphere
    return std::numeric_limits<Real>::infinity();
  }
}

// Sphere

//// public

Sphere::Sphere(const pugi::xml_node& sphere_node) noexcept
    : CSGSurface{sphere_node}, center{sphere_node.child("center")},
      radius{sphere_node.child("radius").attribute("r").as_double()} {}

Real Sphere::Distance(
    const Point& origin, const Point& direction) const noexcept {
  const Point oc{origin - center};
  return SolveQuadratic(1, 2 * (oc * direction), oc * oc - radius * radius);
}

bool Sphere::Contains(const Point& p) const noexcept {
  const Point pc{p - center};
  return pc * pc < radius * radius;
}
