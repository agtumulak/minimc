#include "Reaction.hpp"

#include <stdexcept>
#include <string>

Reaction ToReaction(const std::string& name) {
  if (name == "capture") {
    return Reaction::capture;
  }
  else if (name == "scatter") {
    return Reaction::scatter;
  }
  else if (name == "fission") {
    return Reaction::fission;
  }
  else {
    throw std::runtime_error("Unrecognized reaction name: " + name);
  };
}
