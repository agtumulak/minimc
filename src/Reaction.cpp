#include "Reaction.hpp"

#include <cassert>

Reaction ToReaction(const std::string& name) noexcept {
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
    assert(false); // only mutually exclusive reactions are valid
  };
}
