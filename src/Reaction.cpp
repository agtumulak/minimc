#include "Reaction.hpp"

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
    assert(false); // this should have been caught by the validator
    return {};
  };
}
