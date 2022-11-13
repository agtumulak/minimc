#include "Reaction.hpp"

#include <cassert>
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
    // this should have been caught by the validator
    // there should be no code path where `name` is "birth" or "leak"
    assert(false);
    return {};
  };
}
