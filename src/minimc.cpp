#include "Driver.hpp"
#include "Estimator.hpp"

#include <iostream>
#include <memory>
#include <stdexcept>

int main(int argc, char* argv[]) {
  std::cout << "Welcome to MiniMC!" << std::endl;
  if (argc != 2) {
    throw std::runtime_error("MiniMC accepts exactly one argument");
  }
  auto driver = Driver::Create(argv[1]);
  const auto& result = driver->Solve();
  std::cout << result.to_string();
  return 0;
}
