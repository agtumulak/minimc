#include "Driver.hpp"

#include <iostream>
#include <stdexcept>

int main(int argc, char* argv[]) {
  std::cout << "Welcome to MiniMC!" << std::endl;
  if (argc != 2) {
    throw std::runtime_error("MiniMC accepts exactly one argument");
  }
  auto driver = Driver::Create(argv[1]);
  driver->Solve();
  return 0;
}
