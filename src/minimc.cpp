#include "Driver.hpp"

#include <filesystem>
#include <iostream>
#include <memory>
#include <stdexcept>

int main(int argc, char* argv[]) {
  std::cout << "Welcome to MiniMC!" << std::endl;
  if (argc != 2) {
    throw std::runtime_error("MiniMC accepts exactly one argument");
  }
  const std::filesystem::path input_filepath{argv[1]};
  // output filepath is same as input filepath but with `.out` extension
  auto output_filepath = input_filepath;
  output_filepath.replace_extension(".out");
  Driver::Create(input_filepath, output_filepath)->Solve();
  return 0;
}
