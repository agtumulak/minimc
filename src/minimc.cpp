#include "Driver.hpp"
#include "Estimator.hpp"

#include <filesystem>
#include <fstream>
#include <iostream>
#include <memory>
#include <stdexcept>

int main(int argc, char* argv[]) {
  std::cout << "Welcome to MiniMC!" << std::endl;
  if (argc != 2) {
    throw std::runtime_error("MiniMC accepts exactly one argument");
  }
  const std::filesystem::path input_filepath{argv[1]};
  auto driver = Driver::Create(input_filepath);
  auto output_filepath = input_filepath;
  std::ofstream output_file {output_filepath.replace_extension(".out")};
  const auto& result = driver->Solve();
  output_file << driver->batchsize << std::endl;
  output_file << result.to_string();
  output_file.close();
  std::cout << "Output written to "
            << std::filesystem::absolute(output_filepath) << std::endl;
  return 0;
}
