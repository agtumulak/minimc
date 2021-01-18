#include "FixedSource.hpp"
#include "XMLDocument.hpp"

#include <iostream>
#include <stdexcept>

int main(int argc, char* argv[]) {
  std::cout << "Welcome to MiniMC!" << std::endl;
  if (argc != 2){
    throw std::runtime_error("MiniMC accepts exactly one argument");
  }
  XMLDocument doc{argv[1]};
  FixedSourceStandalone(doc.root).Solve();
  return 0;
}
