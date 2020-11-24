#include "xercesc/util/PlatformUtils.hpp"

#include <iostream>

int main(int argc, char* argv[]) {
  std::cout << "Welcome to MiniMC!" << std::endl;
  xercesc::XMLPlatformUtils::Initialize();
  xercesc::XMLPlatformUtils::Terminate();
  return 0;
}
