#include "Driver.hpp"

// Driver

//// public

Driver::Driver(const pugi::xml_node& root)
    : world{root}, batchsize(std::stoi(
                       root.child("general").child("histories").child_value())),
      threads{
          std::stoul(root.child("general").child("histories").child_value())},
      chunksize{
          std::stoul(root.child("general").child("chunksize").child_value())} {}
