#include "Driver.hpp"
#include "Estimator.hpp"

#include <iostream>
#include <memory>
#include <stdexcept>


#include "pybind11/pybind11.h"
#include "pybind11/embed.h"

namespace py = pybind11;

PYBIND11_EMBEDDED_MODULE(fast_calc, m){
  m.def("add", [](int i, int j){
    return i + j;
  });
}

int main() {
  py::scoped_interpreter guard{};
  auto fast_calc = py::module_::import("fast_calc");
  auto result = fast_calc.attr("add")(1,2).cast<int>();
  py::print(result);
  // py::exec("import sys; sys.path.append('/Users/atumulak/.local/share/virtualenvs/pyminimc-tdgjIrhJ/lib/python3.10/site-packages/'); import numpy as np; print(np.array([1,2,3]))");
   py::exec("import sys; sys.path.append('/Users/atumulak/.local/share/virtualenvs/pyminimc-tdgjIrhJ/lib/python3.10/site-packages/'); import numpy as np; print(np.array([1,2,3]))");
  return 0;
}
