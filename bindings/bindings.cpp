#include <vector>

#include "../include/pybind11/stl.h"
#include "../include/pybind11/pybind11.h"
#include "../include/pybind11/eigen.h"

#include "../lgmres.h"

namespace py = pybind11;


PYBIND11_PLUGIN(LGMRES)
{
    py::module m("LGMRES",
                 "Python bindings to the LGMRES C++ implementation");

    m.def("solve",
          &lgmres_wrapper,
          py::arg("A"),
          py::arg("b"),
          py::arg("x0")=Vec(0),
          py::arg("M")=Mat(0,0),
          py::arg("tol")=1e-5,
          py::arg("maxiter")=1000,
          py::arg("inner_m")=30,
          py::arg("outer_k")=3,
          py::arg("outer_v")=std::vector<Vec>());

    return m.ptr();
}
