//
// Created by Joseph Thompson on 2025-05-28.
//Bindings for python

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include "ALBP.h"

namespace py = pybind11;

PYBIND11_MODULE(albp_module, m) {
    py::class_<PrecedenceRelation>(m, "PrecedenceRelation")
        .def(py::init<>())
        .def_readwrite("pred", &PrecedenceRelation::pred)
        .def_readwrite("succ", &PrecedenceRelation::succ);

    py::class_<ALBP>(m, "ALBP")
        .def(py::init<>())  // default
        .def(py::init<const std::string&>())  // from filename
        .def(py::init<int, int, const std::vector<int>&, const std::vector<std::vector<int>>&>(),  // custom constructor
             py::arg("C"),
             py::arg("N"),
             py::arg("task_times"),
             py::arg("raw_precedence"))
        .def("print", &ALBP::print)
        .def("loadFromFile", &ALBP::loadFromFile)
        .def_readwrite("name", &ALBP::name)
        .def_readwrite("C", &ALBP::C)
        .def_readwrite("N", &ALBP::N)
        .def_readwrite("S", &ALBP::S)
        .def_readwrite("task_time", &ALBP::task_time)
        .def_readwrite("prec_mat", &ALBP::prec_mat)
        .def_readwrite("t_close_mat", &ALBP::t_close_mat)
        .def_readwrite("dir_suc", &ALBP::dir_suc)
        .def_readwrite("dir_pred", &ALBP::dir_pred)
        .def_readwrite("suc", &ALBP::suc)
        .def_readwrite("pred", &ALBP::pred)
        .def_readwrite("precedence_relations", &ALBP::precedence_relations)
        .def_readwrite("task_assignment", &ALBP::task_assignment);
}
