//
// Created by Joseph Thompson on 2025-05-28.
//Bindings for python

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include "ALBP.h"
#include "ils.h"
#include "albp_solution.h"
#include "mhh.h"
#include "vdls.h"
#include "salbp_basics.h"
namespace py = pybind11;

PYBIND11_MODULE(ILS_ALBP, m) {
    py::class_<PrecedenceRelation>(m, "PrecedenceRelation")
            .def(py::init<>())
            .def_readwrite("parent", &PrecedenceRelation::parent)
            .def_readwrite("child", &PrecedenceRelation::child);

    py::class_<ALBP>(m, "ALBP")
            .def(py::init<>()) // default
            .def(py::init<const std::string &>()) // from filename
                // Bind static factory methods as static methods
            .def_static("type_1", &ALBP::type_1,
                py::arg("C"), py::arg("N"), py::arg("task_times"), py::arg("raw_precedence"),
                "Factory constructor for type_1")

            .def_static("type_2", &ALBP::type_2,
                py::arg("S"), py::arg("N"), py::arg("task_times"), py::arg("raw_precedence"),
                "Factory constructor for type_2")

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


    // Bind the ALBPSolution class
    py::class_<ALBPSolution>(m, "ALBPSolution")
            // Constructor
            .def(py::init<int>(), "Create ALBPSolution with number of tasks", py::arg("n_tasks"))

            // Public attributes (read/write)
            .def_readwrite("task_assignment", &ALBPSolution::task_assignment,
                           "Task assignment solution vector")
            .def_readwrite("station_assignments", &ALBPSolution::station_assignments,
                           "Station assignments solution (vector of vectors)")
            .def_readwrite("ranking", &ALBPSolution::ranking,
                           "Tasks in order of ranking")
            .def_readwrite("task_ranking", &ALBPSolution::task_ranking,
                           "Ranking for each task")
            .def_readwrite("n_stations", &ALBPSolution::n_stations,
                           "Number of stations")
            .def_readwrite("cycle_time", &ALBPSolution::cycle_time ,
                            "Cycle time")
            .def_readwrite("loads", &ALBPSolution::loads,
                            "loads of the stations")
            .def_readwrite("n_violations", &ALBPSolution::n_violations,
                           "Number of violations")
            .def_readwrite("n_ranking_violations", &ALBPSolution::n_ranking_violations,
                           "Number of violations from ranking")
            .def_readwrite("method", &ALBPSolution::method)
            .def_readwrite("elapsed_ms", &ALBPSolution::elapsed_ms)

            // Read-only property for n_tasks (since it's private with getter)
            .def_property_readonly("n_tasks", &ALBPSolution::get_n_tasks,
                                   "Number of tasks (read-only)")

            // Public methods
            .def("print", &ALBPSolution::print,
                 "Print the solution")
            .def("task_to_station", &ALBPSolution::task_to_station,
                 "Convert task assignment to station assignment")
            .def("station_to_task", &ALBPSolution::station_to_task,
                 "Convert station assignment to task assignment")
            .def("station_to_ranking", &ALBPSolution::station_to_ranking,
                 "Convert station assignment to ranking")
            .def("ranking_to_task_ranking", &ALBPSolution::ranking_to_task_ranking,
                 "Convert ranking to task ranking")


            // String representation for Python
            .def("__repr__", [](const ALBPSolution &sol) {
                return "<ALBPSolution: " + std::to_string(sol.get_n_tasks()) +
                       " tasks, " + std::to_string(sol.n_stations) +
                       " stations, " + std::to_string(sol.n_violations) + " violations>";
            })

            // Optional: Add a method to get solution summary as dict
            .def("to_dict", [](const ALBPSolution &sol) {
                py::dict d;
                d["n_tasks"] = sol.get_n_tasks();
                d["n_stations"] = sol.n_stations;
                d["n_violations"] = sol.n_violations;
                d["n_ranking_violations"] = sol.n_ranking_violations;
                d["task_assignment"] = sol.task_assignment;
                d["station_assignments"] = sol.station_assignments;
                d["ranking"] = sol.ranking;
                d["task_ranking"] = sol.task_ranking;
                d["solution_time"]= sol.elapsed_ms;
                d["method"] = sol.method;
                return d;
            }, "Convert solution to dictionary");

    m.def("calc_salbp_2_lbs", &calc_salbp_2_lbs,"Calculates several lower bounds for salbp_2 and returns their maximum",
        py::arg("task_time"),
        py::arg("n_stations"));
    m.def("calc_salbp_2_ub", &calc_salbp_2_ub,"Calculates upper bound for salbp_2",
        py::arg("task_time"),
        py::arg("n_stations"));
    // Binding main solver function
    m.def("ils_solve_SALBP1", &ils_solve_SALBP1,
          "Solve SALBP1 using Iterated Local Search",
          py::arg("C"), py::arg("N"), py::arg("task_times"),
          py::arg("raw_precedence"), py::arg("max_iter"),
          py::arg("time_limit"),
          py::arg("op_probs"), py::arg("verbose"),
          py::arg("initial_solution") = std::vector<int>(),
          R"pbdoc(
          Solve SALBP1 using Iterated Local Search

          Parameters:
          -----------
          C : int
              Cycle time
          N : int
              Number of tasks
          task_times : list of int
              Task processing times
          raw_precedence : list of list of int
              Precedence relationships
          max_iter : int
              Maximum iterations
          time_limit : int
            Number of seconds to run
          op_probs : float
              Operation probabilities
          verbose : bool
              Enable verbose output
          initial_solution : list of int, optional
              Initial solution (default: empty)

          Returns:
          --------
          ALBPSolution
              The solved ALBP solution
          )pbdoc");

    m.def("hoff_solve_salbp1",
          [](int C, int N, const std::vector<int> &task_times, const std::vector<std::vector<int> > &raw_precedence) {
              return hoff_solve_salbp1(C, N, task_times, raw_precedence);
          }, py::arg("C"),
          py::arg("N"),
          py::arg("task_times"),
          py::arg("raw_precedence"),
          R"pbdoc(
          Solve SALBP1 using Hoffman heuristic

          Parameters:
          -----------
          C : int
              Cycle time
          N : int
              Number of tasks
          task_times : list of int
              Task processing times
          raw_precedence : list of list of int
              Precedence relationships
          Returns:
          --------
          ALBPSolution
              The solved ALBP solution
          )pbdoc"
    );

    // Second overload using lambda
    m.def("hoff_solve_salbp1",
          [](const ALBP &albp) {
              return hoff_solve_salbp1(albp);
          },
          R"pbdoc(
          Solve SALBP1 using Hoffman heuristic

          Parameters:
          -----------
            albp : ALBP (struct)

          Returns:
          --------
          ALBPSolution
              The solved ALBP solution
          )pbdoc");
    m.def("vdls_solve_salbp1",
          [](int C,
             int N,
             const std::vector<int> &task_times,
             const std::vector<std::vector<int> > &raw_precedence,
             const std::vector<int> &initial_solution,
             const std::optional<int> &max_attempts,
             const std::optional<int> &time_limit
  ) {
              return vdls_solve_salbp1(C, N, task_times, raw_precedence,initial_solution, max_attempts, time_limit);
          }, py::arg("C"),
          py::arg("N"),
          py::arg("task_times"),
          py::arg("raw_precedence"),
          py::arg("initial_solution")= std::vector<int>(),
          py::arg("max_attempts") = std::nullopt,
          py::arg("time_limit") = std::nullopt,
          R"pbdoc(
                          Solve SALBP1 using vdls heuristic

                  Parameters:
                  -----------
                  C : int
                      Cycle time
                  N : int
                      Number of tasks
                  task_times : list of int
                      Task processing times
                  raw_precedence : list of list of int
                      Precedence relationships
                  Returns:
                  --------
                  ALBPSolution
                      The solved ALBP solution
                  )pbdoc"
    );

    m.def("vdls_solve_salbp1",
          [](const ALBP &albp,
             const std::optional<int> &max_attempts,
             const std::optional<int> &time_limit) {
              return vdls_solve_salbp1(albp,  max_attempts, time_limit);
          },
          py::arg("albp"),
          py::arg("max_attempts") = std::nullopt,
          py::arg("time_limit") = std::nullopt,
          R"pbdoc(
              Solve SALBP1 using vdls heuristic

              Parameters:
              -----------
                albp : ALBP (struct)

              Returns:
              --------
              ALBPSolution
                  The solved ALBP solution
              )pbdoc");
    m.def("vdls_solve_salbp2",
          [](int S,
             int N,
             const std::vector<int> &task_times,
             const std::vector<std::vector<int> > &raw_precedence,
             const std::vector<int> &initial_solution,
             const std::optional<int> &max_attempts,
             const std::optional<int> &time_limit
  ) {
              return vdls_solve_salbp2(S, N, task_times, raw_precedence,initial_solution, max_attempts, time_limit);
          }, py::arg("S"),
          py::arg("N"),
          py::arg("task_times"),
          py::arg("raw_precedence"),
          py::arg("initial_solution")= std::vector<int>(),
          py::arg("max_attempts") = std::nullopt,
          py::arg("time_limit") = std::nullopt,
          R"pbdoc(
                          Solve SALBP1 using vdls heuristic

                  Parameters:
                  -----------
                  S : int
                      Number of stations
                  N : int
                      Number of tasks
                  task_times : list of int
                      Task processing times
                  raw_precedence : list of list of int
                      Precedence relationships
                  Returns:
                  --------
                  ALBPSolution
                      The solved ALBP solution
                  )pbdoc"
    );
    m.def("priority_solve_salbp1", &priority_solve_salbp_1,
      "Solve SALBP-1 using priority methods",
      py::arg("C"),
      py::arg("N"),
      py::arg("task_times"),
      py::arg("raw_precedence"),
      py::arg("n_random"),
      R"pbdoc(
                          Solve SALBP1 using different priority methods and station oriented task assignment

                  Parameters:
                  -----------
                  S : int
                      Number of stations
                  N : int
                      Number of tasks
                  task_times : list of int
                      Task processing times
                  raw_precedence : list of list of int
                      Precedence relationships
                  Returns:
                  --------
                  list(ALBPSolution)
                      The solved ALBP solution
                  )pbdoc");
}



