#!/usr/bin/env python3
"""
Test script for the SALBP1_heuristics compiled extension.

Organized into three parts:
    1. Setup            - locate/import the compiled module, solver call wrappers
    2. Unit tests        - exercise individual ALBP / ALBPSolution behaviors in isolation
    3. Regression tests  - run each heuristic solver end-to-end on a benchmark instance

Run directly:
    python test_salbp_heuristics.py
"""
import importlib.util
import os
import sys
import time
from copy import deepcopy


# ===========================================================================
# 1. Setup: locate/import the compiled module, and thin wrappers around
#    each solver's keyword-argument plumbing.
# ===========================================================================

def import_salbp():
    """Locate and import the compiled SALBP1_heuristics extension.

    Returns the module, or None if it couldn't be imported (details are
    printed either way).
    """
    for build_dir in ("cmake-build-python_interface/", "build/"):
        if build_dir not in sys.path:
            sys.path.insert(0, build_dir)

    print(f"Python version: {sys.version}")
    print(f"Python path: {sys.path[:3]}...")
    print(f"current interpreter: {sys.executable}")

    build_dir = "cmake-build-python_interface/"
    print(f"Looking for module in: {build_dir}")
    print(f"Directory exists: {os.path.exists(build_dir)}")
    if os.path.exists(build_dir):
        print("Available files:")
        for f in sorted(os.listdir(build_dir)):
            if f.endswith(".so") or f.endswith(".pyd"):
                full_path = os.path.join(build_dir, f)
                print(f"  {f} (size: {os.path.getsize(full_path)} bytes)")

    try:
        import SALBP1_heuristics
        print("✅ Module imported successfully!")
        print(dir(SALBP1_heuristics))
        return SALBP1_heuristics
    except ImportError as e:
        print(f"❌ Import failed: {e}")
        spec = importlib.util.find_spec("SALBP1_heuristics")
        print(f"Module spec: {spec}")
        return None


def ils_call(salbp, cycle_time, task_times_list, precedence_list,
             max_iterations=1000, operation_probs=0.5, show_verbose=False, init_sol=None):
    """Solve SALBP-1 with Iterated Local Search. Returns None on error."""
    if init_sol is None:
        init_sol = []
    try:
        solution = salbp.ils_solve_SALBP1(
            C=cycle_time,
            N=len(task_times_list),
            task_times=task_times_list,
            raw_precedence=precedence_list,
            max_iter=max_iterations,
            time_limit=20,
            op_probs=operation_probs,
            verbose=show_verbose,
            initial_solution=init_sol,
        )
        if show_verbose:
            print(f"Successfully solved SALBP1 with {len(task_times_list)} tasks and cycle time {cycle_time}")
        return solution
    except Exception as e:
        print(f"Error solving SALBP1 (ils): {e}")
        return None


def hoff_call(salbp, cycle_time, task_times_list, precedence_list):
    """Solve SALBP-1 with the Hoffmann heuristic. Returns None on error."""
    try:
        return salbp.hoff_solve_salbp1(
            C=cycle_time,
            N=len(task_times_list),
            task_times=task_times_list,
            raw_precedence=precedence_list,
            alpha_iter=2,
            alpha_size=0.005,
            beta_iter=int(len(task_times_list) / 2),
            beta_size=0.005,
            reverse=True,
        )
    except Exception as e:
        print(f"Error solving SALBP1 (hoff): {e}")
        return None


def mhh_call(salbp, cycle_time, task_times_list, precedence_list):
    """Solve SALBP-1 with MHH. Returns None on error."""
    try:
        return salbp.mhh_solve_salbp1(
            C=cycle_time,
            N=len(task_times_list),
            task_times=task_times_list,
            raw_precedence=precedence_list,
            alpha_schedule = [0.0,0.5],
            beta_schedule = [0, 0.5]
        )
    except Exception as e:
        print(f"Error solving SALBP1 (mhh): {e}")
        return None

def poke_mhh(salbp, cycle_time, task_times_list, precedence_list, alpha=None, beta=None, ranking=None, gamma=None, seed=None):
    """Solve SALBP-1 with MHH, using alpha and beta schedules. Returns None on error."""
    try:
        sol_1 = salbp.mhh_solve_salbp1(
            C=cycle_time,
            N=len(task_times_list),
            task_times=task_times_list,
            raw_precedence=precedence_list,
            alpha_schedule = alpha,
            beta_schedule = beta,
            gamma = gamma,
            task_priorities = ranking,
            seed = seed,
        )


        return sol_1

    except Exception as e:
        print(f"Error solving SALBP1 (mhh): {e}")
        return None

def vdls_call(salbp, cycle_time, task_times_list, precedence_list):
    """Solve SALBP-1 with VDLS. Returns None on error."""
    try:
        return salbp.vdls_solve_salbp1(
            C=cycle_time,
            N=len(task_times_list),
            task_times=task_times_list,
            raw_precedence=precedence_list,
            time_limit=1.5,
        )
    except Exception as e:
        print(f"Error solving SALBP1 (vdls): {e}")
        return None


def vdls_type2_call(salbp, S, task_times_list, precedence_list):
    """Solve SALBP-2 with VDLS. Returns None on error."""
    try:
        return salbp.vdls_solve_salbp2(
            S=S,
            N=len(task_times_list),
            task_times=task_times_list,
            raw_precedence=precedence_list,
            time_limit=20,
        )
    except Exception as e:
        print(f"Error solving SALBP2 (vdls): {e}")
        return None


def priority_type1_call(salbp, C, task_times_list, precedence_list, n_random=3, time_limit=1):
    """Solve SALBP-1 with priority-rule heuristics. Returns None on error."""
    try:
        return salbp.priority_solve_salbp1(
            C=C,
            N=len(task_times_list),
            task_times=task_times_list,
            raw_precedence=precedence_list,
            n_random=n_random,
            seed=42,
            time_limit=time_limit,
        )
    except Exception as e:
        print(f"Error solving SALBP1 (priority): {e}")
        return None


def priority_type2_call(salbp, S, task_times_list, precedence_list, move_target=False):
    """Solve SALBP-2 with priority-rule heuristics. Returns None on error."""
    try:
        return salbp.priority_solve_salbp2(
            S=S,
            N=len(task_times_list),
            task_times=task_times_list,
            raw_precedence=precedence_list,
            n_random=3,
            move_target=move_target,
            seed=42,
        )
    except Exception as e:
        print(f"Error solving SALBP2 (priority): {e}")
        return None


def tabu_call(salbp, cycle_time, task_times_list, precedence_list, time_limit=1.0):
    """Solve SALBP-1 with Tabu search. Returns None on error."""
    try:
        return salbp.tabu_solve_salbp1(
            C=cycle_time,
            N=len(task_times_list),
            task_times=task_times_list,
            raw_precedence=precedence_list,
            time_limit=time_limit,
        )
    except Exception as e:
        print(f"Error solving SALBP1 (tabu): {e}")
        return None


# ===========================================================================
# 2. Unit tests: individual ALBP / ALBPSolution behaviors, in isolation.
# ===========================================================================

def test_solution_manipulation(salbp):
    """ALBPSolution: task_to_station, station_to_ranking, reverse."""
    solution = salbp.ALBPSolution(5)
    solution.n_stations = 4
    print(f"✅ Created ALBPSolution with {solution.n_tasks} tasks")

    print("Testing solution manipulation")

    solution.task_assignment = [1, 0, 3, 2, 1]
    solution.task_to_station()
    solution.station_to_ranking()
    actual = solution.station_assignments
    expected = [[1], [0, 4], [3], [2]]
    assert len(actual) == len(expected), (
        f"Length mismatch between station assignments and task assignments: {len(actual)} vs {len(expected)}"
    )
    for i, (actual_sublist, expected_sublist) in enumerate(zip(actual, expected)):
        assert actual_sublist == expected_sublist, (
            f"Mismatch at index {i} between station assignments and task assignments: {actual_sublist} vs {expected_sublist}"
        )

    assert len(solution.ranking) == solution.n_tasks, (
        f"Length mismatch between ranking and n_tasks: {len(solution.ranking)} vs {solution.n_tasks}"
    )
    correct_ranking = [1, 0, 4, 3, 2]
    for i in range(5):
        assert solution.ranking[i] == correct_ranking[i], "bad ranking"

    # Reverse functionality
    #[1, 0, 3, 2, 1]
    solution.reverse()
    assert solution.task_assignment == [2, 3, 0, 1, 2], f"Solution not reversed correctly {solution.task_assignment}"
    actual = solution.station_assignments
    expected.reverse()
    assert len(actual) == len(expected), (
        f"Length mismatch between station assignments and task assignments: {len(actual)} vs {len(expected)}"
    )
    for i, (actual_sublist, expected_sublist) in enumerate(zip(actual, expected)):
        assert actual_sublist == expected_sublist, (
            f"Mismatch at index {i} between station assignments and task assignments: {actual_sublist} vs {expected_sublist}"
        )
    solution.loads= [10,20,30,23, 45]
    solution.reverse()
    assert solution.loads[1] == 23, 'load in wrong spot'

    print("✅ basic solution manipulation passed!")


def test_lower_bounds(salbp):
    """calc_salbp_1_lb6 sanity check."""
    task_times = [8, 8, 7, 6, 5, 5, 3, 4, 3]
    C = 10
    lb_6 = salbp.calc_salbp_1_lb6(task_times, C)
    assert lb_6 == 6, f"LB 6 is not acting as expected. We wanted 6, we got {lb_6}"
    print("✅ lb tests passed!")


def test_reverse_precedence_matrix(salbp, C, task_times_list, precedence_list):
    """The 'reverse' flag on ALBP.type_1 should produce the transpose of the precedence matrix."""
    n = len(task_times_list)
    albp = salbp.ALBP.type_1(C, n, task_times_list, precedence_list, False)
    albp_rev = salbp.ALBP.type_1(C, n, task_times_list, precedence_list, True)
    for i in range(n):
        for j in range(n):
            assert albp.prec_mat[i * n + j] == albp_rev.prec_mat[j * n + i], (
                f"precedence matrix is not the transpose, see {i},{j}"
            )
    print("✅ Tested reverse function of ALBP")


def test_add_precedence_relation(salbp):
    """Adding a precedence relation should update both the adjacency matrix
    and the transitive closure matrix."""
    task_times = [12, 7, 15, 9, 18]
    N = len(task_times)
    # Precedence relations [parent, child], one-indexed: 1 -> 3 -> 5, 2 -> 4
    raw_precedence = [[1, 3], [2, 4], [3, 5]]

    albp = salbp.ALBP.type_1(C=1000, N=N, task_times=task_times, raw_precedence=raw_precedence)
    albp.add_precedence_relation([3, 2])

    assert albp.t_close_mat[0 + 1] == 1, "1 did not get 2 as a successor in transitive closure"
    assert albp.t_close_mat[0 + 3] == 1, "1 did not get 4 as a successor in transitive closure"
    assert albp.t_close_mat[2 * 5 + 3] == 1, "3 did not get 4 as a successor in transitive closure"
    assert albp.prec_mat[2 * 5 + 1] == 1, "3 did not get 2 as a successor is adjacency matrix"
    print("✅ tested precedence constraint insertion")

def test_get_critical_paths(salbp):
    """Depth/height should reflect the longest path (in edges) from any source
    to each node and from each node to any sink, with depth_pred/height_suc
    pointing to the neighbor that produced that longest path."""
    task_times = [12, 7, 15, 9, 18]
    N = len(task_times)
    # Precedence relations [parent, child], one-indexed: 1 -> 3 -> 5, 2 -> 4
    raw_precedence = [[1, 3], [2, 4], [3, 5]]

    albp = salbp.ALBP.type_1(C=1000, N=N, task_times=task_times, raw_precedence=raw_precedence)
    cp = albp.get_critical_paths()

    # depth[i] = # edges on the longest path from any source down to task i+1
    # Task1, Task2: sources -> depth 0
    # Task3: via Task1 -> depth 1
    # Task4: via Task2 -> depth 1
    # Task5: via Task1->Task3 -> depth 2
    assert cp.depth == [0, 0, 1, 1, 2], "Depth values incorrect"

    # height[i] = # edges on the longest path from task i+1 down to any sink
    # Task1: ->Task3->Task5 -> height 2
    # Task2: ->Task4 -> height 1
    # Task3: ->Task5 -> height 1
    # Task4, Task5: sinks -> height 0
    assert cp.height == [2, 1, 1, 0, 0], "Height values incorrect"

    # depth_pred[i]: predecessor that gave task i+1 its depth (-1 = source)
    assert cp.depth_pred == [-1, -1, 0, 1, 2], "Depth predecessor pointers incorrect"

    # height_suc[i]: successor that gave task i+1 its height (-1 = sink)
    assert cp.height_suc == [2, 3, 4, -1, -1], "Height successor pointers incorrect"

    print("✅ tested critical path (depth/height) calculation")

def test_ranking_solution(salbp):
    """Depth/height should reflect the longest path (in edges) from any source
    to each node and from each node to any sink, with depth_pred/height_suc
    pointing to the neighbor that produced that longest path."""
    task_times = [12, 7, 15, 9, 18]
    N = len(task_times)
    C = 18
    # Precedence relations [parent, child], one-indexed: 1 -> 3 -> 5, 2 -> 4
    raw_precedence = [[1, 3], [2, 4], [3, 5]]
    ranking = [3, 1, 0, 4, 2]
    sol = salbp.ranking_to_solution(C, N, task_times, raw_precedence, ranking )
    # depth[i] = # edges on the longest path from any source down to task i+1
    # Task1, Task2: sources -> depth 0
    # Task3: via Task1 -> depth 1
    # Task4: via Task2 -> depth 1
    # Task5: via Task1->Task3 -> depth 2
    assert sol.task_assignment == [1,0,2,0,3], "unexpected assignment"
    assert sol.ranking == [1,3,0,2,4], f"ranking is off: {ranking} "
    assert max(sol.loads) <= C, "Cycle time exceeded"
    # height[i] = # edges on the longest path from task i+1 down to any sink
    # Task1: ->Task3->Task5 -> height 2
    # Task2: ->Task4 -> height 1
    # Task3: ->Task5 -> height 1
    # Task4, Task5: sinks -> height 0


    print("✅ tested ranking solutions")
def test_get_path_stats(salbp):
    """Check summary statistics for a subset of tasks."""
    task_times = [12, 7, 15, 9, 18]
    N = len(task_times)
    raw_precedence = [[1, 3], [2, 4], [3, 5]]

    albp = salbp.ALBP.type_1(
        C=1000,
        N=N,
        task_times=task_times,
        raw_precedence=raw_precedence
    )

    # Tasks 1, 3, 5 (0-indexed: 0, 2, 4)
    stats = albp.get_path_stats([0, 2, 4])

    assert stats.total == 45
    assert stats.total_sq == 693
    assert stats.min == 12
    assert stats.max == 18
    assert stats.mean == 15
    assert stats.variance == 6

    print("✅ tested path statistics calculation")

def test_get_positional_weight(salbp):
    """Positional weight should equal each task's own time plus the times
    of all its transitive successors."""
    task_times = [12, 7, 15, 9, 18]
    N = len(task_times)
    # Precedence relations [parent, child], one-indexed: 1 -> 3 -> 5, 2 -> 4
    raw_precedence = [[1, 3], [2, 4], [3, 5]]

    albp = salbp.ALBP.type_1(C=1000, N=N, task_times=task_times, raw_precedence=raw_precedence)
    weights = salbp.get_positional_weight(albp)

    # Task 1: own(12) + successors 3(15) and 5(18)
    assert weights[0] == 45, "Task 1 positional weight incorrect"
    # Task 2: own(7) + successor 4(9)
    assert weights[1] == 16, "Task 2 positional weight incorrect"
    # Task 3: own(15) + successor 5(18)
    assert weights[2] == 33, "Task 3 positional weight incorrect"
    # Task 4: own(9), no successors
    assert weights[3] == 9, "Task 4 positional weight incorrect"
    # Task 5: own(18), no successors
    assert weights[4] == 18, "Task 5 positional weight incorrect"
    print("✅ tested positional weight calculation")


def test_get_reverse_positional_weight(salbp):
    """Reverse positional weight should equal each task's own time plus the
    times of all its transitive predecessors."""
    task_times = [12, 7, 15, 9, 18]
    N = len(task_times)
    # Precedence relations [parent, child], one-indexed: 1 -> 3 -> 5, 2 -> 4
    raw_precedence = [[1, 3], [2, 4], [3, 5]]

    albp = salbp.ALBP.type_1(C=1000, N=N, task_times=task_times, raw_precedence=raw_precedence)
    weights = salbp.get_reverse_positional_weight(albp)

    # Task 1: own(12), no predecessors
    assert weights[0] == 12, "Task 1 reverse positional weight incorrect"
    # Task 2: own(7), no predecessors
    assert weights[1] == 7, "Task 2 reverse positional weight incorrect"
    # Task 3: own(15) + predecessor 1(12)
    assert weights[2] == 27, "Task 3 reverse positional weight incorrect"
    # Task 4: own(9) + predecessor 2(7)
    assert weights[3] == 16, "Task 4 reverse positional weight incorrect"
    # Task 5: own(18) + predecessors 1(12) and 3(15) (transitive)
    assert weights[4] == 45, "Task 5 reverse positional weight incorrect"
    print("✅ tested reverse positional weight calculation")

def test_deepcopy(salbp):
    """deepcopy of an ALBP should be independent of the original."""
    task_times = [12, 7, 15, 9, 18]
    N = len(task_times)
    raw_precedence = [[1, 3], [2, 4], [3, 5]]

    albp = salbp.ALBP.type_1(C=1000, N=N, task_times=task_times, raw_precedence=raw_precedence)
    albp2 = deepcopy(albp)
    assert len(albp.task_time) == len(albp2.task_time), "Tasks do not line up with copy"

    albp.add_precedence_relation([3, 2])
    assert len(albp.precedence_relations) == len(albp2.precedence_relations) + 1, (
        f"new number of precedence relations appear to be incorrect: "
        f"{len(albp.precedence_relations)} vs {len(albp2.precedence_relations)}"
    )
    assert len(albp.dir_suc[2]) > len(albp2.dir_suc[2]), (
        "Wrong number of precedence constraints for copy after adding edge to original"
    )
    print("✅ tested copying")


def test_heads_and_tails(salbp, cycle_time, task_times_list, precedence_list):
    """get_heads / get_tails / get_topo_sort sanity checks."""
    albp = salbp.ALBP.type_1(cycle_time, len(task_times_list), task_times_list, precedence_list, False, False, False)
    heads = salbp.get_heads(albp, False)
    tails = salbp.get_tails(albp, False)
    topo_sort = salbp.get_topo_sort(albp.dir_pred, albp.dir_suc)

    assert sum(albp.task_time) == albp.total_time, (
        f"task times {sum(albp.task_time)} should equal total time {albp.total_time}"
    )
    assert len(topo_sort) == len(task_times_list), "length of topological sort does not match"
    assert len(heads) == len(task_times_list), "number of heads does not match the number of tasks"
    assert len(tails) == len(task_times_list), "number of tails does not match the number of tasks"
    print("✅ Tested heads and tails function of ALBP")


# ===========================================================================
# 3. Regression tests: run each solver end-to-end on the benchmark instance.
# ===========================================================================

def test_ils(salbp, C, t_times, precs):
    start = time.time()
    results = ils_call(salbp, cycle_time=C, task_times_list=t_times, precedence_list=precs, max_iterations=1000)
    assert sum(results.loads) == sum(t_times), 'Task times mismatched with loads'
    assert max(results.loads) <= C, 'Max load greater than cycle time'
    print(f"✅ Created ALBPSolution using ils with {results.n_stations} stations in {time.time() - start} seconds")


def test_tabu(salbp, C, t_times, precs):
    start = time.time()
    results = tabu_call(salbp, cycle_time=C, task_times_list=t_times, precedence_list=precs, time_limit=3)
    assert sum(results.loads) == sum(t_times), 'Task times mismatched with loads'
    assert max(results.loads) <= C, 'Max load greater than cycle time'
    print(f"✅ Created ALBPSolution using tabu with {results.n_stations} stations in {time.time() - start} seconds")


def test_priority_type1(salbp, C, t_times, precs):
    start = time.time()
    results = priority_type1_call(salbp, C, task_times_list=t_times, precedence_list=precs)
    print(f"✅ Created ALBPSolution using priority in {time.time() - start} seconds")


def test_priority_type1_timed(salbp, C, t_times, precs):
    start = time.time()
    results = priority_type1_call(salbp, C, task_times_list=t_times, precedence_list=precs,
                                  n_random=1_000_000, time_limit=0.90)
    for result in results:
        assert sum(result.loads) == sum(t_times), 'Task times mismatched with loads'
        assert max(result.loads) <= C, 'Max load greater than cycle time'
    print(f"✅ Created timed ALBPSolution using priority in {time.time() - start} seconds "
          f"and {len(results)} results")


def test_priority_type2(salbp, t_times, precs):
    start = time.time()
    results = priority_type2_call(salbp, 20, task_times_list=t_times, precedence_list=precs, move_target=True)
    print(f"✅ Created SALBP2 solutions using priority in {time.time() - start} seconds")
    # for result in results:
    #     print(f"    here is the method {result.method} with the cycle time {result.cycle_time}")


def test_hoff(salbp, C, t_times, precs):
    start = time.time()
    results = hoff_call(salbp, cycle_time=C, task_times_list=t_times, precedence_list=precs)
    assert sum(results.loads) == sum(t_times), 'Task times mismatched with loads'
    assert max(results.loads) <= C, 'Max load greater than cycle time'
    print(f"✅ Created ALBPSolution using hoff with {results.n_stations} stations in {time.time() - start} seconds")


def test_mhh(salbp, C, t_times, precs):
    start = time.time()
    results = mhh_call(salbp, cycle_time=C, task_times_list=t_times, precedence_list=precs)
    assert sum(results.loads) == sum(t_times), 'Task times mismatched with loads'
    assert max(results.loads) <= C, 'Max load greater than cycle time'
    print(f"✅ Created ALBPSolution using mhh with {results.n_stations} stations in {time.time() - start} seconds")
    # print("here are the station loads", results.loads)

def test_alpha_beta_mhh(salbp, C, t_times, precs):
    start = time.time()
    results1 = poke_mhh(salbp, cycle_time=C, task_times_list=t_times, precedence_list=precs)
    results2 = poke_mhh(salbp, cycle_time=C, task_times_list=t_times, precedence_list=precs, alpha = [0.2], beta=[0.2])
    results3 = poke_mhh(salbp, cycle_time=C, task_times_list=t_times, precedence_list=precs, alpha = [0], beta=[1.0])
    assert sum(results1.loads) == sum(results2.loads), "sum of loads does not match"
    print(f"✅ Created ALBPSolution using mhh (ALPHA BETA test) with {results1.n_stations} , {results2.n_stations} , {results3.n_stations}  stations in {time.time() - start} seconds")
    # print(f"here are the station assignments {results1.task_assignment} , {results2.task_assignment} , {results3.task_assignment}")
    # print(f"here are the station loads {results1.loads} , {results2.loads} , {results3.loads}")

def test_gamma_mhh(salbp, C, t_times, precs):
    start = time.time()
    results1 = poke_mhh(salbp, cycle_time=C, task_times_list=t_times, precedence_list=precs,gamma=0.5, seed=42)
    results2 = poke_mhh(salbp, cycle_time=C, task_times_list=t_times, precedence_list=precs, gamma=0.5, seed=42)
    results3 = poke_mhh(salbp, cycle_time=C, task_times_list=t_times, precedence_list=precs, gamma = 0.01, seed=42)
    assert sum(results1.loads) == sum(t_times), 'Task times mismatched with loads'
    assert max(results1.loads) <= C, 'Max load greater than cycle time'
    assert sum(results1.loads) == sum(results2.loads) ==sum(results3.loads) , "sum of loads does not match"
    assert results1.task_assignment == results2.task_assignment, "Seeding seems to not control the rng"
    assert results1.task_assignment != results3.task_assignment, "perturbations not working"
    print(f"✅ Created ALBPSolution using mhh (GAMMA TEST) with {results1.n_stations} , {results2.n_stations} , {results3.n_stations}  stations in {time.time() - start} seconds")
    # print(f"here are the station assignments {results1.task_assignment} , {results2.task_assignment} , {results3.task_assignment}")
    # print(f"here are the station loads {results1.loads} , {results2.loads} , {results3.loads}")

def test_priority_change_mhh(salbp, C, t_times, precs):
    start = time.time()
    results1 = poke_mhh(salbp, cycle_time=C, task_times_list=t_times, precedence_list=precs)
    results2 = poke_mhh(salbp, cycle_time=C, task_times_list=t_times, precedence_list=precs, ranking = list(range(len(t_times), -1, -1)))
    assert sum(results1.loads) == sum(results2.loads), "sum of loads does not match"
    assert sum(results1.loads) == sum(t_times), 'Task times mismatched with loads'
    assert max(results1.loads) <= C, 'Max load greater than cycle time'
    assert results2.task_assignment != results1.task_assignment, 'task time change seems to fail'
    print(f"✅ Created ALBPSolution using mhh with {results1.n_stations} , {results2.n_stations} ,  stations in {time.time() - start} seconds")
    # print(f"here are the station assignments {results1.task_assignment} , {results2.task_assignment} , ")
    # print(f"here are the station loads {results1.loads} , {results2.loads} ,")


def test_vdls(salbp, C, t_times, precs):
    start = time.time()
    results = vdls_call(salbp, cycle_time=C, task_times_list=t_times, precedence_list=precs)
    assert sum(results.loads) == sum(t_times), 'Task times mismatched with loads'
    assert max(results.loads) <= C, 'Max load greater than cycle time'
    print(f"✅ Created ALBPSolution using vdls with {results.n_stations} stations in {time.time() - start} seconds")
    print("here are the station loads", results.loads)


def test_vdls_type2(salbp, t_times, precs):
    start = time.time()
    S = 4
    my_lb = salbp.calc_salbp_2_lbs(t_times, S)
    my_ub = salbp.calc_salbp_2_ub(t_times, S)
    results = vdls_type2_call(salbp, S=S, task_times_list=t_times, precedence_list=precs)
    assert sum(results.loads) == sum(t_times), 'Task times mismatched with loads'
    print(f"✅ Created ALBPSolution using vdls with {results.cycle_time} cycle time in {time.time() - start} "
          f"seconds. The lower bound is {my_lb}, upperbound is {my_ub}")
    print("here are the station loads", results.loads)
    print("Testing results to dict", results.to_dict())


# ===========================================================================
# Benchmark instance used by the regression tests
# ===========================================================================

def load_n50_instance():
    """50-task benchmark instance ('instance_n=50_210')."""
    alb_dict = {
        'num_tasks': 50, 'cycle_time': 1000,
        'task_times': {'1': 141, '2': 137, '3': 51, '4': 439, '5': 125, '6': 330, '7': 255, '8': 62, '9': 33,
                       '10': 490, '11': 58, '12': 91, '13': 115, '14': 211, '15': 392, '16': 158, '17': 537,
                       '18': 66, '19': 345, '20': 563, '21': 211, '22': 466, '23': 215, '24': 228, '25': 568,
                       '26': 477, '27': 88, '28': 41, '29': 482, '30': 92, '31': 136, '32': 174, '33': 523,
                       '34': 125, '35': 52, '36': 26, '37': 516, '38': 533, '39': 123, '40': 617, '41': 503,
                       '42': 263, '43': 528, '44': 106, '45': 172, '46': 110, '47': 39, '48': 108, '49': 76,
                       '50': 323},
        'precedence_relations': [['1', '4'], ['2', '5'], ['2', '8'], ['2', '9'], ['2', '10'], ['3', '6'],
                                 ['3', '7'], ['3', '9'], ['3', '11'], ['4', '12'], ['5', '13'], ['6', '14'],
                                 ['8', '16'], ['8', '18'], ['8', '28'], ['9', '15'], ['10', '17'], ['12', '20'],
                                 ['13', '21'], ['14', '19'], ['15', '22'], ['18', '23'], ['19', '24'],
                                 ['20', '28'], ['21', '26'], ['22', '25'], ['22', '27'], ['22', '33'],
                                 ['24', '31'], ['25', '32'], ['26', '29'], ['26', '30'], ['26', '33'],
                                 ['27', '34'], ['29', '35'], ['30', '36'], ['31', '39'], ['32', '37'],
                                 ['33', '38'], ['33', '40'], ['33', '41'], ['33', '44'], ['34', '42'],
                                 ['34', '43'], ['35', '48'], ['36', '48'], ['37', '45'], ['38', '46'],
                                 ['39', '48'], ['40', '47'], ['41', '49'], ['42', '50']],
        'instance': 'instance_n=50_210',
    }
    C = alb_dict['cycle_time']
    t_times = [val for _, val in alb_dict['task_times'].items()]
    precs = [[int(child), int(parent)] for child, parent in alb_dict['precedence_relations']]
    return C, t_times, precs


# ===========================================================================
# Runner: executes every test, isolating failures so one broken test
# doesn't prevent the rest of the suite from reporting results.
# ===========================================================================

def run_test(name, func):
    try:
        func()
        return True
    except Exception as e:
        print(f"❌ {name} failed: {e}")
        import traceback
        traceback.print_exc()
        return False


def main():
    salbp = import_salbp()
    if salbp is None:
        print("❌ Cannot run tests: SALBP1_heuristics module is not available")
        return

    C, t_times, precs = load_n50_instance()

    unit_tests = [
        ("solution_manipulation", lambda: test_solution_manipulation(salbp)),
        ("lower_bounds", lambda: test_lower_bounds(salbp)),
        ("reverse_precedence_matrix", lambda: test_reverse_precedence_matrix(salbp, C, t_times, precs)),
        ("add_precedence_relation", lambda: test_add_precedence_relation(salbp)),
        ("deepcopy", lambda: test_deepcopy(salbp)),
        ("heads_and_tails", lambda: test_heads_and_tails(salbp, C, t_times, precs)),
        ("positional weight", lambda: test_get_positional_weight(salbp)),
        ("test reverse positional weight", lambda: test_get_reverse_positional_weight(salbp)),
        ("critical paths", lambda: test_get_critical_paths(salbp)),
        ("test_ranking_solution", lambda: test_ranking_solution(salbp)),
        ("path stats", lambda: test_get_path_stats(salbp)),

    ]

    regression_tests = [
        ("ils", lambda: test_ils(salbp, C, t_times, precs)),
        ("tabu", lambda: test_tabu(salbp, C, t_times, precs)),
        ("priority_type1", lambda: test_priority_type1(salbp, C, t_times, precs)),
        ("priority_type1_timed", lambda: test_priority_type1_timed(salbp, C, t_times, precs)),
        ("priority_type2", lambda: test_priority_type2(salbp, t_times, precs)),
        ("hoff", lambda: test_hoff(salbp, C, t_times, precs)),
        ("mhh", lambda: test_mhh(salbp, C, t_times, precs)),
        ("alpha beta mhh", lambda: test_alpha_beta_mhh(salbp, C, t_times, precs)),
        ("gamma mhh", lambda: test_gamma_mhh(salbp, C, t_times, precs)),
        ("priority change mhh", lambda: test_priority_change_mhh(salbp, C, t_times, precs)),
        ("vdls", lambda: test_vdls(salbp, C, t_times, precs)),
        ("vdls_type2", lambda: test_vdls_type2(salbp, t_times, precs)),

    ]

    print("\n=== Unit tests ===")
    unit_results = [(name, run_test(name, func)) for name, func in unit_tests]

    print("\n=== Regression tests ===")
    regression_results = [(name, run_test(name, func)) for name, func in regression_tests]

    all_results = unit_results + regression_results
    passed = sum(1 for _, ok in all_results if ok)
    print(f"\n=== Summary: {passed}/{len(all_results)} tests passed ===")
    for name, ok in all_results:
        print(f"  {'✅' if ok else '❌'} {name}")


if __name__ == "__main__":
    main()