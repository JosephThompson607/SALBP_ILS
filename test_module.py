#!/usr/bin/env python3
import sys
import os



# Add the build directory to Python path
build_dir = 'cmake-build-python_interface/'
sys.path.insert(0, build_dir)
build_dir_2 = 'build/'
sys.path.insert(0, build_dir_2)
print(f"Python version: {sys.version}")
print(f"Python path: {sys.path[:3]}...")  # Show first few paths
print(f"current interpreter: {sys.executable}")
print(f"Looking for module in: {build_dir}")
print(f"Directory exists: {os.path.exists(build_dir)}")

if os.path.exists(build_dir):
    print(f"Available files:")
    for f in sorted(os.listdir(build_dir)):
        if f.endswith('.so') or f.endswith('.pyd'):
            full_path = os.path.join(build_dir, f)
            print(f"  {f} (size: {os.path.getsize(full_path)} bytes)")

sys.path.insert(0, 'cmake-build-python_interface/')
try:


    # Try importing
    import ILS_ALBP
    print("✅ Module imported successfully!")

    # Test the module
    solution = ILS_ALBP.ALBPSolution(5)
    solution.n_stations =  4
    print(f"✅ Created ALBPSolution with {solution.n_tasks} tasks")
    solution.task_assignment=[1,0,3,2,1]
    solution.task_to_station()
    solution.station_to_ranking()
    actual = solution.station_assignments
    expected =  [[1], [0, 4], [3], [2]]
    assert len(actual) == len(expected), f"Length mismatch between station assignments and task assignments: {len(actual)} vs {len(expected)}"
    for i, (actual_sublist, expected_sublist) in enumerate(zip(actual, expected)):
        assert actual_sublist == expected_sublist, f"Mismatch at index between station assignmentts and task assignments{i}: {actual_sublist} vs {expected_sublist}"
    assert len(solution.ranking) ==solution.n_tasks,  f"Length mismatch between ranking and n_tasks: {len(solution.ranking)} vs {len(solution.n_tasks)}"
    correct_ranking = [1, 0, 4, 3, 2]
    for i in range(5):
        assert solution.ranking[i] == correct_ranking[i], "bad ranking"
    print("✅ basic solution manipulation passed!")
except ImportError as e:
    print(f"❌ Import failed: {e}")

    # Try to get more details
    import importlib.util
    spec = importlib.util.find_spec('ILS_ALBP')
    print(f"Module spec: {spec}")

try:


    def ils_call(cycle_time, tasks_times_list, precedence_list,
                 max_iterations=1000, operation_probs=0.5,
                 show_verbose=False, init_sol=None):
        """
        Uses ils to generate a SALBP solution

        Parameters:
        -----------
        cycle_time : int
            Maximum time allowed per workstation
        tasks_times_list : list of int
            Processing time for each task
        precedence_list : list of list of int
            Precedence constraints as [predecessor, successor] pairs
        max_iterations : int, optional
            Maximum iterations for the algorithm (default: 1000)
        operation_probs : float, optional
            Operation probabilities parameter (default: 0.5)
        show_verbose : bool, optional
            Whether to show verbose output (default: False)
        init_sol : list of int, optional
            Initial solution if available (default: None)

        Returns:
        --------
        ALBPSolution or None
            The solution object if successful, None if error occurs
        """

        if init_sol is None:
            init_sol = []

        N = len(tasks_times_list)

        try:
            solution = ILS_ALBP.ils_solve_SALBP1(
                C=cycle_time,
                N=N,
                task_times=tasks_times_list,
                raw_precedence=precedence_list,
                max_iter=max_iterations,
                time_limit = 20,
                op_probs=operation_probs,
                verbose=show_verbose,
                initial_solution=init_sol
            )

            if show_verbose:
                print(f"Successfully solved SALBP1 with {N} tasks and cycle time {cycle_time}")

            return solution

        except Exception as e:
            print(f"Error solving SALBP1: {e}")
            return None
    def mhh_call(cycle_time, tasks_times_list, precedence_list,
      ):


        N = len(tasks_times_list)

        try:
            mhh_sol = ILS_ALBP.hoff_solve_salbp1(
                C=cycle_time,
                N=N,
                task_times=tasks_times_list,
                raw_precedence=precedence_list,

            )

            return mhh_sol

        except Exception as e:
            print(f"Error solving SALBP1: {e}")
            print(f"Here are the function types:C {type(C)}, N {type(N)}, task_times{type(tasks_times_list)} (task_times[0] {type(tasks_times_list[0])}, raw_precedence{type(precedence_list)}")
            return None

    def vdls_call(cycle_time, tasks_times_list, precedence_list,
                     ):


        N = len(tasks_times_list)

        try:
            vdls_sol = ILS_ALBP.vdls_solve_salbp1(
                C=cycle_time,
                N=N,
                task_times=tasks_times_list,
                raw_precedence=precedence_list,
                time_limit = 20,

            )

            return vdls_sol

        except Exception as e:
            print(f"Error solving SALBP1: {e}")
            print(f"Here are the function types:C {type(C)}, N {type(N)}, task_times{type(tasks_times_list)} (task_times[0] {type(tasks_times_list[0])}, raw_precedence{type(precedence_list)}")
            return None

    def vdls_type2_call(S, tasks_times_list, precedence_list,
                      ):


        N = len(tasks_times_list)

        try:
            vdls_sol = ILS_ALBP.vdls_solve_salbp2(
                S=S,
                N=N,
                task_times=tasks_times_list,
                raw_precedence=precedence_list,
                time_limit = 20,

            )

            return vdls_sol

        except Exception as e:
            print(f"Error solving SALBP2: {e}")
            print(f"Here are the function types:S {type(S)}, N {type(N)}, task_times{type(tasks_times_list)} (task_times[0] {type(tasks_times_list[0])}, raw_precedence{type(precedence_list)}")
            return None


    import time
    alb_dict =  {'num_tasks': 50, 'cycle_time': 1000, 'task_times': {'1': 141, '2': 137, '3': 51, '4': 439, '5': 125, '6': 330, '7': 255, '8': 62, '9': 33, '10': 490, '11': 58, '12': 91, '13': 115, '14': 211, '15': 392, '16': 158, '17': 537, '18': 66, '19': 345, '20': 563, '21': 211, '22': 466, '23': 215, '24': 228, '25': 568, '26': 477, '27': 88, '28': 41, '29': 482, '30': 92, '31': 136, '32': 174, '33': 523, '34': 125, '35': 52, '36': 26, '37': 516, '38': 533, '39': 123, '40': 617, '41': 503, '42': 263, '43': 528, '44': 106, '45': 172, '46': 110, '47': 39, '48': 108, '49': 76, '50': 323}, 'precedence_relations': [['1', '4'], ['2', '5'], ['2', '8'], ['2', '9'], ['2', '10'], ['3', '6'], ['3', '7'], ['3', '9'], ['3', '11'], ['4', '12'], ['5', '13'], ['6', '14'], ['8', '16'], ['8', '18'], ['8', '28'], ['9', '15'], ['10', '17'], ['12', '20'], ['13', '21'], ['14', '19'], ['15', '22'], ['18', '23'], ['19', '24'], ['20', '28'], ['21', '26'], ['22', '25'], ['22', '27'], ['22', '33'], ['24', '31'], ['25', '32'], ['26', '29'], ['26', '30'], ['26', '33'], ['27', '34'], ['29', '35'], ['30', '36'], ['31', '39'], ['32', '37'], ['33', '38'], ['33', '40'], ['33', '41'], ['33', '44'], ['34', '42'], ['34', '43'], ['35', '48'], ['36', '48'], ['37', '45'], ['38', '46'], ['39', '48'], ['40', '47'], ['41', '49'], ['42', '50']], 'instance': 'instance_n=50_210'}
    C = alb_dict['cycle_time']
    precs = alb_dict['precedence_relations']
    t_times = [val for _, val in alb_dict['task_times'].items()]
    precs = [[int(child), int(parent)]  for child, parent in alb_dict['precedence_relations']]

    start  = time.time()
    results = ils_call(cycle_time=C, tasks_times_list= t_times, precedence_list=precs,show_verbose=False, max_iterations=1000)
    end = time.time()- start
    print(f"✅ Created ALBPSolution using ils with {results.n_stations} stations in {end} seconds")
    start  = time.time()
    results = mhh_call(cycle_time=C, tasks_times_list= t_times, precedence_list=precs)
    end = time.time()- start
    print(f"✅ Created ALBPSolution using mhh with {results.n_stations} stations in {end} seconds")
    start  = time.time()
    results = vdls_call(cycle_time=C, tasks_times_list= t_times, precedence_list=precs)
    end = time.time()- start
    print(f"✅ Created ALBPSolution using vdls with {results.n_stations} stations in {end} seconds")
    start  = time.time()
    S = 4
    my_lb = ILS_ALBP.calc_salbp_2_lbs(t_times, S)
    my_ub = ILS_ALBP.calc_salbp_2_ub(t_times, S)
    results = vdls_type2_call(S=S, tasks_times_list= t_times, precedence_list=precs)
    end = time.time()- start
    print(f"✅ Created ALBPSolution using vdls with {results.cycle_time} cycle time in {end} seconds. The lower bound is {my_lb}, upperbound is {my_ub}")





except Exception as e:
    print(f"❌ Error testing module: {e}")
    import traceback
    traceback.print_exc()



