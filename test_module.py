#!/usr/bin/env python3
import sys
import os

# Add the build directory to Python path
build_dir = '/Users/letshopethisworks2/CLionProjects/SALBP_ILS/cmake-build-python_interface'
sys.path.insert(0, build_dir)

print(f"Python version: {sys.version}")
print(f"Python path: {sys.path[:3]}...")  # Show first few paths
print(f"Looking for module in: {build_dir}")
print(f"Directory exists: {os.path.exists(build_dir)}")

if os.path.exists(build_dir):
    print(f"Available files:")
    for f in sorted(os.listdir(build_dir)):
        if f.endswith('.so') or f.endswith('.pyd'):
            full_path = os.path.join(build_dir, f)
            print(f"  {f} (size: {os.path.getsize(full_path)} bytes)")

# Test if we can import directly by changing to the directory
original_cwd = os.getcwd()
try:
    os.chdir(build_dir)
    print(f"\nChanged to directory: {os.getcwd()}")

    # Try importing
    import ILS_ALBP
    print("✅ Module imported successfully!")

    # Test the module
    solution = ILS_ALBP.ALBPSolution(5)
    print(f"✅ Created ALBPSolution with {solution.n_tasks} tasks")

except ImportError as e:
    print(f"❌ Import failed: {e}")

    # Try to get more details
    import importlib.util
    spec = importlib.util.find_spec('ILS_ALBP')
    print(f"Module spec: {spec}")

except Exception as e:
    print(f"❌ Error testing module: {e}")
    import traceback
    traceback.print_exc()

finally:
    os.chdir(original_cwd)