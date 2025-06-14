import sys
import os

# Enable verbose import debugging
sys.path.insert(0, 'cmake-build-python_interface/')

# Try importing with more detailed error info
try:
    import ILS_ALBP
    print("✅ Success!")
except Exception as e:
    print(f"❌ Exception type: {type(e).__name__}")
    print(f"❌ Exception message: {e}")

    # Check if it's a symbol/initialization issue
    import traceback
    traceback.print_exc()