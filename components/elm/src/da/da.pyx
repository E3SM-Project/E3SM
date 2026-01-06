# da.pyx
# cython: language_level=3
import pandas as pd
import numpy as np
import os
from datetime import datetime
import time

# Global variables for debugging
cdef dict debug_data = {}         # Stores all inputs for debugging
cdef object debug_df = None       # DataFrame for saving inputs
cdef int debug_counter = 0

cdef public void python_init():
    global debug_df
    import time
    time1 = time.time()
    import numpy as np
    
    print("Initialization started in python_init")

