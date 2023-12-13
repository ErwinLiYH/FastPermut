from .correctness import *
from .load_data import *
from .VEP import *

from .c_conf import *

import time
from functools import wraps

conf = {"timer": False}

def timing_decorator(func):
    @wraps(func)
    def wrapper(*args, **kwargs):
        if conf["timer"]:
            start_time = time.time()  # record the start time
            result = func(*args, **kwargs)  # call the original function
            end_time = time.time()  # record the end time
            print(f"Function {func.__name__} took {end_time - start_time} seconds to run.", flush=True)
            return result
        else:
            return func(*args, **kwargs)
    return wrapper