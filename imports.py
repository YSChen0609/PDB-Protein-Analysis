# imports.py

"""
The neccessary packages and some self-defined tools (functions)
"""

import time
import requests
import random
import re
from string import Template 
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

# tools
def log_method(func):
    def wrapper(*args, **kwargs):
        print(f"Starting {func.__name__}...")
        start_time = time.time()
        result = func(*args, **kwargs)
        end_time = time.time()
        print(f"Finished {func.__name__}. Time taken: {end_time - start_time} seconds.")
        return result
    return wrapper

def sampleWithConstraints(samplePool, pattern):
    """
    randomly return ONE sample from the samplePool using the constraint string (regex pattern string)
    """
    random.shuffle(samplePool)
    
    # case: match string (attempt)
    for item in samplePool:
        if pattern.match(item):
            return item
    
    # case: no match string
    return False


def normVecCross(x,y,z):
    """
    Get the unit norm of a plane defined by vectors (v x w), where
        - v = y-x
        - w = z-y
        - x, y, z are (1*3)d pandas series
    """
    x = x.to_numpy(dtype=float)
    y = y.to_numpy(dtype=float)
    z = z.to_numpy(dtype=float)
    
    cross_product = np.cross(y-x, z-y)
    
    return cross_product/np.linalg.norm(cross_product)